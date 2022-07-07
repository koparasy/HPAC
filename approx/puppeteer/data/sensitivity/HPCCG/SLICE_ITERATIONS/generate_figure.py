#!/g/g15/vanover1/GSA/.venv/bin/python

import sys
import os
import json
import warnings
import subprocess

import numpy as np

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio; pio.kaleido.scope.mathjax = None


def visualize(all_data, outfile_base_name, fn):

    kernel_labels = ["DOT_r", "DOT_\u03B1", "MATVEC", "WXPY_p", "WXPY_x", "WXPY_r"]
    section_labels = ["1", "2", "3"]
    method_names = ['SALib-Morris', 'SALib-Sobol', 'IVARS', 'VARS-Morris-ABE', 'VARS-Sobol-TO']
    method_n = [v.replace("-TO", "").replace("-ABE", "") for v in method_names]

    fig = make_subplots(rows = 1, cols = 5, shared_yaxes=True, horizontal_spacing=0.01)

    for trace_i  in range(len(all_data.items())):
        method_name = method_names[trace_i]
        data_dict = all_data[method_name]
        method_name = method_name.replace("-TO", "").replace("-ABE", "")

        # gather all section names and kernel names
        section_names = set()
        kernel_names = set()
        for name in data_dict["names"]:
            section_name = " ".join(name.split("_")[:2])
            kernel_name = "_".join(name.split("_")[2:])
            section_names.add(section_name)
            kernel_names.add(kernel_name)
        section_names = sorted(section_names)
        kernel_names = sorted(kernel_names)

        # gather all sensitivities and store them based on section
        # name and kernel name
        sensitivities = np.zeros((len(section_names), len(kernel_names)))
        for i, sensitivity in enumerate(data_dict["Sensitivity"]):
            name = data_dict["names"][i]
            section_name = " ".join(name.split("_")[:2])
            kernel_name = "_".join(name.split("_")[2:])
            sensitivities[section_names.index(section_name)][kernel_names.index(kernel_name)] = sensitivity

        # give sensitivities a logarithmic scale
        warnings.filterwarnings("ignore", "divide by zero encountered in log")
        sensitivities_log = np.log(sensitivities)
        warnings.resetwarnings()

        # if any of the sensitivities were 0, impute the -inf values
        # given by np.log(0) with twice the magnitude of the log of the
        # smallest nonzero sensitivity
        if not sensitivities.all():
            sorted_flattened_sensitivities = sorted(sensitivities.flatten(), reverse=True)
            second_smallest = sorted_flattened_sensitivities.pop()
            while second_smallest == 0:
                second_smallest = sorted_flattened_sensitivities.pop()
            sensitivities_log[np.where(sensitivities == 0)] = np.log(second_smallest) * 2
            zero_flag = True
        else:
            zero_flag = False

        # normalize sensitivities by scaling them so the min is 0 and
        # then dividing by the maximum value
        sensitivities_log_norm = sensitivities_log + np.abs(np.min(sensitivities_log))
        sensitivities_log_norm = sensitivities_log_norm / np.max(sensitivities_log_norm)

        # plot each heatmap
        fig.add_trace(
            go.Heatmap(
                x = section_labels,
                y = kernel_labels,
                z = sensitivities_log_norm.T,
                name = method_name,
                colorbar = dict(
                    tickmode = "linear",
                    tick0 = 0,
                    dtick = 1.0,
                ),
            ),
            row = 1,
            col = trace_i + 1
        )

    # for subplot titles
    fig.update_annotations(font_size=40)

    fig.add_annotation(
        text = 'Sensitivity',
        x = 1.1,
        y = 0.4,
        xref = 'paper',
        yref = 'paper',
        textangle = 90,
        font_size = 40,
    )

    x_start = -0.20
    for method_name in method_n:
        x_start += 0.2

        if method_name == "IVARS":
            xx = x_start -0.05
        else:
            xx = x_start

        fig.add_annotation(
            text = method_name,
            x = xx + 0.2,
            y = 1.0,
            xref = 'paper',
            yref = 'paper',
            xanchor = 'right',
            yanchor = 'bottom', 
            textangle = 20,
            font_size = 40,
            showarrow = False,
        )

    fig.add_annotation(
        text = 'Section',
        x = 0.5,
        y = -0.35,
        xref = 'paper',
        yref = 'paper',
        font_size = 40,
    )

    fig.update_layout(
        font_family="serif",
        font_size = 36,
        width = 1200,
        margin = dict(
            r = 150,
            b = 100,
            t = 115,
        ),
    )

    fig.write_image(f"{outfile_base_name}/{fn}.pdf")

if __name__ == "__main__":
    assert len(sys.argv) > 1, "Expected a CL arg giving the path to the results file to be visualized"

    outfile_base_name= sys.argv[1]
    fn = sys.argv[2]
    json_results_path = sys.argv[3]

    assert os.path.isfile(json_results_path), "File not found"



    with open(json_results_path) as f:
        visualize(json.load(f), outfile_base_name, fn)
