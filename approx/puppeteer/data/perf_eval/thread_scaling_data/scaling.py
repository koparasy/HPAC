#!/g/g15/vanover1/GSA/.venv/bin/python

import sys
import os

import numpy as np

import plotly.graph_objects as go
import plotly.io as pio; pio.kaleido.scope.mathjax = None
import plotly.express as px
from plotly.subplots import make_subplots

def visualize_scaling(all_data, outfile_base_name, fn):

    aggregated_data = {}
    for input_file_name, data in all_data.items():
        # store aggregated data, associating it with the name that
        # will be used for the legend
        if "accurate" in input_file_name:
            technique = "Original"
        else:
            technique = "Puppeteer"

        if "hpccg" in input_file_name:
            benchmark = "HPCCG"
        else:
            benchmark = "Blackscholes"

        if benchmark not in aggregated_data.keys():
            aggregated_data[benchmark] = {}

        temp = []
        i = 0
        while i < data.shape[0]:
            execution_times = [data[i][2]]
            j = i + 1
            while j < data.shape[0] and data[i][0] == data[j][0]:
                execution_times.append(data[j][2])
                j += 1
            temp.append([data[i][0], np.mean(execution_times)])
            i = j
        aggregated_data[benchmark][technique] = np.array(temp)

    # # normalize wrt single thread
    # for benchmark in aggregated_data.keys():

    #     normalizer = aggregated_data[benchmark]['Original'][0][1]

    #     for technique in aggregated_data[benchmark].keys():

    #         normalized_data = []

    #         for data in aggregated_data[benchmark][technique]:
    #             normalized_data.append([data[0], data[1]/normalizer])
            
    #         aggregated_data[benchmark][technique] = np.array(normalized_data)

    fig = make_subplots(
        rows = 2,
        cols = 1,
        vertical_spacing = 0.05,
        shared_xaxes = True,
    )

    for benchmark_i, benchmark_name in enumerate(aggregated_data.keys()):
        for technique_i, technique_name in enumerate(aggregated_data[benchmark_name].keys()):
            
            if benchmark_i != 0:
                showlegend = False
            else:
                showlegend = True

            fig.add_trace(
                go.Scatter(
                    x = aggregated_data[benchmark_name][technique_name][:,0],
                    y = np.log(aggregated_data[benchmark_name][technique_name][:,1]),
                    marker_color = px.colors.qualitative.Plotly[technique_i],
                    marker_size = 12,
                    name = technique_name,
                    showlegend = showlegend,
                ),
                row = 3 - 2**benchmark_i,
                col = 1,
            )

    fig.update_yaxes(
        gridcolor = "#A9A9A9",
        showline = True,
        linewidth = 1,
        linecolor = "#2a3f5f",
        zeroline = True,
    )
    fig.update_xaxes(
        gridcolor = "#A9A9A9",
        showline = True,
        linewidth = 1,
        linecolor = "#2a3f5f",
        tickmode = "array",
        tickvals = [1, 2, 4, 8, 16, 32],
        zeroline = True,
        zerolinecolor = "#A9A9A9",
        zerolinewidth = 0.8,
    )
    fig.update_layout(
        font_family = "serif",
        font_size = 18,
        plot_bgcolor = "#fff",
        legend = dict(
            x = 0.65,
            y = 0.98,
            bordercolor = "#2a3f5f",
            borderwidth = 1,
            font_size = 25,
        ),
        xaxis2 = dict(
            title = dict(
                text = 'Number of Threads',
                font_size = 25
            )
        ),
        yaxis2 = dict(
            tickmode = 'array',
            tickvals = np.log([1/16, 1/8, 1/4, 1/2, 1, 2, 4]),
            ticktext = ['0.0625', '0.125', '0.25', '0.5', '1.0', '2.0', '4.0'],
            zeroline = True,
            zerolinecolor = "#A9A9A9",
            zerolinewidth = 0.8,
        ),
        yaxis1 = dict(
            tickmode = 'array',
            tickvals = np.log([4, 8, 16, 32, 64]),
            ticktext = ['4', '8', '16', '32', '64'],
        )
    )

    fig.add_annotation(
        text = 'Blackscholes',
        x = 0.87,
        y = 0.28,
        xref = 'paper',
        yref = 'paper',
        font_size = 32,
        bgcolor = "#fff",
        showarrow = False,
    )

    fig.add_annotation(
        text = 'HPCCG',
        x = 0.376,
        y = 0.93,
        bgcolor = "#fff",
        xref = 'paper',
        yref = 'paper',
        font_size = 32,
        showarrow = False,
    )

    fig.add_annotation(
        text = 'Execution Time (seconds)',
        x = -0.1,
        y = 0.5,
        xref = 'paper',
        yref = 'paper',
        textangle = 270,
        font_size = 25,
    )

    fig.write_image(f"{outfile_base_name}/{fn}.pdf")

if __name__ == "__main__":
    assert len(sys.argv) > 2, "Expected CL args giving the paths to the results files to be visualized"

    input_file_paths = []
    for arg_i in range(3, len(sys.argv)):
        assert(os.path.isfile(sys.argv[arg_i]))
        input_file_paths.append(sys.argv[arg_i])


    all_data = {}
    for input_file_path in input_file_paths:
        input_file_name = os.path.basename(input_file_path)
        with open(input_file_path) as f:
            all_data[input_file_name] = np.genfromtxt(input_file_path, delimiter = ":", skip_header = 1)

    if all(["scaling" in x for x in all_data.keys()]):
        visualize_scaling(all_data, sys.argv[1], sys.argv[2])
    else:
        assert False, "Only scaling log files allowed"
