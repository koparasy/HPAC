#!/g/g15/vanover1/GSA/.venv/bin/python

import sys
import os
import json

import numpy as np

import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio; pio.kaleido.scope.mathjax = None

def getTechniqueOrder():
  return ['SALib-Morris', 'SALib-Sobol','IVARS', 'VARS-Morris-ABE', 'VARS-Sobol-TO']


def visualize(all_data, outfile_base_name, normalize, fn=''):

    var_names = ["CNDF1", "CNDF2", "EXP", "LOG", "SQRT"]
    techniques = list(all_data.keys())

    if normalize:
        fn = fn + "_normalized"
        for technique in techniques:
            unnormalized_sens = all_data[technique]['Sensitivity']
            normalized_sens = unnormalized_sens / np.max(unnormalized_sens)
            t = technique.replace("-TO", "").replace("-ABE", "")
            all_data[t] = all_data.pop(technique)
            all_data[t]["Sensitivity"] = normalized_sens
    else:
        fn = fn + "_unnormalized"
        for technique in getTechniqueOrder():
            t = technique.replace("-TO", "").replace("-ABE", "")
            all_data[t] = all_data.pop(technique)

    sensitivities = np.array([np.array([x['Sensitivity'] for x in all_data.values()])[:,i] for i in range(len(var_names))])
    techniques = list(all_data.keys())

    np.random.seed(2)
    point_offsets = np.linspace(-0.2,0.2, num=len(var_names))
    np.random.shuffle(point_offsets)

    fig = go.Figure()

    for i in range(len(var_names)):
        fig.add_trace(
            go.Box(
                name = var_names[i],
                x = techniques,
                y = sensitivities[i],
                boxpoints = "all",
                pointpos = point_offsets[i],
                marker = dict(
                    color = px.colors.qualitative.G10[i],
                    size = 13),
                line = dict(color = 'rgba(0,0,0,0)'),
                fillcolor = 'rgba(0,0,0,0)',
            )
        )

    fig.update_yaxes(
        type="log",
        showline = True,
        linewidth = 1,
        linecolor = "#2a3f5f",
        gridcolor = "#A9A9A9",
    )
    fig.update_xaxes(
        showline = True,
        linewidth = 1,
        linecolor = "#2a3f5f",
    )
    fig.update_layout(
        yaxis_title="Sensitivity",
        font_family="serif",
        font_size = 24,
        plot_bgcolor = "#fff",
        legend = dict(
            orientation="h",
            x = -0.01,
            y = 1.20,
            bordercolor = "#2a3f5f",
            borderwidth = 1,
        ),
        yaxis_tickformat = "~e",
    )

    fig_no = 0
    fig.write_image(f"{outfile_base_name}/{fn}.pdf")


if __name__ == "__main__":
    assert len(sys.argv) > 3, "Expected a CL arg giving the path to the JSON results file to be visualized"

    outfile_base_name= sys.argv[1]
    fn = sys.argv[2]
    json_results_path = sys.argv[3]
    assert os.path.isfile(json_results_path), "File not found"

    with open(json_results_path) as f:
        all_data = json.load(f)
        visualize(all_data, outfile_base_name, False, fn)
#        visualize(all_data, outfile_base_name, normalize=True)
