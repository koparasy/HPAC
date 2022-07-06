#!/g/g15/vanover1/GSA/.venv/bin/python

import sys
import os

import numpy as np

import plotly.graph_objects as go
import plotly.io as pio; pio.kaleido.scope.mathjax = None
from plotly.subplots import make_subplots
import plotly.express as px

def visualize(all_data, outfile_base_name, fn):

    technique_labels = {
        'adapt' : 'ADAPT',
        'Puppeteer' : 'Puppeteer',
    }
    benchmark_labels = {
        'DCT' : 'DCT',
        'Blackscholes' : 'Blackscholes',
        'HPCCG' : 'HPCCG',
    }
    technique_colors = {
        'adapt' : px.colors.qualitative.Plotly[2],
        'Puppeteer' : px.colors.qualitative.Plotly[1],
    }

    data_dict = {}
    row_i = 0
    while row_i < len(all_data):

        benchmark_name = all_data[row_i][0]
        framework_name = all_data[row_i][1]

        data_dict[benchmark_name] = {
            framework_name : {
                'Peak Memory' : all_data[row_i][2],
                'Execution Time' : all_data[row_i][3],
            }
        }

        row_j = row_i + 1
        while row_j < len(all_data) and all_data[row_j][0] == benchmark_name:

            framework_name = all_data[row_j][1]

            data_dict[benchmark_name][framework_name] = {
                        'Peak Memory' : all_data[row_j][2],
                        'Execution Time' : all_data[row_j][3],
            }

            row_j += 1

        row_i = row_j

    aggregated_data = {}
    for benchmark_name, data in data_dict.items():

        original = data['original']

        temp = {}
        for framework, ddata in data.items():
            if framework == 'original':
                continue

            temp[framework] = {
                'Peak Memory' : ddata['Peak Memory']/original['Peak Memory'],
                'Execution Time' : ddata['Execution Time']/original['Execution Time'],
            }

        aggregated_data[benchmark_name] = temp

    fig = make_subplots(
        rows = 1,
        cols = 2,
        shared_yaxes = True,
        subplot_titles = ['Execution Time', 'Peak Memory Usage'],
        horizontal_spacing = 0.05,
    )

    for technique_name in list(aggregated_data.values())[0].keys():
        benchmark_names = [benchmark_name for benchmark_name in aggregated_data.keys()]
        fig.add_trace(
            go.Bar(
                y = [benchmark_labels[x] for x in benchmark_names],
                x = np.log([y['Execution Time'] for y in [x[technique_name] for x in [aggregated_data[benchmark_name] for benchmark_name in benchmark_names]]]),
                orientation = 'h',
                name = technique_labels[technique_name],
                marker_color = technique_colors[technique_name],
            ),
            row = 1,
            col = 1,
        )

    for technique_name in list(aggregated_data.values())[0].keys():
        benchmark_names = [benchmark_name for benchmark_name in aggregated_data.keys()]
        fig.add_trace(
            go.Bar(
                y = [benchmark_labels[x] for x in benchmark_names],
                x = np.log([y['Peak Memory'] for y in [x[technique_name] for x in [aggregated_data[benchmark_name] for benchmark_name in benchmark_names]]]),
                orientation = 'h',
                name = technique_labels[technique_name],
                showlegend = False,
                marker_color = technique_colors[technique_name],
            ),
            row = 1,
            col = 2,
        )

    # for subplot titles
    fig.update_annotations(font_size=40)

    fig.update_yaxes(
        gridcolor = "#A9A9A9",
        showline = True,
        linewidth = 1,
        linecolor = "#2a3f5f",
    )
    fig.update_xaxes(
        gridcolor = "#A9A9A9",
        showline = True,
        linewidth = 1,
        linecolor = "#2a3f5f",
    )

    fig.update_layout(
        font_family = "serif",
        font_size = 30,
        plot_bgcolor = "#fff",
        height = 450,
        width = 1000,
        legend = dict(
            x = 0.25,
            y = 0.93,
            bordercolor = "#2a3f5f",
            borderwidth = 1,
        ),
        xaxis = dict(
            tickmode = 'array',
            tickvals = np.log([1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64, 100.0]),
            ticktext = [' ', '2x', ' ', '8x', ' ', '32x', ' ', '100x'],
            tickangle = 45,
        ),
        xaxis2 = dict(
            tickmode = 'array',
            tickvals = np.log([1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0]),
            ticktext = [' ', '2x', ' ', ' ', '16x', ' ', '64x', ' ', '256x'],
            tickangle = 45,
        )
    )

    fig.write_image(f"{outfile_base_name}/{fn}.pdf")


if __name__ == "__main__":
    assert len(sys.argv) > 3, "Expected CL args giving the path to the results file to be visualized"
    assert(os.path.isfile(sys.argv[3]))

    outfile_base_name = sys.argv[1]
    fn = sys.argv[2]
    input_file_path = sys.argv[3]

    with open(input_file_path) as f:
        visualize(np.genfromtxt(input_file_path, delimiter = ":", skip_header = 1, dtype = None, encoding = None), outfile_base_name, fn)
