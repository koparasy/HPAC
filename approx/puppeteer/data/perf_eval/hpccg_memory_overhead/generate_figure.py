#!/g/g15/vanover1/GSA/.venv/bin/python

import sys
import os

import numpy as np

import plotly.graph_objects as go
import plotly.io as pio; pio.kaleido.scope.mathjax = None
import plotly.express as px

def eprint(*args, **kwargs):
      print(*args, file=sys.stderr, **kwargs)


def visualize(all_data, outfile_base_name, fn):

    # first three entries in each row give dimensions of the problem;
    # take their product for the total size
    # last entry in each row gives the peak memory consumption, the
    # only quantity we are concerned with
    aggregated_data = {}
    for input_file_name, data in all_data.items():

        # store aggregated data, associating it with the name that
        # will be used for the legend
        if "accurate" in input_file_name:
            descriptor = "Original"
        else:
            descriptor = "Puppeteer"

        aggregated_data[descriptor] = np.array([[row[:3].prod(), row[-1]] for row in data])
        
    # add ADAPT
    aggregated_data['ADAPT'] = np.array([[64**3, 28999514479], [128**3, 0], [192**3, 0], [256**3, 0], [320**3, 0]])
    
    fig = go.Figure()
    for descriptor, data in aggregated_data.items():
        fig.add_trace(
            go.Bar(
                x = ["64<sup>3</sup>", "128<sup>3</sup>", "192<sup>3</sup>", "256<sup>3</sup>", "320<sup>3</sup>"],
                y = data[:,-1]/1073741824, # divide by GB
                name = descriptor,
            )
        )

    x_start = 0.128
    for i in range(4):
        x_start += 0.2
        fig.add_shape(
            type = "rect",
            xref = "paper",
            yref = "paper",
            x0 = x_start,
            y0 = 0,
            x1 = x_start + 0.05,
            y1 = 0.95,
            line_dash = "dot",
            line_color = "black",
            fillcolor = px.colors.qualitative.Plotly[2],
            opacity = 0.3,
        )
        if i == 2:
            text = 'Out of Memory'
            y = 0.175
        else:
            text = 'Out of Memory (>128GB)'
            y = 0.365
        fig.add_annotation(
            text = text,
            x = x_start + 0.04,
            y = y,
            xref = 'paper',
            yref = 'paper',
            textangle = 270,
            font_size = 24,
        )

    fig.update_yaxes(
        title_text = "Peak Memory Usage (GB)",
        gridcolor = "#A9A9A9",
        showline = True,
        linewidth = 1,
        linecolor = "#2a3f5f",
    )
    fig.update_xaxes(
        title_text = "Problem Dimension",
        gridcolor = "#A9A9A9",
        showline = True,
        linewidth = 1,
        linecolor = "#2a3f5f",
    )    
    fig.update_layout(
        font_family = "serif", 
        font_size = 20,
        height = 450,
        plot_bgcolor = "#fff",
        legend = dict(
            x = 0.65,
            y = 0.94,
            bordercolor = "#2a3f5f",
            borderwidth = 1,
        )
    )

    fig_no = 0
#    while os.path.isfile(f"{outfile_base_name}_peak_mem_figure{fig_no}.png"):
#        fig_no += 1
    fig.write_image(f"{outfile_base_name}/{fn}.pdf")


def visualize_scaling(all_data, outfile_base_name, fn):

    aggregated_data = {}
    for input_file_name, data in all_data.items():
        print(input_file_name)
        # store aggregated data, associating it with the name that
        # will be used for the legend
        if "accurate" in input_file_name:
            descriptor = "Original"
        else: 
            descriptor = "Puppeteer"
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
        aggregated_data[descriptor] = np.array(temp)
        print(descriptor)
        print(aggregated_data)
    print(aggregated_data)
    fig = go.Figure()
    for descriptor, data in aggregated_data.items():
        fig.add_trace(
            go.Scatter(
                x = data[:,0],
                y = data[:,1],
                name = descriptor,
            )
        )

    fig.update_yaxes(
        title_text = "Execution Time (Seconds)",
        gridcolor = "#A9A9A9",
        showline = True,
        linewidth = 1,
        linecolor = "#2a3f5f",
    )
    fig.update_xaxes(
        title_text = "Number of Threads",
        gridcolor = "#A9A9A9",
        showline = True,
        linewidth = 1,
        linecolor = "#2a3f5f",
        tickmode = "array",
        tickvals = [1, 2, 4, 8, 16, 32]
    )    
    fig.update_layout(
        font_family = "serif", 
        font_size = 18,
        plot_bgcolor = "#fff",
        legend = dict(
            x = 0.60,
            y = 0.75,
            bordercolor = "#2a3f5f",
            borderwidth = 1,
        )
    )

    fig_no = 0
#    while os.path.isfile(f"{outfile_base_name}_figure{fig_no}.png"):
#        fig_no += 1
    # fig.write_image(f"{outfile_base_name}_thread_scaling_figure{fig_no}.png")
    fig.write_image(f"{outfile_base_name}/{fn}.pdf")


if __name__ == "__main__":
    assert len(sys.argv) > 1, "Expected CL args giving the paths to the results files to be visualized"

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
    elif all(["memory" in x for x in all_data.keys()]):
        visualize(all_data, sys.argv[1], sys.argv[2])
    else:
        assert False, "Mix of scaling and non-scaling log files not allowed"
