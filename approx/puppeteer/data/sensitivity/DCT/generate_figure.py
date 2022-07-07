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
import plotly.express as px

def visualize_json(all_data, outfile_base_name, fn):

    ordered_kernels = np.argsort([int(x.split("_")[-1]) for x in all_data["Puppet-With-Quant"]["names"]])

    #method_names = ( "SALib-Morris", "SALib-Sobol","IVARS")
    method_names = ("Puppet-Without-Quant","Puppet-With-Quant","ADAPT-With-Quant")


    names = ["Puppeteer <br> Without-Quantization","Puppeteer <br> With-Quantization", "ADAPT  <br> With-Quantization"]
    fig = make_subplots(
        rows = 1,
        cols = 3,
        subplot_titles = names,
        shared_yaxes = True,
        horizontal_spacing = 0.01,
    )

    min_values = 0.0
    for trace_i, method_name in enumerate((method_names)):
        data_dict = all_data[method_name]

        # gather all sensitivities and store them in an 8x8 grid
        sensitivities = []
        for i in ordered_kernels:
            sensitivities.append(data_dict["Sensitivity"][i])

        sensitivities = np.reshape(sensitivities, (8,8))

        # give sensitivities a logarithmic scale
        warnings.filterwarnings("ignore", "divide by zero encountered in log")
        sensitivities_log = np.log(sensitivities)
        warnings.resetwarnings()

        # if any of the sensitivities were 0, impute the -inf values
        # given by np.log(0) with twice the magnitude of the log of the
        # smallest nonzero sensitivity
        if np.all(sensitivities == 0):
#          sensitivities_log.zeros(sensitivitie.shape)
          sensitivities_log_norm = np.zeros(sensitivities.shape) #sensitivities_log
        elif (not sensitivities.all()):
            sorted_flattened_sensitivities = sorted(sensitivities.flatten(), reverse=True)
            second_smallest = sorted_flattened_sensitivities.pop()
            while second_smallest == 0:
                second_smallest = sorted_flattened_sensitivities.pop()
            sensitivities_log[np.where(sensitivities == 0)] = np.log(second_smallest) * 2
            zero_flag = True
            sensitivities_log_norm = sensitivities_log + np.abs(np.min(sensitivities_log))
            sensitivities_log_norm = sensitivities_log_norm / np.max(sensitivities_log_norm)
            if min_values > np.min(sensitivities_log_norm):
              min_values = np.min(sensitivities_log_norm)
        else:
            zero_flag = False
            sensitivities_log_norm = sensitivities_log + np.abs(np.min(sensitivities_log))
            sensitivities_log_norm = sensitivities_log_norm / np.max(sensitivities_log_norm)
            if min_values > np.min(sensitivities_log_norm):
              min_values = np.min(sensitivities_log_norm)

        # normalize sensitivities by scaling them so the min is 0 and
        # then dividing by the maximum value

        # plot each heatmap
        if 'ADAPT' in method_name:
          fig.add_trace(
              go.Heatmap(
                  x = [str(x) for x in range(8)],
                  y = [str(x) for x in range(8)],
                  z = sensitivities_log_norm,
                  zmin=0,
                  zmax=1,
#                  name =  names[trace_i],
                  colorbar = dict(
                      tick0 = 0,
                      dtick = 1.0,
                  ),

              ),
              row = 1,
              col = trace_i + 1
          )
        else:
          fig.add_trace(
              go.Heatmap(
                  x = [str(x) for x in range(8)],
                  y = [str(x) for x in range(8)],
                  z = sensitivities_log_norm,
#                  name =  names[trace_i],
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
    fig.update_annotations(font_size=34)

    fig.add_annotation(
        text = 'Sensitivity',
        x = 1.1,
        y = 0.4,
        xref = 'paper',
        yref = 'paper',
        textangle = 90,
        font_size = 40,
    )

    tick_descr = dict(
            ticks = "outside",
            tickwidth = 3,
            tickcolor = "#2a3f5f",
            ticklen = 10,
            tickvals = list(range(8)),
            ticktext = ["0", " ", "2", " ", "4", " ", "6", " ", "7"],
        )

    fig.update_layout(
        font_family="serif",
        font_size = 40,
        width = 1250,
        height = 500,
        margin = dict(
            r = 150,
        ),
        xaxis = tick_descr,
        xaxis2 = tick_descr,
        xaxis3 = tick_descr,
        yaxis = tick_descr,
    )

    fig.write_image(f"{outfile_base_name}/{fn}.pdf")


def visualize_dat(all_data, outfile_base_name):

    NUM_RUNS_PER_FRONTIER = 5

    aggregated_data = []
    for i in range(0, all_data.shape[0], NUM_RUNS_PER_FRONTIER):
        aggregated_data.append([all_data[i:i + NUM_RUNS_PER_FRONTIER,2].mean(), all_data[i,3]])
    aggregated_data = np.array(aggregated_data)

    # normalize
    baseline_runtime = aggregated_data[-1][0]
    runtime_to_speedup_logspace = lambda x : np.array([np.log(baseline_runtime / x[0]), x[1]])
    aggregated_data = np.apply_along_axis(runtime_to_speedup_logspace, axis=1, arr=aggregated_data)

    fig = make_subplots(specs=[[{"secondary_y" : True}]])
    fig.add_trace(
        go.Scatter(
            x = list(range(2,len(aggregated_data) + 1)),
            y = aggregated_data[1:,0],
            name = "Speedup",
        ),
        secondary_y = True,
    )

    fig.add_trace(
        go.Scatter(
            x = list(range(2,len(aggregated_data) + 1)),
            y = aggregated_data[1:,1],
            name = "Quality (PSNR)<br><sup>(Larger is better)</sup>",
        ),
        secondary_y = False,
    )

    # Set y-axes titles
    fig.update_yaxes(
        showline = True,
        linewidth = 1,
        linecolor = "#2a3f5f",
        color = px.colors.qualitative.Plotly[0],
        secondary_y = True,
        ticks = "outside",
        tickwidth = 1,
        tickcolor = "#2a3f5f",
        ticklen = 5,
        tickmode = 'array',
        tickvals = np.log([1.0, 2.0, 4.0, 8.0, 20.0]),
        ticktext = ['1x', '2x', '4x', '8x', '20x'],
    )
    fig.update_yaxes(
        showline = True,
        linewidth = 1,
        linecolor = "#2a3f5f",
        color = px.colors.qualitative.Plotly[1],
        secondary_y = False,
        ticks = "outside",
        tickwidth = 1,
        tickcolor = "#2a3f5f",
        ticklen = 5,
    )
    fig.update_xaxes(
        ticks = "outside",
        tickwidth = 1,
        tickcolor = "#2a3f5f",
        ticklen = 5,
        showline = True,
        linewidth = 1,
        linecolor = "#2a3f5f",
        title_text = "Diagonal",
        tick0 = 2,
        dtick = 2,
    )
    fig.update_layout(
        font_family="serif",
        font_size = 24,
        plot_bgcolor = "#fff",
        height = 420,
        legend = dict(
            x = 0.5,
            y = 0.8,
            bordercolor = "#2a3f5f",
            borderwidth = 1,
        )
    )

    fig.write_image(f"{outfile_base_name}_frontier_figure{fig_no}.pdf")

if __name__ == "__main__":
    assert len(sys.argv) > 1, "Expected a CL arg giving the path to the results file to be visualized"

    outfile_base_name= sys.argv[1]
    fn = sys.argv[2]
    json_results_path = sys.argv[3]
    assert os.path.isfile(json_results_path), "File not found"

    with open(json_results_path) as f:
      visualize_json(json.load(f), outfile_base_name, fn)

