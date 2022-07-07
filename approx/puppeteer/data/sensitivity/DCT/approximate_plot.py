import pandas as pd
import numpy as np
import seaborn as sns
import sys
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio; pio.kaleido.scope.mathjax = None
import plotly.express as px

outfile_base_name= sys.argv[1]
fn = sys.argv[2]
data = sys.argv[3]

df = pd.read_csv(data)
xAxisName = 'Number of perforated Kernels'
df[xAxisName] = df[df.columns[1]]
df=df.drop(df.columns[1], axis=1)
df = df.groupby(['Analysis',xAxisName]).mean().reset_index()
accurate=df[df['Analysis'] == 'Accurate']['Execution Time'].values[0]
df['Speedup'] = accurate/df['Execution Time']
yAxisName = 'PSNR'
df['PSNR'] = df['PSNR'].astype(float)
df[xAxisName] = df[xAxisName].astype(int)
df =df[df['Analysis'] != 'Accurate']

fig = make_subplots(
    rows = 2,
    cols = 1,
    vertical_spacing = 0.05,
    shared_xaxes = True,
)

markers=[123, 138]
ncolumns = {'WQ':'Puppeteer With-Quantization', 'WOQ':'Puppeteer Without-Quantization'}
df = df.replace({'Analysis':ncolumns})
for i, metric in enumerate(['Speedup', 'PSNR']):
  for j, analysis in enumerate(df['Analysis'].unique()):
    showlegend = True
    if i != 0:
      showlegend = False

    fig.add_trace(
        go.Scatter(
            x = df[df['Analysis'] == analysis][xAxisName],
            y = df[df['Analysis'] == analysis][metric],
            marker_color = px.colors.qualitative.Plotly[2+j],
            marker = dict(
              symbol=markers[j],
              size = 9,
              line=dict(width=2)
            ),
            mode='markers',
            name = analysis,
            showlegend = showlegend,
        ),
        row = 3 - 2**i,
        col = 1,
    )

fig.update_yaxes(
    title_standoff = 0,
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
    zeroline = True,
    zerolinecolor = "#A9A9A9",
    zerolinewidth = 0.8,
)

fig.update_layout(
    font_family = "serif",
    font_size = 18,
    plot_bgcolor = "#fff",
    legend = dict(
        title='Analysis Method',
        x = 0.1,
        y = 0.34,
        bordercolor = "#2a3f5f",
        borderwidth = 1,
        font_size = 25,
    ),
    xaxis2 = dict(
        title = dict(
            text = 'Number of Perforated Kernels',
            font_size = 25
        )
    ),
    yaxis2 = dict(
        title = dict(
            text = 'Speedup',
            font_size = 25
        ),
        tickmode = 'array',
        zeroline = True,
        zerolinecolor = "#A9A9A9",
        zerolinewidth = 0.8,
    ),
    yaxis1 = dict(
        title = dict(
            text = 'PSNR (db)',
            font_size = 25
        ),
        tickmode = 'array',
    )
)

fig.write_image(f"{outfile_base_name}/{fn}.pdf")

