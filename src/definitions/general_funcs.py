from dash import html
import dash_bootstrap_components as dbc

import plotly.graph_objects as go

import json
import pandas as pd
import re
import textwrap

import definitions.layout_styles as styles
from definitions.variable_names import labels


# -- Frontend --------------------


def bold_it(text):
    return html.Span(f' {text} ', style={'font-weight': 'bold'})


def underline_it(text):
    return html.Span(f'{text}', style={'text-decoration': 'underline'})


def badge_it(text, color, fs=styles.TEXT['font-size'], pad='4px 5px'):
    return dbc.Badge(text, color=color, style={'padding': pad, 'font-size': fs})


def wrap_it(n_spaces=1):
    return html.Span([html.Br()] * n_spaces)


# -- Backend --------------------

def get_label(var):
    # Remove year from var name if present
    vs = var.split('_')
    if bool(re.match('[0-9]', vs[-1])):
        var_name = '_'.join(vs[:-1])
    else:
        var_name = '_'.join(vs)
    # Remove self or maternal report from depression item
    if 'DEP' in var_name:
        var_name = var_name[1:]
    # get label
    lab = labels[var_name]
    return lab


# Opening Descriptives JSON file
with open('./assets/descrip.json') as jf:
    desc = json.load(jf)


def desc_plot(var):
    # Read in descriptives from file
    d = desc[var]
    summ, count, dens = (pd.DataFrame.from_dict(x) for x in d)

    # Initialize figure
    fig = go.Figure()

    if bool(re.match(r'.DEP[0-9]+|alcohol|canabis|smoking', var)):

        likert = ['Not true', 'Sometimes', 'True']
        colors = ['lightblue', 'salmon', 'crimson']

        x_lab = '<br>'.join(textwrap.wrap(get_label(var), width=20))  # Wrap label text

        # Stacked histogram
        for i, c in enumerate(count['count'].iloc[:-1]):
            fig.add_trace(go.Bar(x=[0], y=[c], name=likert[i], marker_color=colors[i], opacity=0.6,
                                 marker_line_width=1.5, marker_line_color='white',
                                 customdata=[count['prop_noNA'].iloc[i] * 100],
                                 hovertemplate="""<br> n = <b>%{y}</b> (%{customdata:.1f}%) <br><extra></extra>"""))

            fig.update_layout(barmode='stack', autosize=False, width=400, height=500,
                              margin=dict(l=20, r=20, t=20, b=20),
                              yaxis=dict(title_text='Count'),
                              xaxis=dict(ticktext=['', x_lab, ''], tickvals=[-0.5, 0, 0.5], tickfont=dict(size=15)))
        return fig

    else:
        x_lab = get_label(var)

        # Distribution plot
        fig.add_trace(go.Scatter(x=dens['x'], y=dens['dens'], fill='tozeroy', mode='lines', line_color='salmon'))
        # add mean and median
        fig.add_vline(x=summ['mean'][0], line_width=1.5, line_color='crimson', name='meanvalue')
        fig.add_vline(x=summ['50%'][0], line_width=1, line_dash='dash', line_color='crimson', name='medianvalue')
        # Adjust x-axis
        fig.update_xaxes(title_text=x_lab, range=[round(summ['min']), round(summ['max'])],
                         ticks='outside', showline=True, linecolor='black', gridcolor='lightgrey')
        fig.update_layout(autosize=False, width=600, height=300, margin=dict(l=20, r=20, t=20, b=20),
                          yaxis=dict(showticklabels=False, showline=True, linecolor='black'),
                          )
        return fig
