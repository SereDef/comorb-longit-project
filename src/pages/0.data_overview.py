import dash_bootstrap_components as dbc
from dash import dcc, html, register_page

from plotly.subplots import make_subplots
import plotly.graph_objects as go

import pandas as pd
import numpy as np

import definitions.layout_styles as styles
from definitions.variable_names import labels
from definitions.general_funcs import bold_it, badge_it, wrap_it

# Quick descriptives file
all_desc = pd.read_csv('./assets/descdf.csv', index_col=0)

# Only age variables
age_desc = all_desc[list(all_desc.columns[all_desc.columns.str.contains('age')])]


def plot_overview(width, exclude_cmr=('BMI', 'FMI', 'LMI', 'TFI', 'alcohol', 'canabis', 'smoking')):
    """Input: dimensions (optional).
       Creates an interactive scatterplot figure with time of measurement on the x-axis and all measures on the y-axis.
       Markers indicate the median age at measurement [+/- 95% CI in age range].
       These values + number of observations are also displayed when hovering over the points.
    """
    fig = make_subplots(2, 1, vertical_spacing=0.05, row_heights=[13 * 2, 25 * 2])

    colors = {'sDEP': '71, 107, 237', 'mDEP': '148, 103, 189', 'CMR': '214, 39, 40'}  # blue, purple and red
    legend = {'sDEP': 'Depression\nself-reports', 'mDEP': 'Depression\nmaternal reports',
              'CMR': 'Cardio-metabolic\nrisk markers'}

    def scat(vartype):
        # Get ages of measurement from age_desc
        df = age_desc[list(age_desc.columns[age_desc.columns.str.contains(vartype)])]

        if vartype == 'CMR':
            # Extract variable names from labels (only keeping those of interest)
            cmrkeys = [i for i in list(labels.keys())[14:] if i not in exclude_cmr]

            # Specify their label (including manually combining some measures into one row)
            cmrlab = {k: labels[k] for k in cmrkeys if k in labels}
            cmrlab['weight'] = 'Weight / BMI'
            cmrlab['total_fatmass'] = 'Total fat mass / FMI'
            cmrlab['total_leanmass'] = 'Total lean mass / LMI'

            names = list(cmrlab.values())
        else:
            names = list(labels.values())[0:13]  # depression items 1-13

        # Create grid of values and ranges
        x, y = [i.flatten() for i in np.meshgrid(df.loc['50%'], names)]
        x_min, _ = np.meshgrid(df.loc['5%'], names)
        x_max, _ = np.meshgrid(df.loc['95%'], names)

        # The CMR grid will need to contain some "NA" (and so do the counts)
        if vartype == 'CMR':
            ages = [i[-1][:-1] for i in df.columns.str.split('_')]  # get age without 'y'

            cmrmat = pd.DataFrame(columns=ages, index=cmrlab.keys())
            counts = cmrmat.copy()

            for age in ages:
                # what measures are available at that age
                t = all_desc[list(
                    all_desc.columns[all_desc.columns.str.contains(age) & ~(all_desc.columns.str.contains('DEP|age'))])]
                meas = [x.split('_' + age)[0] for x in t.columns if x.split('_' + age)[0] not in exclude_cmr]
                nobs = [t.loc['count', x] for x in t.columns if x.split('_' + age)[0] not in exclude_cmr]

                cmrmat.loc[meas, age] = [cmrlab[m] for m in meas]
                counts.loc[meas, age] = [int(c) for c in nobs]

            y = list(cmrmat.stack(future_stack=True))
            counts = list(counts.stack(future_stack=True))

        else:
            counts = []
            for i in [vartype[0] + s for s in list(labels.keys())[0:13]]:
                counts.extend([int(c) for c in
                               list(all_desc.loc['count', list(all_desc.columns[all_desc.columns.str.contains(i)])])])

        p = go.Scatter(x=x, y=y, mode='markers',
                       marker=dict(size=10, symbol='square', color=f'rgb({colors[vartype]})', opacity=.9),
                       error_x=dict(type='data', symmetric=False,
                                    array=x_max.flatten() - x,
                                    arrayminus=x - x_min.flatten(),
                                    color=f'rgba({colors[vartype]}, 0.3)',  # use this for opacity
                                    thickness=10, width=0),
                       name=legend[vartype],
                       text=[f'{round(i, 1)}-{round(j, 1)}' for i, j in zip(x_min.flatten(), x_max.flatten())],
                       customdata=[f'{i}' for i in counts],
                       hovertemplate="""<b>%{y}</b> <br> %{x:.1f} [%{text}] years <br> <i>N =</i> %{customdata} <br><extra></extra>""")
        return p

    fig.add_trace(scat('sDEP'), 1, 1)
    fig.add_trace(scat('mDEP'), 1, 1)
    fig.add_trace(scat('CMR'), 2, 1)

    cmrgroups = {'Anthropometry': [4, 'lightblue'],
                 'Fat distribution': [5, 'yellow'],
                 'Arteries': [4, 'orange'],
                 'Heart': [4, 'red'],
                 'Blood-based    <br>metabolic    <br>markers': [6, 'brown'],
                 'Inflammation': [2, 'blue']}

    cmrgrouplabels = []
    y0 = -0.5

    for g in cmrgroups:
        y1 = y0 + cmrgroups[g][0]
        cmrgrouplabels.append(dict(x0=-width/10000-.010, x1=-width/10000-.005, xref='paper',
                                   y0=y0, y1=y1 - 0.25, yref='y2',
                                   type='rect', fillcolor=cmrgroups[g][1], opacity=.3, line_width=0,
                                   label=dict(text=f'<b>{g}    </b>', textposition='top right', padding=2,
                                              font=styles.OVERVIEW_LABELS)
                                   ))
        y0 = y1 + 0.05

    fig.update_xaxes(range=[9, 26], mirror=True, ticks='outside', linecolor='black', gridcolor='lightgrey',
                     tickmode='linear', tick0=9, dtick=1)

    fig.update_yaxes(mirror=True, ticks='inside', linecolor='black')

    fig.update_layout(autosize=False, width=width, height=width*.55,
                      yaxis1=dict(range=(12.7, -0.7), tickfont=styles.OVERVIEW_LABELS),
                      yaxis2=dict(range=(24.7, -0.7), tickfont=styles.OVERVIEW_LABELS),
                      xaxis1=dict(tickfont=styles.OVERVIEW_LABELS),
                      xaxis2=dict(tickfont=styles.OVERVIEW_LABELS,
                                  title='Child age (years)', titlefont=styles.OVERVIEW_LABELS),
                      hoverlabel=dict(font_size=styles.OVERVIEW_LABELS['size']),
                      shapes=cmrgrouplabels, plot_bgcolor='white', margin=dict(l=20, r=20, t=20, b=20),
                      legend=dict(orientation='h', entrywidth=300, font=dict(size=16),
                                  yanchor='bottom', y=1.02, xanchor='left', x=0))

    return fig


register_page(__name__, path='/')

layout = dbc.Row([
    dbc.Col(width={'size': 10, 'offset': 1},
            children=[
                html.Div(['Overview of the data available from the', bold_it('ALSPAC'), 'cohort.', wrap_it(2),
                          'This includes 13 symptoms of depression from the Short Mood and Feelings Questionnaire (SMFQ), \
                          completed by mothers (', badge_it(' ', color='rgb(148, 103, 189)', pad='6px 7px'),
                          ') or participants themselves (', badge_it(' ', color='rgb(71, 107, 237)', pad='6px 7px'),
                          ') at several time points across 15 years.', wrap_it(), 'We also measured a total of 25 physical \
                          health markers (', badge_it(' ', color='rgb(214, 39, 40)', pad='6px 7px'), '), spanning from \
                          anthropometry to inflammation. Hover over graph below for more information about these data points.'],
                         style=styles.TEXT),
                html.Hr(),

                dcc.Graph(id='overview-fig',  # style={'width': '1000vw', 'height': '100vw'},
                          figure=plot_overview(width=styles.OVERVIEW_WIDTH))
                ])
])
