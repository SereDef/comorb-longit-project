from dash import html
import dash_bootstrap_components as dbc
from dash_iconify import DashIconify

from plotly.subplots import make_subplots
import plotly.graph_objects as go

import json
import pyreadr
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


def info_icon(c=styles.NAVBAR_COLOR, w=15):
    return DashIconify(icon='entypo:info-with-circle', color=c, width=w)


def dep_var_checklist(disable_maternal=False):
    """Defines the radio buttons used to select self vs. maternal reports of depression."""
    opts = [{'label': 'Self-reported', 'value': 'sDEP', 'disabled': False},
            {'label': 'Maternal report', 'value': 'mDEP', 'disabled': disable_maternal}]
    return opts


def cmr_var_checklist(reporter='s'):
    """Defines the dropdown menu used to select the CMR marker to model."""

    # Manually define the order of the markers to display
    sortd = ['FMI', 'LMI', 'BMI', 'waist_circ', 'android_fatmass', 'total_fatmass', 'total_leanmass',
             'tot_chol', 'HDL_chol', 'LDL_chol', 'insulin', 'triglyc', 'CRP']

    opts = [{'value': v, 'label': get_label(v) + f' ({v})' if len(v) < 4 else get_label(v), 'disabled': False} for v in
            sortd]

    # Disable options for which no maternal depression model was estimated
    if reporter == 'm':
        sortd_m = ['FMI', 'LMI', 'BMI', 'waist_circ', 'total_fatmass', 'total_leanmass']
        for i in opts:
            if i['value'] not in sortd_m:
                i['disabled'] = True

    return opts

def stat_checklist(disable_etas=False):
    """Defines the checklist used to set stationarity assumptions."""
    opts = [{'label': 'Between-person (latent) effect','value': 'l', 'disabled': disable_etas},
            {'label': 'Within-person parameters','value': 'p'}]
    return opts

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


def temp_plot(depname, cmrname):
    """Input: names of the depression report (sDEP = self or mDEP = parental reports) and cardio-metabolic risk (CMR)
       marker.
       Creates an interactive scatterplot figure with time of measurement on the x-axis and values of
       dep and cmr markers on the (double) y-axis. Markers indicate the median [IQR] of each measure. Median values and
       time points are also displayed when hovering over the points.
    """
    # load summary dataframe
    summ = pyreadr.read_r(f'./assets/results/mod1/g_lstat_pstat/{depname}_{cmrname}.RData')['dat_summ']

    # extract timepoints
    t_dep = [float(x.split('_')[-1][:-1]) for x in summ.columns[:summ.shape[1] // 2]]
    t_cmr = [float(x.split('_')[-1][:-1]) for x in summ.columns[summ.shape[1] // 2:]]

    # scatterplot function
    def scat(t, name, fullname, shortname):
        means = summ.loc['Median', summ.columns.str.contains(name)]
        p = go.Scatter(x=t, y=means,
                       error_y=dict(type='data', symmetric=False,  # visible=True,
                                    array=summ.loc['3rd Qu.', summ.columns.str.contains(name)] - means,
                                    arrayminus=means - summ.loc['1st Qu.', summ.columns.str.contains(name)]),
                       name=fullname, text=[f'{shortname} {n}' for n in range(1, len(t) + 1)],
                       marker=dict(size=10, symbol='square', opacity=.8), opacity=.7,
                       hovertemplate="""<b>%{text}</b> <br> Median: %{y:.2f} <br> Timepoint: %{x} years \
                       <br><extra></extra>""")
        return p

    fig = make_subplots(
        specs=[[{'secondary_y': True}]])  # Specify a double y axis to allow for different scale of depression and CMR

    fig.add_trace(scat(t_dep, depname, 'Depression score', 'DEP'), secondary_y=False)
    fig.add_trace(scat(t_cmr, cmrname, cmrname, cmrname), secondary_y=True)

    # Set y-axes
    def yrange(name):
        sub = summ.filter(like=name)
        ymin = sub.min(axis=1)['1st Qu.']
        ymax = sub.max(axis=1)['3rd Qu.']
        rang = ymax - ymin
        y_max_lower = ymax
        ymin = ymin - (rang / 10)
        ymax = ymax + (rang / 10)
        return [ymin, ymax, y_max_lower]

    fig.update_yaxes(title_text='<b>Depression score</b>', secondary_y=False, range=yrange(depname)[:2],
                     mirror=True, ticks='outside', showline=True, linecolor='black', gridcolor='lightgrey')
    fig.update_yaxes(title_text=f'<b>{cmrname}</b>', secondary_y=True, range=yrange(cmrname)[:2],
                     mirror=True, ticks='outside', showline=True, linecolor='black', gridcolor='lightgrey')
    # Set x-axis
    fig.update_xaxes(title_text='Years', mirror=True, ticks='outside', showline=True, linecolor='black',
                     gridcolor='lightgrey')

    # Group "cross-sectional" time points using grey background
    crosspoints = []
    for i in range(len(t_dep)):
        xmin = min(t_dep[i], t_cmr[i]) - .2
        xmax = max(t_dep[i], t_cmr[i]) + .2
        # rectangles
        crosspoints.append(dict(x0=str(xmin), x1=str(xmax), y0=0, y1=1, xref='x', yref='paper',
                                type='rect', fillcolor='lightgray', opacity=.3, line_width=0, layer='below'))
        # text
        fig.add_trace(go.Scatter(x=[xmin + (xmax - xmin) / 2], y=[yrange(depname)[2]], mode='text', text=i + 1,
                                 textposition='top center',
                                 textfont_size=13, textfont_color='dimgray', showlegend=False))
    # Background
    fig.update_layout(  # title = dict(text='Included measures\n', font=dict(size=15), automargin=True, yref='paper'),
        plot_bgcolor='white', shapes=crosspoints, margin=dict(l=20, r=20, t=20, b=20))

    return fig


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
