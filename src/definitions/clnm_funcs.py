
from dash import html, dash_table
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto

import numpy as np
import pyreadr
import textwrap

import definitions.layout_styles as styles
import definitions.elements_ids as ids

from definitions.general_funcs import get_label


# -- Frontend --------------------

def display_column2(which_net, position):
    specs = {'t': ['Temporal (within-person) network', ids.TEMP_NET, ids.TEMP_CI_TAB],
             'c': ['Contemporaneous (within-person) network', ids.CONT_NET, ids.CONT_CI_TAB],
             'b': ['Contemporaneous (between-person) network', ids.BETW_NET, ids.BETW_CI_TAB]}

    col = dbc.Col(width=4,
                  children=[
                      html.H5(specs[which_net][0], style=styles.SUB_TITLE2),
                      dbc.Stack(gap=2,
                                children=[
                                    html.Div(children=[cyto.Cytoscape(id=specs[which_net][1],
                                                                     layout={'name': 'preset'},
                                                                     style={'width': '30%', 'height': '60%',
                                                                            'position': 'absolute', 'left': position,
                                                                            'top': styles.CLNM_V_POS},
                                                                     minZoom=1, maxZoom=1,  # disable user zooming
                                                                     elements=make_net2(which_net),
                                                                     stylesheet=make_style_net2(which_net))
                                                      ]),
                                    html.Div(dbc.Accordion(start_collapsed=True,
                                                           children=[
                                                               dbc.AccordionItem(
                                                                   title='Inspect centrality indices',
                                                                   children=[
                                                                       dash_table.DataTable(id=specs[which_net][2],
                                                                           data=make_table2(which_net).to_dict(
                                                                               'records'),
                                                                           columns=[{'name': i, 'id': i} for i in
                                                                                    make_table2(which_net).columns],
                                                                           sort_action='custom', sort_mode='single',
                                                                           sort_by=[],
                                                                           fixed_columns={'headers': True, 'data': 0},
                                                                           # Fix node name column
                                                                           style_header={'fontWeight': 'bold'},
                                                                           style_cell={'fontSize': 18,
                                                                                       'font-family': 'sans-serif'},
                                                                           style_cell_conditional=[
                                                                               {'if': {'column_id': 'Node'},
                                                                                'width': '300px'}],
                                                                           style_data={'whiteSpace': 'normal',
                                                                                       'height': 'auto',
                                                                                       'lineHeight': '20px',
                                                                                       'minWidth': '100px',
                                                                                       'width': '100px',
                                                                                       'maxWidth': '100px'},
                                                                           style_table={'overflowX': 'auto',
                                                                                        'minWidth': '100%'})
                                                                   ], style=styles.TEXT)]))
                                ])
                  ])
    return col


# -- Backend ---------------------


def read_res2(which_net, path='./assets/results/mod2/'):
    """Input: temporal / contemporaneous / between person network input.
       Open the .RData file created by Rscript 2.CLNM. This contains the following elements:
       - fit: fit measures + number of observations the network is based on.
       - layout: the spring graphical disposition computed by `qgraph`
       - t_net/c_net/b_net: dataframe with all edge weights
       - t_cent/c_cent/b_cent:: 95% confidence intervals for those weights
       Use: fit, lay, wm, ci = read_res2(t')
    """
    res = pyreadr.read_r(f'{path}unpruned_norm_info.RData')

    fit = res['fit']

    # layout computed by qgraph (spring algorithm)
    lay = res['layout'].rename(columns={0: 'x_og', 1: 'y_og'})
    # rescale to pixels
    lay['x'] = np.interp(lay['x_og'], (lay['x_og'].min(), lay['x_og'].max()), (0, styles.CLNM_WIDTH))
    lay['y'] = np.interp(lay['y_og'], (lay['y_og'].min(), lay['y_og'].max()), (0, styles.CLNM_WIDTH))
    # add class
    lay['class'] = ['dep' if t else 'cmr' for t in lay.index.str.contains('DEP')]

    nw = res[f'{which_net}_net']
    nw['from'] = nw['from'].map({x: lay.index[x - 1] for x in nw['from']})
    nw['to'] = nw['to'].map({x: lay.index[x - 1] for x in nw['to']})
    nw['dir'] = ['neg' if x < 0 else 'pos' for x in nw.weight]

    ci = res[f'{which_net}_cent']

    return fit, lay, nw, ci


def make_net2(which_net):
    """Input: network type (temporal, contemporaneous or between person).
       Creates a network structure.
    """
    fit, lay, nw, _ = read_res2(which_net)  # read in data

    # already trimmed estimates
    nodes = [{'data': {'id': node, 'label': '\n'.join(textwrap.wrap(get_label(node), width=20))},
              'classes': lay.loc[node, 'class'],
              'position': {'x': lay.loc[node, 'x'], 'y': lay.loc[node, 'y']}}
             for node in lay.index]

    edges = [{'data': {'source': a, 'target': b, 'weight': round(abs(w) * 10, 3), 'width': round(abs(w) * 10, 2)},
              'classes': c}
             for a, b, w, c in nw.itertuples(index=False)]

    network = nodes + edges

    return network


def make_style_net2(which_net):
    sn2 = [{'selector': 'node', 'style': {'height': styles.CLNM_NODE_SIZE, 'width': styles.CLNM_NODE_SIZE,
                                          'label': 'data(label)', 'text-wrap': 'wrap',
                                          'font-size': styles.CLNM_NODE_LABEL}},
           # Edge opacity and width
           {'selector': 'edge', 'style': {'opacity': 'data(weight)', 'width': 'data(width)'}},
           # Color nodes by group
           {'selector': '.dep', 'style': {'background-color': 'lightblue'}},
           {'selector': '.cmr', 'style': {'background-color': 'pink'}},
           # Color edges by positive/negative weights
           {'selector': '.neg', 'style': {'line-color': 'red', 'target-arrow-color': 'red'}},
           {'selector': '.pos', 'style': {'line-color': 'blue', 'target-arrow-color': 'blue'}}]

    if which_net == 't':
        sn2 = sn2 + [
            {'selector': 'edge', 'style': {'curve-style': 'straight', 'target-arrow-shape': 'vee', 'arrow-scale': .8}}]

    return sn2


def make_table2(which_net):
    """Input: network type (temporal, contemporaneous or between person).
       Creates the table with centrality indices tab.
    """
    ci = read_res2(which_net)[3]  # read in data

    # Centrality indices tab
    ci_tab = ci.round(2)
    ci_tab.insert(0, 'Node', [get_label(name) for name in ci_tab.index])  # Replace node name with its label

    return ci_tab
