import numpy as np
import os
import pyreadr
import textwrap

import definitions.layout_styles as styles
from definitions.general_funcs import get_label

PATH = './assets/results/mod3/'

# -- Frontend --------------------

times = [file.split('_')[1][:-7] for file in os.listdir(PATH) if file.endswith('.RData')]
# times.sort()  # just for readability, not necessary

# create marks dict to use in slider, also provides a map for make_net3
time_marks3 = dict()
for t in times:
    where = round(np.mean([float(n) for n in t.split('-')]), 1)
    time_marks3[where] = {'label': f'\n{t} years',
                          'style': {'transform': 'rotate(45deg)', 'whitespace': 'nowrap',
                                    'font-size': styles.CSNM_NODE_LABEL*.8}}  # 'color': '#f50'


# -- Backend ---------------------


def read_res3(time, path=PATH):
    """Input: time point of interest.
       Open the .RData file created by Rscript 2.CLNM (one for each time point). This contains the following elements:
       - wm: dataframe with all edge weights
       - ci: 95% confidence intervals for those weights
       - fit: fit measures + number of observations the network is based on.
       - layout: the spring graphical disposition computed by `qgraph`
       Use: wm, ci, fit, lay = read_res3('9.8y-9.8y')
    """
    res = pyreadr.read_r(f'{path}crosnet_{time}y.RData')
    # weight matrix
    wm = res['wm']
    wm['link'] = wm.index
    wm[['a', 'b']] = wm.link.str.split(' ', expand=True)
    wm = wm.loc[wm.a != wm.b, ]  # remove links to between an edge and itself
    wm = wm.reset_index()[['a', 'b', 'V1']].rename(columns={'a': 'node1', 'b': 'node2', 'V1': 'weight'})
    wm['dir'] = ['neg' if x < 0 else 'pos' for x in wm.weight]
    # centrality indices
    ci = res['ci']
    ci['class'] = ['dep' if i else 'cmr' for i in ci.node.str.contains('DEP')]
    # fit measures and number of observations
    fit = res['fit'].T.round(3)
    # layout computed by qgraph (spring algorithm)
    lay = res['layout'].rename(columns={0: 'x_og', 1: 'y_og'})
    # rescale to pixels
    lay['x'] = np.interp(lay['x_og'], (lay['x_og'].min(), lay['x_og'].max()), (0, styles.CSNM_WIDTH))
    lay['y'] = np.interp(lay['y_og'], (lay['y_og'].min(), lay['y_og'].max()), (0, styles.CSNM_WIDTH))

    return wm, ci, fit, lay


def make_net3(timepoint):
    """Input: time point of interest.
       Creates a network structure and the table with centrality indices.
    """
    lab = time_marks3[timepoint]['label'][1:-6]  # get label and remove '\n' in the beginning and ' years' at the end

    wm, ci, _, lay = read_res3(lab)  # read in data

    # tim estimates?
    wm_trim = wm.loc[abs(wm.weight) > 0.01, ].reset_index(drop=True)

    nodes = [{'data': {'id': node, 'label': '\n'.join(textwrap.wrap(get_label(node), width=20))}, 'classes': group,
              'position': {'x': lay.loc[node, 'x'], 'y': lay.loc[node, 'y']}}
             for node, group in ci[['node', 'class']].itertuples(index=False)]

    edges = [{'data': {'source': a, 'target': b, 'weight': w, 'width': round(abs(w) * 20, 2)}, 'classes': c}
             for a, b, w, c in wm_trim.itertuples(index=False)]

    network = nodes + edges

    # Centrality indices tab
    ci_tab = ci.round(2).drop(columns='class')
    ci_tab.insert(1, 'Node', [get_label(name) for name in ci_tab['node']])  # Replace node name with its label
    ci_tab = ci_tab.drop(columns='node')

    return network, ci_tab


style_net3 = [{'selector': 'node', 'style': {'height': styles.CSNM_NODE_SIZE, 'width': styles.CSNM_NODE_SIZE,
                                             'label': 'data(label)', 'text-wrap': 'wrap',
                                             'font-size': styles.CSNM_NODE_LABEL}},
              # Edge opacity and width
              {'selector': 'edge', 'style': {'opacity': 'data(weight)', 'width': 'data(width)'}},
              # Color nodes by group
              {'selector': '.dep', 'style': {'background-color': 'lightblue'}},
              {'selector': '.cmr', 'style': {'background-color': 'pink'}},
              # Color edges by positive/negative weights
              {'selector': '.neg', 'style': {'line-color': 'red'}},
              {'selector': '.pos', 'style': {'line-color': 'blue'}}]
