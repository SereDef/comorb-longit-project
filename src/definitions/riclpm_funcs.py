import pandas as pd
import pyreadr
import textwrap

import definitions.layout_styles as styles
from definitions.general_funcs import badge_it, get_label

# -- Backend ---------------------

def read_res_riclpm(depname, cmrname, params='stat', sex='', path='./assets/results/mod1/'):
    """Input: names of the depression report (sDEP = self or mDEP = parental reports) and cardio-metabolic risk (CMR)
       marker, + stationarity assumptions for long term parameters ('free' vs. 'stat').
       Open the .RData file created by Rscript 1.RICLPM (one for dep-cmr marker pair).
       This contains the following elements:
       - fm: fit measures
       - estimates: unstandardized estimates (+ robust SE, pvalues and CIs)
    """
    sexdict = {'m': '_sex_1', 'f': '_sex_2', '': ''}

    res = pyreadr.read_r(f'{path}ri_p{params}/{depname}_{cmrname}{sexdict[sex]}.RData')

    cors = res['c']
    fitm = res['fm'].T

    esti = res['estimates']
    esti['sign'] = [0 if low <= 0 <= upp else 1 for low, upp in
                    zip(esti['ci.lower'], esti['ci.upper'])]  # based on CI because pvalue is sometimes NA
    # Add standardized estimates
    esti.insert(5, 'std_est', res['stad_esti'])
    # Make covariance easier to query
    if params == 'free':
        esti['label'].replace('cor1', 'rcov1', inplace=True)

    return fitm, esti, cors


def make_net_riclpm(depname, cmrname, params='stat', sex='', net_width=styles.CLPM_WIDTH):
    """Input: names of the depression report (sDEP = self or mDEP = parental reports) and cardio-metabolic
       risk (CMR) marker; model stationarity and size of the graph.
       Creates a list of dictionaries to use as input for the elements arg of cytoscape graph.
       This only creates the core structure, style parameters are defined in style_sheet1.
    """
    # read data
    _, esti, cor = read_res_riclpm(depname, cmrname, params, sex)

    # load summary dataframe
    summ = pyreadr.read_r(f'./assets/results/mod1/g_lstat_pstat/{depname}_{cmrname}.RData')['dat_summ']

    # Extract the estimated parameters from the result files (returns estimates and significance [= '*' or '']
    def extr_est(name):
        if params == 'free':
            name = name.replace('cor', 'rcov')

        e, l, u, s = esti.loc[esti.label.str.contains(name)][
            ['std_est', 'ci.lower', 'ci.upper', 'sign']].round(2).iloc[0]
        return e, l, u, s

    # Ready to draw
    nt = summ.shape[1] // 2  # Number of time points

    # Get time points of measurement from summary
    timedic = dict(
        zip([f'dep{i + 1}' if 'DEP' in n else f'cmr{i - summ.shape[1] // 2 + 1}' for i, n in enumerate(summ.columns)],
            [float(i.split('_')[-1][:-1]) for i in summ.columns]))

    pos_top = 30
    pos_bot = 550  # Vertical coordinates (in pixel)
    rshift = 90
    obs_pos = 90
    imp_pos = 150  # Distances from left edge and top (in pixels)
    vs = ['dep', 'cmr']

    elm = list()  # Initialize

    # Variable names
    elm.append({'data': {'id': f'name_dep', 'label': 'Depression\nscore'}, 'classes': 'notes',
                'position': {'x': 120, 'y': pos_top + obs_pos}})
    elm.append({'data': {'id': f'name_cmr', 'label': '\n'.join(textwrap.wrap(get_label(cmrname), width=10))},
                'classes': 'notes',
                'position': {'x': 120, 'y': pos_bot - obs_pos}})

    for ri, pos in enumerate([pos_top, pos_bot]):
        # Eta factors nodes
        elm.append({'data': {'id': f'ri{ri}', 'label': 'RI'},
                    'classes': 'latent',
                    'position': {'x': net_width / 2 + rshift,
                                 'y': pos}
                    })

    #  Eta factors correlations
    elm.append({'data': {'source': 'ri0',
                         'target': 'ri1',
                         'firstname': 'covRI',
                         'label': '%.2f' % extr_est('covRI')[0]}
                })

    for i in range(1, nt + 1):

        elm.append({'data': {'source': f'imp_dep{i}', 'target': f'imp_cmr{i}', 'firstname': 'imp_cov',
                             'label': '%.2f' % extr_est(f'cor{i}')[0]}})

        for ri, v in enumerate(vs):

            # define vertical position
            p = [pos_top + obs_pos, pos_top + imp_pos] if v == 'dep' else [pos_bot - obs_pos, pos_bot - imp_pos]

            # Observed variables and impulses
            elm.extend([
                # Observed variables
                {'data': {'id': f'{v}{i}',
                          'firstname': f'{v.upper()}\n{timedic[v + str(i)]} yrs',
                          'label': summ.columns[(i - 1) + (nt * ri)]},
                 'classes': 'observed',
                 'position': {'x': ((net_width / nt) * i) - (net_width / nt) / 2 + rshift, 'y': p[0]}},
                # Impulses
                {'data': {'id': f'imp_{v}{i}',
                          'label': f'impulse {i}'},
                 'classes': 'latent',
                 'position': {'x': ((net_width / nt) * i) - (net_width / nt) / 2 + rshift, 'y': p[1]}},

                # impulses links
                {'data': {'source': f'imp_{v}{i}',
                          'target': f'{v}{i}',
                          'firstname': 'imp_link'}},
                # lambdas
                {'data': {'source': f'ri{ri}',
                          'target': f'{v}{i}',
                          'firstname': 'lambda',
                          'label': '1'}}
            ])

            if i < nt:
                otherv = abs(ri - 1)
                # time adjustment
                t1 = timedic[f'{v}{i + 1}'] - timedic[f'{v}{i}']
                t2 = timedic[f'{v}{i + 1}'] - timedic[f'{vs[otherv]}{i}']

                def calc_weight(term, vname, index):
                    if term == 'AR':
                        t = timedic[f'{vname}{index + 1}'] - timedic[f'{vname}{index}']
                    elif term == 'CL':
                        t = timedic[f'{vname}{index + 1}'] - timedic[f'{vs[otherv]}{index}']

                    if params == 'free':
                        t = 1 / t  # divide by time if not stationary

                    weight = '%.2f\n[%.2f; %.2f]' % (extr_est(f'{term}_{vname}{index}')[0] * t,
                                                     extr_est(f'{term}_{vname}{index}')[1] * t,
                                                     extr_est(f'{term}_{vname}{index}')[2] * t)
                    return weight

                # AR and CL terms
                elm.extend([{'data': {'source': f'imp_{v}{i}', 'target': f'imp_{v}{i + 1}',
                                      'weight': calc_weight('AR', v, i),
                                      'sign': extr_est(f'AR_{v}{i}')[3],
                                      'label': f'AR{i}',
                                      'firstname': 'direct'}},

                            {'data': {'source': f'imp_{vs[otherv]}{i}', 'target': f'imp_{v}{i + 1}',
                                      'weight': calc_weight('CL', v, i),
                                      'sign': extr_est(f'CL_{v}{i}')[3],
                                      'label': f'CL{i}',
                                      'firstname': 'direct'}}
                            ])

    return elm


# Define the stile of the graph
width = styles.CLPM_WIDTH

style_net_riclpm = [
    # Noting down names
    {'selector': '.notes',
     'style': {'background-color': 'white', 'font-size': styles.CLPM_NODE_LABEL,
               'content': 'data(label)', 'color': 'k', 'font-weight': 'bold', 'text-halign': 'left',
               'text-valign': 'center', 'text-wrap': 'wrap'}},

    # Nodes - shape & color
    {'selector': '.observed',
     'style': {'shape': 'rectangle', 'height': styles.CLPM_NODE_SIZE * 2, 'width': styles.CLPM_NODE_SIZE * 4,
               'border-width': 2, 'background-color': 'white', 'border-color': 'k',
               'content': 'data(firstname)', 'color': 'grey', 'text-halign': 'center', 'text-valign': 'center',
               'text-wrap': 'wrap'}},

    {'selector': '.latent',
     'style': {'shape': 'round', 'height': styles.CLPM_NODE_SIZE, 'width': styles.CLPM_NODE_SIZE,
               'border-width': 1, 'background-color': 'white', 'border-color': 'silver'}},

    # Edges
    {'selector': 'edge[firstname *= "direct"]',  # directed paths
     'style': {'curve-style': 'straight', 'target-arrow-shape': 'vee', 'width': 3, 'arrow-scale': 1.2}},

    {'selector': 'edge[firstname *= "imp_link"]',  # impulses links
     'style': {'curve-style': 'straight', 'target-arrow-shape': 'vee', 'width': 1, 'arrow-scale': .8}},

    # Correlations
    {'selector': 'edge[firstname *= "covRI"]',
     'style': {'curve-style': 'unbundled-bezier', 'target-arrow-shape': 'vee', 'source-arrow-shape': 'vee', 'width': 1,
               'control-point-distances': [-width * .45, -width * .52, -width * .53, -width * .53, -width * .45],
               'control-point-weights': [0.01, 0.20, 0.5, 0.80, 0.99],
               'label': 'data(label)', 'font-size': styles.CLPM_EDGE_LABEL, 'text-background-color': 'silver',
               'text-background-opacity': .7}},

    {'selector': 'edge[firstname *= "imp_cov"]',
     'style': {'curve-style': 'unbundled-bezier', 'target-arrow-shape': 'vee', 'source-arrow-shape': 'vee', 'width': 1,
               'label': 'data(label)', 'font-size': styles.CLPM_EDGE_LABEL, 'text-background-color': 'silver',
               'text-background-opacity': .7}},

    # Lambdas
    {'selector': 'edge[firstname *= "lambda"]',
     'style': {'curve-style': 'straight', 'target-arrow-shape': 'vee', 'width': 1, 'arrow-scale': .8,
               'label': 'data(label)', 'font-size': styles.CLPM_EDGE_LABEL, 'text-background-color': 'silver',
               'text-background-opacity': .7}},

    # Dashed lines
    # {'selector':'edge[weight < 0.01][weight > -0.01]', 'style':{'line-style':'dashed'}},
    # {'selector':'edge[sign < 1]', 'style':{'line-style':'dashed'}},
]
# Set the color of each edge type and the distance between source and label displaying its weight (to avoid overlapping)
d = {'AR': ['red', styles.CLPM_WIDTH/14],
     'CL': ['green', styles.CLPM_WIDTH/14]}

for c in d.keys():
    style_net_riclpm.extend([{'selector': f'[label *= "{c}"][sign > 0]',
                              'style': {'line-color': d[c][0],
                                        'target-arrow-color': d[c][0],
                                        'source-label': 'data(weight)',
                                        'source-text-offset': d[c][1],
                                        'font-size': styles.CLPM_EDGE_LABEL * 1.25,
                                        'text-wrap': 'wrap',
                                        'font-weight': 'bold',
                                        'text-background-color': d[c][0],
                                        'text-background-opacity': .5}},
                             {'selector': f'[label *= "{c}"][sign < 1]',
                              'style': {'line-style': 'dashed',
                                        'width': 1.5,
                                        'line-color': d[c][0],
                                        'target-arrow-color': d[c][0],
                                        'source-label': 'data(weight)',
                                        'source-text-offset': d[c][1],
                                        'font-size': styles.CLPM_EDGE_LABEL,
                                        'text-wrap': 'wrap',
                                        'text-background-color': d[c][0],
                                        'text-background-opacity': .5}}
                       ])


def make_table_riclpm(depname, cmrname, params='stat', sex=''):
    """Input: names of the depression report (sDEP = self or mDEP = parental reports) and cardio-metabolic
       risk (CMR) marker; model structure.
       Extracts the fit measures for the specified model and stores into a table to be displayed next to the graph.
    """
    fitm = read_res_riclpm(depname, cmrname, params, sex)[0]

    dt = pd.DataFrame(
        fitm[['ntotal','npar', 'df', 'chisq', 'pvalue', 'cfi', 'tli', 'rmsea', 'srmr', 'aic', 'bic']]).T
    dt = dt.rename(columns={'lavaan::fitmeasures(m)': ' '}).round(3).astype('object')
    dt.insert(loc=0, column='Fit measures',
              value=['Sample size', 'Number of parameters', 'Degrees of freedom', '\u03C7\u00b2', 'P-value (\u03C7\u00b2)',
                     'CFI', 'TLI', 'RMSEA', 'SRMR', 'AIC', 'BIC'])

    def format_row(x, form):
        try:
            return form.format(x)
        except ValueError:
            return x

    format_dict = {
        '{:d}': ['ntotal', 'npar', 'df'],  # integers
        '{0:.1f}': ['aic', 'bic'],  # 1 decimal point
        '{0:.2f}': ['chisq', 'cfi', 'tli', 'rmsea', 'srmr'],  # 2 decimal points
        '{0:.3f}': ['pvalue'],  # 3 decimal points
    }
    for fmt in format_dict:
        dt.loc[format_dict[fmt]] = dt.loc[format_dict[fmt]].map(format_row, form=fmt)

    # dt.reset_index(drop=True, inplace=True)

    return dt
