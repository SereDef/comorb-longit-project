from dash import html, dcc
import dash_bootstrap_components as dbc

import pandas as pd
import pyreadr
import textwrap

import definitions.layout_styles as styles

from definitions.general_funcs import badge_it, get_label


# -- Frontend --------------------
def param_checklist(depname, cmrname, p='lt', best=False, badgefont=styles.TEXT['font-size']):
    """Defines the checkboxes used to determine model structure (what parameters to estimate)."""

    pref = '' if p == 'lt' else 'ma'
    disable_ar = True if p == 'lt' else False  # long-term AR terms are always included
    cols = ['crimson', 'green'] if p == 'lt' else ['orange', 'lightblue']
    position = 'left' if p == 'lt' else 'right'

    if best:
        val = best_fit1(depname, cmrname, list1=p)  # outputs param list of the best fitting model
    else:
        val = ['ltCL_dep', 'ltCL_cmr', 'ltAR_dep', 'ltAR_cmr'] if p == 'lt' else []  # classic CLPM

    return html.Div(style={'width': '50%', 'height': '65%', 'float': position},
                    children=[dcc.Checklist(id=f'{p}-checklist',
                                            options=[
                                                {'label': html.Span([
                                                    badge_it(f'{pref}AR', cols[0], badgefont),
                                                    ' depression']),
                                                 'value': f'{p}AR_dep', 'disabled':disable_ar},
                                                {'label': html.Span(
                                                    [badge_it(f'{pref}AR', cols[0], badgefont),
                                                     ' cardio-metabolic risk']),
                                                 'value': f'{p}AR_cmr', 'disabled':disable_ar},
                                                {'label': html.Span(
                                                    [badge_it(f'{pref}CL', cols[1], badgefont),
                                                     ' depression \u290F cardio-metab.']),
                                                 'value': f'{p}CL_dep'},
                                                {'label': html.Span(
                                                    [badge_it(f'{pref}CL', cols[1], badgefont),
                                                     ' cardio-metab. \u290F depression']),
                                                 'value': f'{p}CL_cmr'}],
                                            value=val,
                                            style=styles.TEXT,
                                            inputStyle={'margin-left': styles.MARGIN_CHECKLIST,
                                                        'margin-right': styles.MARGIN_CHECKLIST},
                                            labelStyle={'display': 'block'})])


def make_button(label, id_name, color, margin='15px', fs=styles.TEXT['font-size']):
    return dbc.Button(label, id=id_name, color='secondary', n_clicks=0,
                      style={'font-size': fs, 'font-weight': 'bold', 'background-color': color,
                             'padding': '4px 10px', 'margin-left': margin})


# -- Backend ---------------------

# Matrix with parameter combinations (8x57) used to fit the models in Rscript 1.
model_structure = pd.read_csv('./assets/model_structure.csv').set_index('Unnamed: 0')


def read_res_gclpm(depname, cmrname, lambdas='free', params='stat', path='./assets/results/mod1/'):
    """Input: names of the depression report (sDEP = self or mDEP = parental reports) and cardio-metabolic risk (CMR)
       marker, + stationarity assumptions for lambdas and long term parameters ('free' vs. 'stat').
       Open the .RData file created by Rscript 1.GCLPM (one for dep-cmr marker pair).
       This contains the following elements:
       - summ: a summary dataframe (with information about which marker is used, the timepoints included
         + mean ranges and number of observations)
       - fit_meas: fit measures for every parameter combination, when the model converged.
       - estimates: (unstandardized) estimates (+ robust SE, pvalues and CIs) for every parameter combination,
         when the model converged.
       - failed: list of models that did not converge, with corresponding error or warning message.
       Use: summ, fitm, esti, fail = read_res1('sDEP','FMI') -OR- summ = read_res1('sDEP','BMI')[0]
    """
    res = pyreadr.read_r(f'{path}g_l{lambdas}_p{params}/{depname}_{cmrname}.RData')

    summ = res['dat_summ']
    fitm = res['fit_meas'].T
    esti = res['estimates'].set_index('rep(f, nrow(es))')
    esti['sign'] = [0 if low <= 0 <= upp else 1 for low, upp in
                    zip(esti['ci.lower'], esti['ci.upper'])]  # based on CI because pvalue is sometimes NA
    fail = list(res['failed'].index)  # TODO: report warning message ...

    return summ, fitm, esti, fail


def best_fit_gclpm(depname, cmrname, lambdas, params, list1=None):
    """Input: names of the depression report (sDEP = self or mDEP = parental reports) and cardio-metabolic risk (CMR)
       marker.
       By default, returns a dataframe with the model name (indicating the excluded parameters) and structure (0,1...).
       When list1 is provided, returns the list of long-term (list1='lt') or short-term (list1='ma') parameters
       estimated in the model.
    """
    fitm = read_res1(depname, cmrname, lambdas, params)[1]

    # Best fitting model (lowest BIC)
    mod = fitm.index[fitm.bic == fitm.bic.min()][0]
    if list1:  # Return list of parameters estimated in the model
        return list(model_structure.index[(model_structure[mod] > 0) & (model_structure.index.str.contains(list1))])
    else:  # Return a dataframe with its name and model structure
        return pd.DataFrame(model_structure[mod])


def make_net_gclpm(depname, cmrname, lambdas='free', params='stat', which_model='maAR_dep-maCL_dep-maAR_cmr-maCL_cmr',
              net_width=styles.CLPM_WIDTH):
    """Input: names of the depression report (sDEP = self or mDEP = parental reports) and cardio-metabolic
       risk (CMR) marker; model structure and size of the graph.
       Creates a list of dictionaries to use as input for the elements arg of cytoscape graph.
       This only creates the core structure, style parameters are defined in style_sheet1.
    """
    # read data
    summ, _, esti, fail = read_res_gclpm(depname, cmrname, lambdas, params)

    # Best fitting model
    if which_model == 'best':
        modstr = best_fit_gclpm(depname, cmrname, lambdas, params).T
    elif which_model in fail:
        return 'fail'
    else:
        modstr = pd.DataFrame(model_structure[which_model]).T

    # Extract the estimated parameters from the result files (returns estimates and significance [= '*' or '']
    def extr_est(name=None, which=None, lamb=False, eta_corr=False, imp_corr=False, model_output=esti):
        df = model_output.loc[modstr.index[0]]
        if lamb:
            e, l, u, s = df.loc[(df.lhs == f'eta_{name}') & (df.op == '=~'), ['est', 'ci.lower', 'ci.upper', 'sign'
                                                                              ]].round(2).iloc[which]
        elif eta_corr:
            e, l, u, s = df.loc[(df.lhs == 'eta_dep') & (df.rhs == 'eta_cmr'), ['est', 'ci.lower', 'ci.upper', 'sign'
                                                                                ]].round(2).iloc[0]
        elif imp_corr:
            e, l, u, s = df.loc[df.label.str.contains('comv'), ['est', 'ci.lower', 'ci.upper', 'sign']].round(2).iloc[which]
        else:
            e, l, u, s = df.loc[df.label.str.contains(name)].iloc[::-1][['est', 'ci.lower', 'ci.upper', 'sign'
                                                                         ]].round(2).iloc[which]
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
    imp_pos = 190  # Distances from left edge and top (in pixels)
    vs = ['dep', 'cmr']

    elm = list()  # Initialize

    # Variable names
    elm.append({'data': {'id': f'name_dep', 'label': 'Depression\nscore'}, 'classes': 'notes',
                'position': {'x': 120, 'y': pos_top + obs_pos}})
    elm.append({'data': {'id': f'name_cmr', 'label': '\n'.join(textwrap.wrap(get_label(cmrname), width=10))},
                'classes': 'notes',
                'position': {'x': 120, 'y': pos_bot - obs_pos}})

    for eta, pos in enumerate([pos_top, pos_bot]):
        # Eta factors nodes
        elm.append({'data': {'id': f'eta{eta}', 'label': 'Eta'}, 'classes': 'latent',
                    'position': {'x': net_width / 2 + rshift, 'y': pos}})

    #  Eta factors correlations
    elm.append({'data': {'source': 'eta0', 'target': 'eta1', 'firstname': 'eta_corr',
                         'label': '%.2f' % extr_est(eta_corr=True)[0]}})

    for i in range(1, nt + 1):

        imp_ = 'imp_' if i > 1 else ''
        elm.append({'data': {'source': f'{imp_}dep{i}', 'target': f'{imp_}cmr{i}', 'firstname': 'imp_corr',
                             'label': '%.2f' % extr_est(which=i - 1, imp_corr=True)[0]}})

        for eta, v in enumerate(vs):

            # ===== Other nodes =====
            # define vertical position
            p = [pos_top + obs_pos, pos_top + imp_pos] if v == 'dep' else [pos_bot - obs_pos, pos_bot - imp_pos]

            # Observed variables
            elm.append(
                {'data': {'id': f'{v}{i}', 'firstname': f'{v.upper()}\n{timedic[v + str(i)]} yrs',
                          'label': summ.columns[(i - 1) + (nt * eta)]},
                 'classes': 'observed',
                 'position': {'x': ((net_width / nt) * i) - (net_width / nt) / 2 + rshift, 'y': p[0]}})

            # Impulses
            if i > 1:
                elm.append(
                    {'data': {'id': f'imp_{v}{i}', 'label': f'impulse {i}'},
                     'classes': 'latent',
                     'position': {'x': ((net_width / nt) * i) - (net_width / nt) / 2 + rshift, 'y': p[1]}})
                # impulses links
                elm.append({'data': {'source': f'imp_{v}{i}', 'target': f'{v}{i}', 'firstname': 'imp_link'}})

            # ===== Edges =====
            # lambdas
            elm.append({'data': {'source': f'eta{eta}', 'target': f'{v}{i}', 'firstname': 'lambda',
                                 'label': '%.2f' % extr_est(f'{v}', i - 1, lamb=True)[0]}})

            if i < nt:
                otherv = abs(eta - 1)
                # time adjustment
                t1 = timedic[f'{v}{i + 1}'] - timedic[f'{v}{i}']
                t2 = timedic[f'{v}{i + 1}'] - timedic[f'{vs[otherv]}{i}']
                # AR and CL terms
                if modstr[f'ltAR_{v}'].iloc[0]:
                    elm.append({'data': {'source': f'{v}{i}',
                                         'target': f'{v}{i + 1}',
                                         'weight': '%.2f\n[%.2f; %.2f]' % (extr_est(f'^AR_{v}', i - 1)[0] * t1,
                                                                           extr_est(f'^AR_{v}', i - 1)[1] * t1,
                                                                           extr_est(f'^AR_{v}', i - 1)[2] * t1),
                                         'sign': extr_est(f'^AR_{v}', i - 1)[3],
                                         'label': f'AR{i}', 'firstname': 'direct'}})
                if modstr[f'ltCL_{v}'].iloc[0]:
                    elm.append({'data': {'source': f'{vs[otherv]}{i}',
                                         'target': f'{v}{i + 1}',
                                         'weight': '%.2f\n[%.2f; %.2f]' % (extr_est(f'^CL_{v}', i - 1)[0] * t2,
                                                                           extr_est(f'^CL_{v}', i - 1)[1] * t2,
                                                                           extr_est(f'^CL_{v}', i - 1)[2] * t2),
                                         'sign': extr_est(f'^CL_{v}', i - 1)[3],
                                         'label': f'CL{i}', 'firstname': 'direct'}})
                # ma AR and maCL terms
                if i > 1:
                    if modstr[f'maAR_{v}'].iloc[0]:
                        elm.append({'data': {'source': f'imp_{v}{i}',
                                             'target': f'{v}{i + 1}',
                                             'weight': '%.2f\n[%.2f; %.2f]' % (extr_est(f'^maAR_{v}', i - 2)[0] * t1,
                                                                               extr_est(f'^maAR_{v}', i - 2)[1] * t1,
                                                                               extr_est(f'^maAR_{v}', i - 2)[2] * t1),
                                             'sign': extr_est(f'^maAR_{v}', i - 2)[3],
                                             'label': f'maAR{i}', 'firstname': 'direct'}})
                    if modstr[f'maCL_{v}'].iloc[0]:
                        elm.append({'data': {'source': f'imp_{vs[otherv]}{i}',
                                             'target': f'{v}{i + 1}',
                                             'weight': '%.2f\n[%.2f; %.2f]' % (extr_est(f'^maCL_{v}', i - 2)[0] * t2,
                                                                               extr_est(f'^maCL_{v}', i - 2)[1] * t2,
                                                                               extr_est(f'^maCL_{v}', i - 2)[2] * t2),
                                             'sign': extr_est(f'^maCL_{v}', i - 2)[3],
                                             'label': f'maCL{i}', 'firstname': 'direct'}})
    return elm


# Define the stile of the graph
width = styles.CLPM_WIDTH

style_net_gclpm = [
    # Noting down names
    {'selector': '.notes',
     'style': {'background-color': 'white', 'font-size': styles.CLPM_NODE_LABEL,
               'content': 'data(label)', 'color': 'k', 'font-weight': 'bold', 'text-halign': 'left',
               'text-valign': 'center', 'text-wrap': 'wrap'}},

    # Nodes - shape & color
    {'selector': '.observed',
     'style': {'shape': 'rectangle', 'height': styles.CLPM_NODE_SIZE*2, 'width': styles.CLPM_NODE_SIZE*4,
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
    {'selector': 'edge[firstname *= "eta_corr"]',
     'style': {'curve-style': 'unbundled-bezier', 'target-arrow-shape': 'vee', 'source-arrow-shape': 'vee', 'width': 1,
               'control-point-distances': [-width*.45, -width*.52, -width*.53, -width*.53, -width*.45],
               'control-point-weights': [0.01, 0.20, 0.5, 0.80, 0.99],
               'label': 'data(label)', 'font-size': styles.CLPM_EDGE_LABEL, 'text-background-color': 'silver',
               'text-background-opacity': .7}},

    {'selector': 'edge[firstname *= "imp_corr"]',
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
d = {'AR': ['red', styles.CLPM_WIDTH/25],
     'CL': ['green', styles.CLPM_WIDTH/3.8],
     'maAR': ['orange', styles.CLPM_WIDTH/12],
     'maCL': ['lightblue', styles.CLPM_WIDTH/20]}

for c in d.keys():
    style_net_gclpm.extend([{'selector': f'[label *= "{c}"][sign > 0]',
                        'style': {'line-color': d[c][0],
                                  'target-arrow-color': d[c][0],
                                  'source-label': 'data(weight)',
                                  'source-text-offset': d[c][1],
                                  'font-size': styles.CLPM_EDGE_LABEL*1.25,
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
                                  'text-wrap': 'wrap',
                                  'source-text-offset': d[c][1],
                                  'font-size': styles.CLPM_EDGE_LABEL,
                                  'text-background-color': d[c][0],
                                  'text-background-opacity': .5}}
                       ])


def make_table_gclpm(depname, cmrname, lambdas='free', params='stat', which_model='maAR_dep-maCL_dep-maAR_cmr-maCL_cmr'):
    """Input: names of the depression report (sDEP = self or mDEP = parental reports) and cardio-metabolic
       risk (CMR) marker; model structure.
       Extracts the fit measures for the specified model and stores into a table to be displayed next to the graph.
    """
    fitm = read_res_gclpm(depname, cmrname, lambdas, params)[1]

    if which_model == 'best':
        which_model = fitm.index[fitm.bic == fitm.bic.min()][0]  # Best fitting model (lowest BIC)

    dt = pd.DataFrame(
        fitm.loc[which_model, ['ntotal','npar', 'df', 'chisq', 'pvalue', 'cfi', 'tli', 'rmsea', 'srmr', 'aic', 'bic', ]])
    dt = dt.rename(columns={which_model: ' '}).round(3).astype('object')
    dt.insert(loc=0, column='Fit measures',
              value=['Sample size','Number of parameters', 'Degrees of freedom', '\u03C7\u00b2', 'P-value (\u03C7\u00b2)',
                     'CFI', 'TLI', 'RMSEA', 'SRMR', 'AIC', 'BIC'])

    def format_row(x, form):
        try:
            return form.format(x)
        except ValueError:
            return x

    format_dict = {
        '{:d}': ['ntotal','npar', 'df'],  # integers
        '{0:.1f}': ['aic', 'bic'],  # 1 decimal point
        '{0:.2f}': ['chisq', 'cfi', 'tli', 'rmsea', 'srmr'],  # 2 decimal points
        '{0:.3f}': ['pvalue'],  # 3 decimal points
    }
    for fmt in format_dict:
        dt.loc[format_dict[fmt]] = dt.loc[format_dict[fmt]].map(format_row, form=fmt)

    return dt
