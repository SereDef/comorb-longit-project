from dash import html, dcc
import dash_bootstrap_components as dbc

from plotly.subplots import make_subplots
import plotly.graph_objects as go

import pandas as pd
import pyreadr
import textwrap

from definitions.general_funcs import badge_it, get_label


# -- Frontend --------------------
def dep_var_checklist(disable_maternal=False):
    opts = [{'label': 'Self-reported', 'value': 'sDEP', 'disabled': False},
            {'label': 'Maternal report', 'value': 'mDEP', 'disabled': disable_maternal}]
    return opts


def cmr_var_checklist(reporter='s'):
    # Read all fitted model names
    # allfiles = [x.split('.')[0] for x in os.listdir(PATH+'mod1')]
    # moth = [s[5:]+'_' for s in allfiles if "mDEP" in s]

    sortd = ['FMI', 'LMI', 'BMI', 'waist_circ', 'android_fatmass', 'total_fatmass', 'total_leanmass',
             'tot_chol', 'HDL_chol', 'LDL_chol', 'insulin', 'triglyc', 'CRP']

    opts = [{'value': v, 'label': get_label(v) + f' ({v})' if len(v) < 4 else get_label(v), 'disabled': False} for v in
            sortd]

    if reporter == 'm':
        sortd_m = ['FMI', 'LMI', 'BMI', 'waist_circ', 'total_fatmass', 'total_leanmass']
        for i in opts:
            if i['value'] not in sortd_m:
                i['disabled'] = True

    return opts


def param_checklist(depname, cmrname, p='lt', best=False):
    pref = '' if p == 'lt' else 'ma'
    cols = ['crimson', 'green'] if p == 'lt' else ['orange', 'lightblue']
    position = 'left' if p == 'lt' else 'right'

    if best:
        val = best_fit1(depname, cmrname, list1=p)
    else:
        val = ['ltCL_dep', 'ltCL_cmr', 'ltAR_dep', 'ltAR_cmr'] if p == 'lt' else []

    return html.Div(style={'width': '50%', 'height': '65%', 'float': position},
                    children=[dcc.Checklist(id=f'{p}-checklist',
                                            options=[
                                                {'label': html.Span([badge_it(f'{pref}AR', cols[0]), ' depression']),
                                                 'value': f'{p}AR_dep'},
                                                {'label': html.Span(
                                                    [badge_it(f'{pref}AR', cols[0]), ' cardio-metabolic risk']),
                                                 'value': f'{p}AR_cmr'},
                                                {'label': html.Span([badge_it(f'{pref}CL', cols[1]),
                                                                     ' depression \u290F cardio-metab.']),
                                                 'value': f'{p}CL_dep'},
                                                {'label': html.Span([badge_it(f'{pref}CL', cols[1]),
                                                                     ' cardio-metab. \u290F depression']),
                                                 'value': f'{p}CL_cmr'}],
                                            value=val,
                                            inputStyle={'margin-left': '20px', 'margin-right': '20px'},
                                            labelStyle={'display': 'block'})])


def make_button(label, id_name, color, margin='15px'):
    return dbc.Button(label, id=id_name, color='secondary', n_clicks=0,
                      style={'font-weight': 'bold', 'background-color': color,
                             'padding': '4px 10px', 'margin-left': margin})


# -- Backend ---------------------

# Matrix with parameter combinations (8x57) used to fit the models in Rscript 1.
model_structure = pd.read_csv('./assets/model_structure.csv').set_index('Unnamed: 0')


def read_res1(depname, cmrname, path='./assets/results/mod1/'):
    """Input: names of the depression report (sDEP = self or mDEP = parental reports) and cardio-metabolic risk (CMR) marker.
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
    res = pyreadr.read_r(f'{path}{depname}_{cmrname}.RData')

    summ = res['dat_summ']
    fitm = res['fit_meas'].T
    esti = res['estimates'].set_index('rep(f, nrow(es))')
    esti['sign'] = [0 if low <= 0 <= upp else 1 for low, upp in
                    zip(esti['ci.lower'], esti['ci.upper'])]  # based on CI because pvalue is sometimes NA
    fail = list(res['failed'].index)  # TODO: report warning message ...

    return summ, fitm, esti, fail


def best_fit1(depname, cmrname, list1=None):
    """Input: names of the depression report (sDEP = self or mDEP = parental reports) and cardio-metabolic risk (CMR) marker.
       By default, returns a dataframe with the model name (indicating the excluded parameters) and structure (0,1...).
       When list1 is provided, returns the list of long-term (list1='lt') or short-term (list1='ma') parameters
       estimated in the model.
    """
    fitm = read_res1(depname, cmrname)[1]

    # Best fitting model (lowest AIC)
    mod = fitm.index[fitm.aic == fitm.aic.min()][0]
    if list1:  # Return list of parameters estimated in the model
        return list(model_structure.index[(model_structure[mod] > 0) & (model_structure.index.str.contains(list1))])
    else:  # Return a dataframe with its name and model structure
        return pd.DataFrame(model_structure[mod])


def make_plot1(depname, cmrname):
    """Input: names of the depression report (sDEP = self or mDEP = parental reports) and cardio-metabolic risk (CMR) marker.
       Creates an interactive scatterplot figure with time of measurement on the x-axis and values of
       dep and cmr markers on the (double) y-axis. Markers indicate the median [IQR] of each measure. Median values and
       time points are also displayed when hovering over the points.
    """
    # load summary dataframe
    summ = read_res1(depname, cmrname)[0]

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

    # Group "cross-sectional" timepoints using grey background
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


def make_net1(depname, cmrname, which_model='maCL_dep-maCL_cmr-maAR_dep-maAR_cmr', width=1000):
    """Input: names of the depression report (sDEP = self or mDEP = parental reports) and cardio-metabolic
       risk (CMR) marker; model structure and size of the graph.
       Creates a list of dictionaries to use as input for the elements arg of cytoscape graph.
       This only creates the core structure, style parameters are defined in style_sheet1.
    """
    # read data
    summ, _, esti, fail = read_res1(depname, cmrname)

    # Best fitting model
    if which_model == 'best':
        modstr = best_fit1(depname, cmrname).T
    elif which_model in fail:
        return 'fail'  # modstr = best_fit(depname, cmrname).T
    else:
        modstr = pd.DataFrame(model_structure[which_model]).T

    # Extract the estimated parameters from the result files (returns estimates and significance [= '*' or '']
    def extr_est(name=None, which=None, lamb=False, eta_corr=False, imp_corr=False, model_output=esti):
        df = model_output.loc[modstr.index[0]]
        if lamb:
            e, s = df.loc[(df.lhs == f'eta_{name}') & (df.op == '=~'), ['est', 'sign']].round(2).iloc[which]
        elif eta_corr:
            e, s = df.loc[(df.lhs == 'eta_dep') & (df.rhs == 'eta_cmr'), ['est', 'sign']].round(2).iloc[0]
        elif imp_corr:
            e, s = df.loc[df.label.str.contains('comv'), ['est', 'sign']].round(2).iloc[which]
        else:
            e, s = df.loc[df.label.str.contains(name)].iloc[::-1][['est', 'sign']].round(2).iloc[which]
        return e, s

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

    elm.append({'data': {'id': f'name_dep', 'label': 'Depression\nscore'}, 'classes': 'notes',
                'position': {'x': 120, 'y': pos_top + obs_pos}})
    elm.append({'data': {'id': f'name_cmr', 'label': '\n'.join(textwrap.wrap(get_label(cmrname), width=10))},
                'classes': 'notes',
                'position': {'x': 120, 'y': pos_bot - obs_pos}})

    for eta, pos in enumerate([pos_top, pos_bot]):
        # Eta factors nodes
        elm.append({'data': {'id': f'eta{eta}', 'label': 'Eta'}, 'classes': 'latent',
                    'position': {'x': width / 2 + rshift, 'y': pos}})

    #  Eta factors correlations
    elm.append({'data': {'source': 'eta0', 'target': 'eta1', 'firstname': 'eta_corr',
                         'label': '%.2f' % extr_est(eta_corr=True)[0]}})

    for i in range(1, nt + 1):

        elm.append({'data': {'source': f'imp_dep{i}', 'target': f'imp_cmr{i}', 'firstname': 'imp_corr',
                             'label': '%.2f' % extr_est(which=i - 1, imp_corr=True)[0]}})

        for eta, v in enumerate(vs):

            # ===== Other nodes =====
            # define vertical position
            p = [pos_top + obs_pos, pos_top + imp_pos] if v == 'dep' else [pos_bot - obs_pos, pos_bot - imp_pos]

            elm.extend([
                # Observed variables
                {'data': {'id': f'{v}{i}', 'firstname': f'{v.upper()}\n{timedic[v + str(i)]} yrs',
                          'label': summ.columns[(i - 1) + (nt * eta)]},
                 'classes': 'observed',
                 'position': {'x': ((width / nt) * i) - (width / nt) / 2 + rshift, 'y': p[0]}},
                # Impulses
                {'data': {'id': f'imp_{v}{i}', 'label': f'impulse {i}'},
                 'classes': 'latent',
                 'position': {'x': ((width / nt) * i) - (width / nt) / 2 + rshift, 'y': p[1]}}
            ])

            # ===== Edges: lambdas
            elm.append({'data': {'source': f'eta{eta}', 'target': f'{v}{i}', 'firstname': 'lambda',
                                 'label': '%.2f' % extr_est(f'{v}', i - 1, lamb=True)[0]}})
            # impulses link
            elm.append({'data': {'source': f'imp_{v}{i}', 'target': f'{v}{i}', 'firstname': 'imp_link'}})

            if i < nt:
                otherv = abs(eta - 1)
                # maAR and AR terms
                if modstr[f'ltAR_{v}'].iloc[0]:
                    elm.append({'data': {'source': f'{v}{i}', 'target': f'{v}{i + 1}',
                                         'weight': '%.2f' % (extr_est(f'^AR_{v}', i - 1)[0] * (
                                                 timedic[f'{v}{i + 1}'] - timedic[f'{v}{i}'])),
                                         'sign': extr_est(f'^AR_{v}', i - 1)[1],
                                         'label': f'AR{i}', 'firstname': 'direct'}})
                if modstr[f'maAR_{v}'].iloc[0]:
                    elm.append({'data': {'source': f'imp_{v}{i}', 'target': f'{v}{i + 1}',
                                         'weight': '%.2f' % (extr_est(f'^maAR_{v}', i - 1)[0] / (
                                                 timedic[f'{v}{i + 1}'] - timedic[f'{v}{i}'])),
                                         'sign': extr_est(f'^maAR_{v}', i - 1)[1],
                                         'label': f'maAR{i}', 'firstname': 'direct'}})
                # maCL and CL terms
                if modstr[f'ltCL_{v}'].iloc[0]:
                    elm.append({'data': {'source': f'{v}{i}', 'target': f'{vs[otherv]}{i + 1}',
                                         'weight': '%.2f' % (extr_est(f'^CL_{v}', i - 1)[0] * (
                                                 timedic[f'{vs[otherv]}{i + 1}'] - timedic[f'{v}{i}'])),
                                         'sign': extr_est(f'^CL_{v}', i - 1)[1],
                                         'label': f'CL{i}', 'firstname': 'direct'}})
                if modstr[f'maCL_{v}'].iloc[0]:
                    elm.append({'data': {'source': f'imp_{v}{i}', 'target': f'{vs[otherv]}{i + 1}',
                                         'weight': '%.2f' % (extr_est(f'^maCL_{v}', i - 1)[0] / (
                                                 timedic[f'{vs[otherv]}{i + 1}'] - timedic[f'{v}{i}'])),
                                         'sign': extr_est(f'^maCL_{v}', i - 1)[1],
                                         'label': f'maCL{i}', 'firstname': 'direct'}})
    return elm


# Define the stile of the graph
style_net1 = [
    # Noting down names
    {'selector': '.notes',
     'style': {'background-color': 'white', 'font-size': 18,
               'content': 'data(label)', 'color': 'k', 'font-weight': 'bold', 'text-halign': 'left',
               'text-valign': 'center', 'text-wrap': 'wrap'}},

    # Nodes - shape & color
    {'selector': '.observed',
     'style': {'shape': 'rectangle', 'height': 40, 'width': 80, 'border-width': 2, 'background-color': 'white',
               'border-color': 'k',
               'content': 'data(firstname)', 'color': 'grey', 'text-halign': 'center', 'text-valign': 'center',
               'text-wrap': 'wrap'}},

    {'selector': '.latent',
     'style': {'shape': 'round', 'height': 20, 'width': 20, 'border-width': 1, 'background-color': 'white',
               'border-color': 'silver'}},

    # Edges
    {'selector': 'edge[firstname *= "direct"]',  # directed paths
     'style': {'curve-style': 'straight', 'target-arrow-shape': 'vee', 'width': 3, 'arrow-scale': 1.2}},

    {'selector': 'edge[firstname *= "imp_link"]',  # impulses links
     'style': {'curve-style': 'straight', 'target-arrow-shape': 'vee', 'width': 1, 'arrow-scale': .8}},

    # Correlations
    {'selector': 'edge[firstname *= "eta_corr"]',
     'style': {'curve-style': 'unbundled-bezier', 'target-arrow-shape': 'vee', 'source-arrow-shape': 'vee', 'width': 1,
               'control-point-distances': [-450, -520, -530, -520, -450],
               'control-point-weights': [0.01, 0.20, 0.5, 0.80, 0.99],
               'label': 'data(label)', 'font-size': 15, 'text-background-color': 'silver',
               'text-background-opacity': .7}},

    {'selector': 'edge[firstname *= "imp_corr"]',
     'style': {'curve-style': 'unbundled-bezier', 'target-arrow-shape': 'vee', 'source-arrow-shape': 'vee', 'width': 1,
               'label': 'data(label)', 'font-size': 15, 'text-background-color': 'silver',
               'text-background-opacity': .7}},

    # Lambdas
    {'selector': 'edge[firstname *= "lambda"]',
     'style': {'curve-style': 'straight', 'target-arrow-shape': 'vee', 'width': 1, 'arrow-scale': .8,
               'label': 'data(label)', 'font-size': 15, 'text-background-color': 'silver',
               'text-background-opacity': .7}},

    # Dashed lines
    # {'selector':'edge[weight < 0.01][weight > -0.01]', 'style':{'line-style':'dashed'}},
    # {'selector':'edge[sign < 1]', 'style':{'line-style':'dashed'}},
]
# Set the color of each edge type and the distance between source and label displaying its weight (to avoid overlapping)
d = {'AR': ['red', 40], 'maAR': ['orange', 70], 'CL': ['green', 220], 'maCL': ['lightblue', 40]}

for c in d.keys():
    style_net1.extend([{'selector': f'[label *= "{c}"][sign > 0]',
                       'style': {'line-color': d[c][0], 'target-arrow-color': d[c][0],
                                 'source-label': 'data(weight)', 'source-text-offset': d[c][1], 'font-size': 20,
                                 'font-weight': 'bold',
                                 'text-background-color': d[c][0], 'text-background-opacity': .5}},
                      {'selector': f'[label *= "{c}"][sign < 1]',
                       'style': {'line-style': 'dashed', 'width': 1.5,
                                 'line-color': d[c][0], 'target-arrow-color': d[c][0],
                                 'source-label': 'data(weight)', 'source-text-offset': d[c][1], 'font-size': 18,
                                 'text-background-color': d[c][0], 'text-background-opacity': .5}}
                       ])


def make_table1(depname, cmrname, which_model='maCL_dep-maCL_cmr-maAR_dep-maAR_cmr'):
    """Input: names of the depression report (sDEP = self or mDEP = parental reports) and cardio-metabolic
       risk (CMR) marker; model structure.
       Extracts the fit measures for the specified model and stores into a table to be displayed next to the graph.
    """
    fitm = read_res1(depname, cmrname)[1]

    if which_model == 'best':
        which_model = fitm.index[fitm.aic == fitm.aic.min()][0]  # Best fitting model (lowest AIC)

    dt = pd.DataFrame(
        fitm.loc[which_model, ['npar', 'df', 'chisq', 'pvalue', 'cfi', 'tli', 'rmsea', 'srmr', 'aic', 'bic', ]])
    dt = dt.rename(columns={which_model: ' '}).round(3).astype('object')
    dt.insert(loc=0, column='Fit measures',
              value=['Number of parameters', 'Degrees of freedom', '\u03C7\u00b2', 'P-value (\u03C7\u00b2)',
                     'CFI', 'TLI', 'RMSEA', 'SRMR', 'AIC', 'BIC'])

    def format_row(x, form):
        try:
            return form.format(x)
        except ValueError:
            return x

    format_dict = {
        '{:d}': ['npar', 'df'],  # integers
        '{0:.1f}': ['aic', 'bic'],  # 1 decimal point
        '{0:.2f}': ['chisq', 'cfi', 'tli', 'rmsea', 'srmr'],  # 2 decimal points
        '{0:.3f}': ['pvalue'],  # 3 decimal points
    }
    for fmt in format_dict:
        dt.loc[format_dict[fmt]] = dt.loc[format_dict[fmt]].map(format_row, form=fmt)

    return dt