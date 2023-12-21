from dash import dcc, callback, Output, Input, State, ctx, no_update
import dash_bootstrap_components as dbc

import pandas as pd

import definitions.elements_ids as ids
from definitions.general_funcs import desc_plot, dep_var_checklist, cmr_var_checklist, temp_plot

from definitions.gclpm_funcs import param_checklist, model_structure, \
    best_fit_gclpm, make_net_gclpm, style_net_gclpm, make_table_gclpm


# Control variable selection (if mother reports are selected, only some CMR markers are available)
@callback(
    Output(ids.DEP_SELECTION, 'options'),
    Output(ids.CMR_SELECTION, 'options'),
    Input(ids.DEP_SELECTION, 'value'),
    Input(ids.CMR_SELECTION, 'value')
)
def update_dropdown_options(dep_selection, cmr_selection):
    if cmr_selection not in ['FMI', 'LMI', 'BMI', 'waist_circ', 'total_fatmass', 'total_leanmass']:
        return dep_var_checklist(disable_maternal=True), no_update

    if dep_selection == 'sDEP':
        return dep_var_checklist(), cmr_var_checklist('s')

    return dep_var_checklist(), cmr_var_checklist('m')


# Control variable selection, updating the scatterplot and the ticks on the model selection pane
@callback(
    Output(ids.TIME_GRAPH, 'figure'),
    Input(ids.DEP_SELECTION, 'value'),
    Input(ids.CMR_SELECTION, 'value')
)
def update_time_plot(dep_selection, cmr_selection):
    # Prevent update if the selection is not allowed (i.e., combination of markers not available)
    if dep_selection == 'mDEP' and cmr_selection not in [
                                   'FMI', 'LMI', 'BMI', 'waist_circ', 'total_fatmass', 'total_leanmass']:
        return no_update

    # Default output: only long term parameters
    return temp_plot(dep_selection, cmr_selection)


# Based on variable and model selection display graph (or notify if model didn't converge)
@callback(
    Output(ids.CYTO_GRAPH, 'elements'),
    Output(ids.FITM_TABLE, 'children'),
    Output(ids.FAILED_MODEL, 'children'),
    Output(ids.LT_CHECKLIST, 'value'),
    Output(ids.MA_CHECKLIST, 'value'),

    Input(ids.DEP_SELECTION, 'value'),
    Input(ids.CMR_SELECTION, 'value'),
    Input(ids.STAT_CHECKLIST, 'value'),
    Input(ids.UPDATE_BUTTON, 'n_clicks'),
    Input(ids.BESTFIT_BUTTON, 'n_clicks'),
    State(ids.LT_CHECKLIST, 'value'),
    State(ids.MA_CHECKLIST, 'value')
)
def update_graph(dep_selection, cmr_selection, stationary, click_manual, click_best, lt_checklist, ma_checklist):
    # Prevent update if the selection is not allowed (i.e., combination of markers not available)
    lambdas = 'stat' if 'l' in stationary else 'free'
    params = 'stat' if 'p' in stationary else 'free'

    if dep_selection == 'mDEP' and cmr_selection not in [
                                   'FMI', 'LMI', 'BMI', 'waist_circ', 'total_fatmass', 'total_leanmass']:
        return no_update, no_update, no_update, no_update, no_update

    if ctx.triggered_id == ids.UPDATE_BUTTON:

        checked = lt_checklist + ma_checklist
        series = pd.Series([1 if x in checked else 0 for x in model_structure.index], index=model_structure.index)
        retrieve_model = model_structure.columns[model_structure.eq(series, axis=0).all()]

        if len(retrieve_model) == 0:
            return no_update, no_update, 'Sorry, I did not estimate this model.', no_update, no_update

        model_name = retrieve_model[0]

        updated_graph = make_net_gclpm(dep_selection, cmr_selection, lambdas, params, which_model=model_name)

        if updated_graph == 'fail':  # prevents any single output updating
            return no_update, no_update, 'Sorry, this model did not converge.', no_update, no_update

        tab_df = make_table_gclpm(dep_selection, cmr_selection, lambdas, params, which_model=model_name)
        update_tab = dbc.Table.from_dataframe(df=tab_df, color='light', striped=True, bordered=True, hover=True,
                                              size='lg')

        return updated_graph, update_tab, None, no_update, no_update

    elif ctx.triggered_id == ids.BESTFIT_BUTTON:

        updated_graph = make_net_gclpm(dep_selection, cmr_selection, lambdas, params, which_model='best')

        tab_df = make_table_gclpm(dep_selection, cmr_selection, lambdas, params, which_model='best')
        update_tab = dbc.Table.from_dataframe(df=tab_df, color='light', striped=True, bordered=True, hover=True,
                                              size='lg')

        return updated_graph, update_tab, None, \
            best_fit_gclpm(dep_selection, cmr_selection, lambdas, params, 'lt'), \
               best_fit_gclpm(dep_selection, cmr_selection, lambdas, params, 'ma')

    updated_graph = make_net_gclpm(dep_selection, cmr_selection, lambdas, params)

    if updated_graph == 'fail':  # prevents output updating
        return no_update, no_update, 'Sorry, this model did not converge.', no_update, no_update

    tab_df = make_table_gclpm(dep_selection, cmr_selection, lambdas, params)
    update_tab = dbc.Table.from_dataframe(df=tab_df, color='light', striped=True, bordered=True, hover=True, size='lg')

    return make_net_gclpm(dep_selection, cmr_selection, lambdas, params), update_tab, None, \
        ['ltCL_dep', 'ltCL_cmr', 'ltAR_dep', 'ltAR_cmr'], []


# Display info upon tapping on node
@callback(
    Output(ids.POP1, 'is_open'),
    Output(ids.POP1, 'title'),
    Output(ids.POP1, 'children'),
    Input(ids.CYTO_GRAPH, 'tapNodeData'),
    [State(ids.POP1, 'is_open')]
)
def display_tap_node_data1(node, is_open):
    if node:
        titl = node['label']

        if titl == 'Eta':
            return is_open, no_update, no_update

        plot = dcc.Graph(figure=desc_plot(titl))

        return not is_open, titl, plot

    return is_open, no_update, no_update
