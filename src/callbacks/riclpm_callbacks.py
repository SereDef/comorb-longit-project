from dash import dcc, callback, Output, Input, State, ctx, no_update
import dash_bootstrap_components as dbc

import definitions.elements_ids as ids
from definitions.general_funcs import desc_plot, dep_var_checklist, cmr_var_checklist, temp_plot

from definitions.riclpm_funcs import \
    read_res_riclpm, make_net_riclpm, style_net_riclpm, make_table_riclpm


# Control variable selection (if mother reports are selected, only some CMR markers are available)
@callback(
    Output(ids.DEP_SELECTION_RICLPM, 'options'),
    Output(ids.CMR_SELECTION_RICLPM, 'options'),
    Input(ids.DEP_SELECTION_RICLPM, 'value'),
    Input(ids.CMR_SELECTION_RICLPM, 'value')
)
def update_dropdown_options(dep_selection, cmr_selection):

    if cmr_selection not in ['FMI', 'LMI', 'BMI', 'waist_circ', 'total_fatmass', 'total_leanmass']:
        return dep_var_checklist(disable_maternal=True), no_update

    if dep_selection == 'sDEP':
        return dep_var_checklist(), cmr_var_checklist('s')

    return dep_var_checklist(), cmr_var_checklist('m')


# Control variable selection, updating the scatterplot and the ticks on the model selection pane
@callback(
    Output(ids.TIME_GRAPH_RICLPM, 'figure'),
    Input(ids.DEP_SELECTION_RICLPM, 'value'),
    Input(ids.CMR_SELECTION_RICLPM, 'value')
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
    Output(ids.CYTO_GRAPH_RICLPM, 'elements'),
    Output(ids.FITM_TABLE_RICLPM, 'children'),

    Input(ids.DEP_SELECTION_RICLPM, 'value'),
    Input(ids.CMR_SELECTION_RICLPM, 'value'),
    Input(ids.STAT_CHECKLIST_RICLPM, 'value'),
    Input(ids.SEX_SELECTION_RICLPM, 'value'),
)
def update_graph(dep_selection, cmr_selection, stationary, sex_selection):
    # Prevent update if the selection is not allowed (i.e., combination of markers not available)
    params = 'stat' if 'p' in stationary else 'free'

    sex = '' if len(sex_selection) != 1 else sex_selection[0]

    if dep_selection == 'mDEP' and cmr_selection not in [
                                   'FMI', 'LMI', 'BMI', 'waist_circ', 'total_fatmass', 'total_leanmass']:
        return no_update, no_update

    updated_graph = make_net_riclpm(dep_selection, cmr_selection, params, sex)

    tab_df = make_table_riclpm(dep_selection, cmr_selection, params, sex)
    update_tab = dbc.Table.from_dataframe(df=tab_df, color='light', striped=True, bordered=True, hover=True, size='lg')

    return updated_graph, update_tab


# Display info upon tapping on node
@callback(
    Output(ids.POP1_RICLPM, 'is_open'),
    Output(ids.POP1_RICLPM, 'title'),
    Output(ids.POP1_RICLPM, 'children'),
    Input(ids.CYTO_GRAPH_RICLPM, 'tapNodeData'),
    [State(ids.POP1_RICLPM, 'is_open')]
)
def display_tap_node_data1(node, is_open):
    if node:
        titl = node['label']

        if titl == 'Eta':
            return is_open, no_update, no_update

        plot = dcc.Graph(figure=desc_plot(titl))

        return not is_open, titl, plot

    return is_open, no_update, no_update
