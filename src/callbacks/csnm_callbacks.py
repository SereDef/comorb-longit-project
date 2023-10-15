from dash import dcc, callback, Output, Input, State, ctx, no_update

import definitions.elements_ids as ids

from definitions.general_funcs import desc_plot
from definitions.csnm_funcs import make_net3


# Display cross-sectional network based on selected time point
@callback(
    Output(ids.CROS_NET, 'elements'),
    Output(ids.CROS_CI_TAB, 'data'),
    Input(ids.CROS_NET_SLIDER, 'value'),
    Input(ids.CROS_CI_TAB, 'sort_by')
)
def update_crosnet(timepoint, sort_by):

    if ctx.triggered_id == ids.CROS_CI_TAB:
        dff = make_net3(timepoint)[1].sort_values(sort_by[0]['column_id'],
                                                  ascending=sort_by[0]['direction'] == 'desc', inplace=False)

        return no_update, dff.to_dict('records')

    return make_net3(timepoint)[0], make_net3(timepoint)[1].to_dict('records')


# Display info upon tapping on node
@callback(
    Output(ids.POP3, 'is_open'),
    Output(ids.POP3, 'title'),
    Output(ids.POP3, 'children'),
    Input(ids.CROS_NET, 'tapNodeData'),
    [State(ids.POP3, 'is_open')],
)
def display_tap_node_data3(node, is_open):
    if node:
        titl = node['label']
        plot = dcc.Graph(figure=desc_plot(node['id']))
        return not is_open, titl, plot

    return is_open, no_update, no_update
