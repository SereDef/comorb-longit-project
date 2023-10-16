from dash import callback, Output, Input

import definitions.elements_ids as ids

from definitions.clnm_funcs import make_table2


# Adjust order in CI tables based on selected column
@callback(
    Output(ids.TEMP_CI_TAB, 'data'),
    Input(ids.TEMP_CI_TAB, 'sort_by'), prevent_initial_call=True
)
def update_temp_tab(sort_by):

    dff = make_table2('t').sort_values(
        sort_by[0]['column_id'], ascending=sort_by[0]['direction'] == 'desc', inplace=False)

    return dff.to_dict('records')


@callback(
    Output(ids.CONT_CI_TAB, 'data'),
    Input(ids.CONT_CI_TAB, 'sort_by'), prevent_initial_call=True
)
def update_cont_tab(sort_by):

    dff = make_table2('c').sort_values(
        sort_by[0]['column_id'], ascending=sort_by[0]['direction'] == 'desc', inplace=False)

    return dff.to_dict('records')


@callback(
    Output(ids.BETW_CI_TAB, 'data'),
    Input(ids.BETW_CI_TAB, 'sort_by'), prevent_initial_call=True
)
def update_betw_tab(sort_by):

    dff = make_table2('b').sort_values(
        sort_by[0]['column_id'], ascending=sort_by[0]['direction'] == 'desc', inplace=False)

    return dff.to_dict('records')
