from dash import html, dcc, dash_table, register_page
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto

import definitions.elements_ids as ids
from definitions.general_funcs import bold_it
from definitions.csnm_funcs import time_marks3, make_net3, style_net3

register_page(__name__, path='/csnm')

layout = dbc.Row([
    dbc.Col(width={'size': 10, 'offset': 1},
            children=[
                html.Br(),
                html.Span(['Results of each', bold_it('cross-sectional network model'), 'conducted as follow-up analyses. \
                Please select the timepoint of interest from the slides to visualize the network structure and corresponding \
                centrality indices.']),
                html.Br(), html.Hr(),

                # Slider
                html.Span(children=[html.Div('Select a time point:'), html.Br(style={'line-height': '5'})]),
                dcc.Slider(id=ids.CROS_NET_SLIDER, min=9.7, max=24.2, step=None, value=9.7, marks=time_marks3, included=False),

                # Results
                dbc.Row([
                    # Network
                    dbc.Col(width=6,
                            children=[cyto.Cytoscape(id=ids.CROS_NET,
                                                     layout={'name': 'preset'},
                                                     style={'width': '40%', 'height': '80%',
                                                            'position': 'absolute', 'left': '10vw', 'top': '35vh'},
                                                     minZoom=1, maxZoom=1,  # disable user zooming
                                                     elements=make_net3(9.7)[0],
                                                     stylesheet=style_net3)]),
                    # Table
                    dbc.Col(width={'size': 5, 'offset': 1},
                            children=[html.Br(), html.Br(),
                                      dash_table.DataTable(
                                          id=ids.CROS_CI_TAB,
                                          columns=[{'name': i, 'id': i} for i in make_net3(9.7)[1].columns],
                                          sort_action='custom', sort_mode='single', sort_by=[],
                                          fixed_columns={'headers': True, 'data': 1},  # Fix node name column
                                          style_header={'fontWeight': 'bold'},
                                          style_cell={'fontSize': 20, 'font-family': 'sans-serif'},
                                          style_cell_conditional=[{'if': {'column_id': 'Node'}, 'width': '250px'}],
                                          style_data={'whiteSpace': 'normal', 'height': 'auto', 'lineHeight': '20px',
                                                      'minWidth': '100px', 'width': '100px', 'maxWidth': '100px'},
                                          style_table={'overflowX': 'auto', 'minWidth': '100%'})
                                      ]),

                    # Pop variable descriptives
                    dbc.Offcanvas(id=ids.POP3,
                                  children=[dcc.Graph()],
                                  title='Lab name',
                                  is_open=False,
                                  placement='end',
                                  style={'width': 700}),  # backdrop
                ])
            ])
])