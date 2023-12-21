from dash import dcc, html, register_page
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
import dash_cytoscape as cyto

import definitions.elements_ids as ids
import definitions.layout_styles as styles

from definitions.general_funcs import bold_it, badge_it, underline_it, wrap_it, info_icon, \
    dep_var_checklist, cmr_var_checklist, stat_checklist, temp_plot

from definitions.gclpm_funcs import param_checklist, make_button, \
    make_net_gclpm, style_net_gclpm, make_table_gclpm


register_page(__name__, path='/gclpm')

layout = dbc.Row([
    dbc.Col(width={'size': 10, 'offset': 1},
            children=[
                html.Br(),
                html.Div(['Results of the', bold_it('Generalized Cross-Lag Panel Model'), '(gCLPM) described in the paper [link].',
                          wrap_it(2), 'Using the selection pane below, you can decide which depression report (i.e., self or \
                          parental reports) and cardio-metabolic risk (CMR) marker to model.', wrap_it(),
                          'You can then inspect the variables included in the model by clicking on the ', info_icon(),
                          ' tab or on the graph nodes directly.', wrap_it(),
                          'Check the table on the right for info on the model fit.', wrap_it(),
                          underline_it('Note'), ': by default, the "classic" cross-lag panel model is presented, but the \
                          parameter combination and stationarity assumptions can be customized using the checkboxes on \
                          the right (don\'t forget to hit the ', badge_it('Update model', 'grey'),
                          ' button to see the changes).', wrap_it(), 'Hit the ', badge_it('Best fit', 'silver'),
                          ' button to display the best fitting model (i.e. lowest BIC).', wrap_it(2)],
                         style=styles.TEXT),

                html.Div(style={'background-color': styles.SELECTION_PANEL_COLOR,
                                'border-radius': '30px'},
                         children=[
                            html.Br(),
                            # Input
                            dbc.Row([
                                dbc.Col(width={'size': 3},  # lg=3, md=2, sm=1,
                                        children=[
                                            html.Div(style={'margin-left': f'{int(styles.MARGIN_CHECKLIST[:-2])*2}px'},
                                                     children=[
                                                         html.H5(style=styles.SUB_TITLE1, children='Depression score'),
                                                         dcc.RadioItems(id=ids.DEP_SELECTION,
                                                                        options=dep_var_checklist(), value='sDEP',
                                                                        inputStyle={'margin-left': styles.MARGIN_CHECKLIST,
                                                                                    'margin-right': styles.MARGIN_CHECKLIST},
                                                                        style=styles.TEXT),
                                                         html.Br(),
                                                         html.H5(style=styles.SUB_TITLE1, children='Cardio-metabolic marker'),
                                                         dcc.Dropdown(id=ids.CMR_SELECTION,
                                                                      options=cmr_var_checklist(), value='FMI',
                                                                      clearable=False,
                                                                      style=styles.TEXT)
                                                     ])
                                        ]),
                                dbc.Col(width={'size': 6},  # lg=7, md=8, sm=9,
                                        children=[
                                            html.H5(style=styles.SUB_TITLE1, children='Model structure'),
                                            dbc.Col(width={'size': 'auto'},
                                                    children=[param_checklist('sDEP', 'FMI', p='lt')]),
                                            dbc.Col(width={'size': 'auto'},
                                                    children=[param_checklist('sDEP', 'FMI', p='ma')]),

                                            html.Div(id=ids.FAILED_MODEL)
                                        ]),
                                dbc.Col(width={'size': 'auto'},  # lg=3, md=2, sm=1,
                                        children=[
                                            html.Div(style={'margin-right': f'{int(styles.MARGIN_CHECKLIST[:-2])*2}px'},
                                                     children=[
                                                         html.H5(style=styles.SUB_TITLE1,
                                                                 children='Assume stationary over time:'),
                                                         dbc.Row(children=[
                                                            dcc.Checklist(id=ids.STAT_CHECKLIST,
                                                                          options=stat_checklist(),
                                                                          value=['p'],
                                                                          style=styles.TEXT,
                                                                          inputStyle={'margin-left': styles.MARGIN_CHECKLIST,
                                                                                      'margin-right': styles.MARGIN_CHECKLIST},
                                                                          labelStyle={'display': 'block'}),
                                                                          ], align='start'),
                                                         dbc.Row([html.Br()]),
                                                         dbc.Row(children=[
                                                                     dbc.Col(width={'size': 'auto'},
                                                                             children=[make_button('Best fit', ids.BESTFIT_BUTTON, 'silver')]),
                                                                     dbc.Col(width={'size': 'auto'},
                                                                             children=[make_button('Update model', ids.UPDATE_BUTTON, 'grey')])
                                                                 ], class_name='g-1', justify='end', align='end')
                                                     ])
                                        ])
                                ], class_name='g-0', justify='between'),
                             html.Br()
                         ]),

                # Time plot
                html.Div([
                    html.Br(),
                    dmc.Accordion(  # disableChevronRotation=True,
                        children=[
                            dmc.AccordionItem([
                                dmc.AccordionControl('Inspect the variables included in this model',
                                                     icon=info_icon(),
                                                     style=styles.TEXT),
                                dmc.AccordionPanel(dcc.Graph(id=ids.TIME_GRAPH,
                                                             figure=temp_plot('sDEP', 'FMI'))),
                            ], value='info')
                        ], radius='xl', variant='separated'
                    )
                ]),

                # Results
                dbc.Row([
                    # Network
                    dbc.Col(width=9,
                            children=[
                                cyto.Cytoscape(id=ids.CYTO_GRAPH,
                                               layout={'name': 'preset', 'fit': False},
                                               style={'width': '100%', 'height': '55vh'},
                                               elements=make_net_gclpm('sDEP', 'FMI'),
                                               stylesheet=style_net_gclpm,
                                               minZoom=1, maxZoom=1)  # disable user zooming
                                ]),
                    # Table
                    dbc.Col(width=3,
                            children=[html.Br(),
                                      html.Div(id=ids.FITM_TABLE,
                                               children=[dbc.Table.from_dataframe(
                                                   df=make_table_gclpm('sDEP', 'FMI'),
                                                   style={'font-size': styles.CLPM_TABLE_TEXT},
                                                   color='light', striped=True, bordered=True, hover=True, size='lg')])
                                      ])
                ], justify='between'),

                html.Div(style=styles.TEXT,
                         children=['* ', underline_it('Note'), ': the (unstandardized) estimates above are scaled with \
                         respect to time, so they represent a', bold_it('yearly'), 'change.']),

                # Pop variable descriptives
                dbc.Offcanvas(style={'width': styles.OFFCANVAS_WIDTH},
                              id=ids.POP1,
                              children=[dcc.Graph()],
                              title='',
                              is_open=False,
                              placement='end')
    ])
])
