from dash import dcc, html, register_page
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto

import definitions.elements_ids as ids
import definitions.layout_styles as styles

from definitions.general_funcs import bold_it, badge_it, underline_it, wrap_it

from definitions.clpm_funcs import \
    dep_var_checklist, cmr_var_checklist, param_checklist, make_button, make_plot1, make_net1, style_net1, make_table1


register_page(__name__, path='/clpm')

layout = dbc.Row([
    dbc.Col(width={'size': 10, 'offset': 1},
            children=[
                html.Br(),
                html.Div(['Results of the generalized', bold_it('cross-lag panel model'), 'described as model 1 in the paper.',
                          wrap_it(2), 'Using the selection pane below, you can decide which depression report (i.e., self or \
                          parental reports) and cardio-metabolic risk (CMR) marker you want to model. You can then inspect the \
                          variables included in the model by clicking on the inspect icon or on the graph nodes directly. Check \
                          the table on the right for info on the model fit.', wrap_it(),
                          underline_it('Note'), ': by default, the "classic" cross-lag panel model is presented, but the parameter \
                          conbination can be constumized using the tickboxes on the right (don\'t forget to hit the ',
                          badge_it('Update model', 'grey'), ' button to see the changes). Hit the ', badge_it('Best fit', 'silver'),
                          ' button to display the best fitting model (i.e. lowest BIC).'],
                         style=styles.TEXT),
                html.Hr(),

                # Input
                dbc.Row([
                    dbc.Col(width={'size': 'auto'}, # lg=3, md=2, sm=1,
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
                                                          style=styles.TEXT)
                                         ])
                            ]),
                    dbc.Col(width={'size': 6}, # lg=7, md=8, sm=9,
                            children=[
                                html.H5(style=styles.SUB_TITLE1, children='Model structure'),
                                dbc.Col(width={'size': 'auto'},
                                        children=[param_checklist('sDEP', 'FMI', p='lt')]),
                                dbc.Col(width={'size': 'auto'},
                                        children=[param_checklist('sDEP', 'FMI', p='ma')]),

                                html.Div(id=ids.FAILED_MODEL)
                            ]),
                    dbc.Col(width={'size': 'auto'}, # lg=3, md=2, sm=1,
                            children=[
                                html.Div(style={'margin-right': f'{int(styles.MARGIN_CHECKLIST[:-2])*2}px'},
                                         children=[
                                             html.H5(style=styles.SUB_TITLE1, children='Stationarity assumtion'),
                                             dbc.Row(children=[
                                                dcc.Checklist(id='stat-checklist',
                                                              options=[
                                                                {'label': 'Between-person (latent) effect',
                                                                 'value': 'l'},
                                                                {'label': 'Within-person parameters',
                                                                 'value': 'p'}],
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
                    ], class_name='g-0', justify='between'), html.Hr(),

                # Time plot
                html.Div([dbc.Accordion(start_collapsed=True,
                                        children=[
                                            dbc.AccordionItem(title='Inspect the variables included in this model',
                                                              children=[
                                                                  dcc.Graph(id=ids.TIME_GRAPH,
                                                                            figure=make_plot1('sDEP', 'FMI'))],
                                                              style=styles.TEXT)
                                        ])
                          ]),

                # Results
                dbc.Row([
                    # Network
                    dbc.Col(width=9,
                            children=[
                                cyto.Cytoscape(id=ids.CYTO_GRAPH,
                                               layout={'name': 'preset', 'fit': False},
                                               style={'width': '100%', 'height': '100vh'},
                                               elements=make_net1('sDEP', 'FMI'),
                                               stylesheet=style_net1,
                                               minZoom=1, maxZoom=1)  # disable user zooming
                                ]),
                    # Table
                    dbc.Col(width=3,
                            children=[html.Br(),
                                      html.Div(id=ids.FITM_TABLE,
                                               children=[dbc.Table.from_dataframe(
                                                   df=make_table1('sDEP', 'FMI'),
                                                   style={'font-size': '12px'},
                                                   color='light', striped=True, bordered=True, hover=True, size='lg')])
                                      ])
                ], justify='between'),

                # Pop variable descriptives
                dbc.Offcanvas(style={'width': 700},
                              id=ids.POP1,
                              children=[dcc.Graph()],
                              title='',
                              is_open=False,
                              placement='end')
    ])
])
