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
                          ' button to display the best fitting model (i.e. lowest AIC).'],
                         style=styles.TEXT),
                html.Hr(),

                # Input
                dbc.Row([
                    dbc.Col(width={'size': 4, 'offset': 1}, lg=3, md=2, sm=1,
                            children=[
                                html.H5(style=styles.SUB_TITLE1, children='Depression score'),
                                dcc.RadioItems(id=ids.DEP_SELECTION,
                                               options=dep_var_checklist(), value='sDEP',
                                               inputStyle={'margin-left': '20px', 'margin-right': '20px'},
                                               style=styles.TEXT),
                                html.Br(),
                                html.H5(style=styles.SUB_TITLE1, children='Cardio-metabolic marker'),
                                dcc.Dropdown(id=ids.CMR_SELECTION,
                                             options=cmr_var_checklist(), value='FMI',
                                             style=styles.TEXT)
                            ]),
                    dbc.Col(width={'size': 6}, lg=7, md=8, sm=9,
                            children=[
                                html.H5(style=styles.SUB_TITLE1, children='Model estimation'),
                                param_checklist('sDEP', 'FMI', p='lt'),
                                param_checklist('sDEP', 'FMI', p='ma'),
                                html.Div(style={'width': '40%', 'height': '35%', 'float': 'right'},
                                         children=[
                                             make_button('Best fit', ids.BESTFIT_BUTTON, 'silver'),
                                             make_button('Update model', ids.UPDATE_BUTTON, 'grey')
                                         ]),
                                html.Div(style={'width': '60%', 'height': '25%', 'float': 'left', 'color': 'red'},
                                         id=ids.FAILED_MODEL)
                            ])
                ], justify='between'), html.Hr(),

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
                    dbc.Col(width=10,
                            children=[
                                cyto.Cytoscape(id=ids.CYTO_GRAPH,
                                               layout={'name': 'preset', 'fit': False},
                                               style={'width': '100%', 'height': '100vh'},
                                               elements=make_net1('sDEP', 'FMI'),
                                               stylesheet=style_net1,
                                               minZoom=1, maxZoom=1)  # disable user zooming
                                ]),
                    # Table
                    dbc.Col(width=2,
                            children=[html.Br(),
                                      html.Div(id=ids.FITM_TABLE,
                                               children=[dbc.Table.from_dataframe(
                                                   df=make_table1('sDEP', 'FMI'),
                                                   color='light', striped=True, bordered=True, hover=True, size='lg')])
                                      ])
                ]),

                # Pop variable descriptives
                dbc.Offcanvas(style={'width': 700},
                              id=ids.POP1,
                              children=[dcc.Graph()],
                              title='',
                              is_open=False,
                              placement='end')
    ])
])
