from dash import dcc, html, register_page
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto

import definitions.elements_ids as ids
import definitions.layout_styles as styles

from definitions.general_funcs import bold_it, badge_it, underline_it, wrap_it

from definitions.gclpm_funcs import \
    dep_var_checklist, cmr_var_checklist, make_plot1
from definitions.riclpm_funcs import \
    read_res_riclpm, make_net_riclpm, style_net_riclpm, make_table_riclpm

register_page(__name__, path='/riclpm')

layout = dbc.Row([
    dbc.Col(width={'size': 10, 'offset': 1},
            children=[
                html.Br(),
                html.Div(['Results of the', bold_it('Random-Intercept Cross-Lag Panel Model'), '(ri-CLPM) described in the paper.',
                          wrap_it(2)],
                         style=styles.TEXT),
                html.Hr(),

                # Input
                dbc.Row([
                    dbc.Col(width={'size': 3},  # lg=3, md=2, sm=1,
                            children=[
                                html.Div(style={'margin-left': f'{int(styles.MARGIN_CHECKLIST[:-2])*2}px'},
                                         children=[
                                             html.H5(style=styles.SUB_TITLE1, children='Depression score'),
                                             dcc.RadioItems(id=ids.DEP_SELECTION_RICLPM,
                                                            options=dep_var_checklist(), value='sDEP',
                                                            inputStyle={'margin-left': styles.MARGIN_CHECKLIST,
                                                                        'margin-right': styles.MARGIN_CHECKLIST},
                                                            style=styles.TEXT),
                                             html.Br(),
                                             html.H5(style=styles.SUB_TITLE1, children='Cardio-metabolic marker'),
                                             dcc.Dropdown(id=ids.CMR_SELECTION_RICLPM,
                                                          options=cmr_var_checklist(), value='FMI',
                                                          clearable=False,
                                                          style=styles.TEXT)
                                         ])
                            ]),
                    dbc.Col(width={'size': 4},  # lg=3, md=2, sm=1,
                            children=[
                                html.Div(style={'margin-right': f'{int(styles.MARGIN_CHECKLIST[:-2])*2}px'},
                                         children=[
                                             html.H5(style=styles.SUB_TITLE1, children='Assume stationary over time:'),
                                             dbc.Row(children=[
                                                dcc.Checklist(id=ids.STAT_CHECKLIST_RICLPM,
                                                              options=[
                                                                {'label': 'Between-person (latent) effect',
                                                                 'value': 'l', 'disabled': True},
                                                                {'label': 'Within-person parameters',
                                                                 'value': 'p'}],
                                                              value=['l', 'p'],
                                                              style=styles.TEXT,
                                                              inputStyle={'margin-left': styles.MARGIN_CHECKLIST,
                                                                          'margin-right': styles.MARGIN_CHECKLIST},
                                                              labelStyle={'display': 'block'}),
                                                              ], align='start')
                                         ])
                            ]),
                    dbc.Col(width={'size': 5},  # lg=3, md=2, sm=1,
                            children=[
                                html.Div(['Using this selection pane, you can decide which depression report (i.e., self or \
                                          parental reports) and cardio-metabolic risk (CMR) marker you want to model. You can then inspect the \
                                          variables included in the model by clicking on the inspect icon or on the graph nodes directly. Check \
                                          the table on the right for info on the model fit.', wrap_it(),
                                          underline_it('Note'), ': by default, within-person parameter are assumed stationary over time \
                                          but this can be costumized using the tickboxes on the right'],
                                         style=styles.TEXT)
                            ])
                    ], class_name='g-0', justify='left'), html.Hr(),

                # Time plot
                html.Div([dbc.Accordion(start_collapsed=True,
                                        children=[
                                            dbc.AccordionItem(title='Inspect the variables included in this model',
                                                              children=[
                                                                  dcc.Graph(id=ids.TIME_GRAPH_RICLPM,
                                                                            figure=make_plot1('sDEP', 'FMI'))],
                                                              style=styles.TEXT)
                                        ])
                          ]),

                # Results
                dbc.Row([
                    # Network
                    dbc.Col(width=9,
                            children=[
                                cyto.Cytoscape(id=ids.CYTO_GRAPH_RICLPM,
                                               layout={'name': 'preset', 'fit': False},
                                               style={'width': '100%', 'height': '100vh'},
                                               elements=make_net_riclpm('sDEP', 'FMI'),
                                               stylesheet=style_net_riclpm,
                                               minZoom=1, maxZoom=1)  # disable user zooming
                                ]),
                    # Table
                    dbc.Col(width=3,
                            children=[html.Br(),
                                      html.Div(id=ids.FITM_TABLE_RICLPM,
                                               children=[dbc.Table.from_dataframe(
                                                   df=make_table_riclpm('sDEP', 'FMI'),
                                                   style={'font-size': styles.CLPM_TABLE_TEXT},
                                                   color='light', striped=True, bordered=True, hover=True, size='lg')])
                                      ])
                ], justify='between'),

                # Pop variable descriptives
                dbc.Offcanvas(style={'width': styles.OFFCANVAS_WIDTH},
                              id=ids.POP1_RICLPM,
                              children=[dcc.Graph()],
                              title='',
                              is_open=False,
                              placement='end')
    ])
])
