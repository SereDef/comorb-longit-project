from dash import dcc, html, register_page
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
import dash_cytoscape as cyto

import definitions.elements_ids as ids
import definitions.layout_styles as styles

from definitions.general_funcs import bold_it, badge_it, underline_it, wrap_it, info_icon, \
    dep_var_checklist, cmr_var_checklist, stat_checklist, temp_plot

from definitions.riclpm_funcs import \
    read_res_riclpm, make_net_riclpm, style_net_riclpm, make_table_riclpm

register_page(__name__, path='/riclpm')

layout = dbc.Row([
    dbc.Col(width={'size': 10, 'offset': 1},
            children=[
                html.Br(),
                html.Div(['Results of the', bold_it('Random-Intercept Cross-Lag Panel Model'), '(RI-CLPM) described in the paper [link].',
                          wrap_it(2)],
                         style=styles.TEXT),

                # Selection pane
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
                                dbc.Col(width={'size': 3, 'offset': 1},  # lg=3, md=2, sm=1,
                                        children=[
                                            html.Div(style={'margin-right': f'{int(styles.MARGIN_CHECKLIST[:-2])*2}px'},
                                                     children=[
                                                         html.H5(style=styles.SUB_TITLE1, children='Assume stationary over time:'),
                                                         dbc.Row(children=[
                                                            dcc.Checklist(id=ids.STAT_CHECKLIST_RICLPM,
                                                                          options= stat_checklist(disable_etas=True),
                                                                          value=['l', 'p'],
                                                                          style=styles.TEXT,
                                                                          inputStyle={'margin-left': styles.MARGIN_CHECKLIST,
                                                                                      'margin-right': styles.MARGIN_CHECKLIST},
                                                                          labelStyle={'display': 'block'}),
                                                                          ], align='start'),
                                                         html.Br(),
                                                         html.H5(style=styles.SUB_TITLE1, children='Stratify by sex:'),
                                                         dcc.Checklist(id=ids.SEX_SELECTION_RICLPM,
                                                                       options=[
                                                                            {'label': 'Female', 'value': 'f'},
                                                                            {'label': 'Male', 'value': 'm'}
                                                                        ], value=[],
                                                                        inputStyle={'margin-left': styles.MARGIN_CHECKLIST,
                                                                                    'margin-right': styles.MARGIN_CHECKLIST},
                                                                        style=styles.TEXT),
                                                     ])
                                        ]),
                                dbc.Col(width={'size': 5},  # lg=3, md=2, sm=1,
                                        children=[
                                            html.Div([bold_it('Instructions'), wrap_it(),
                                                      'Using this selection pane, you can decide which depression report \
                                                      (i.e., self or parental reports) and cardio-metabolic risk (CMR) \
                                                      marker to model.', wrap_it(), 'You can then inspect the variables \
                                                      included in the model by clicking on the ', info_icon(),
                                                     ' tab or on the graph nodes directly.', wrap_it(),
                                                      'Check the table on the right for info on the model fit.', wrap_it(),
                                                      underline_it('Note'), ': by default, within-person parameter are \
                                                      assumed stationary over time and the full sample is considered \
                                                      but this can be customized using the checkboxes on the right.'],
                                                     style=styles.TEXT)
                                        ])
                                ], class_name='g-0', justify='left'
                            ),
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
                                dmc.AccordionPanel(dcc.Graph(id=ids.TIME_GRAPH_RICLPM,
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
                                cyto.Cytoscape(id=ids.CYTO_GRAPH_RICLPM,
                                               layout={'name': 'preset', 'fit': False},
                                               style={'width': '100%', 'height': styles.CLPM_HEIGHT},
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
                                                   color='light', striped=True, bordered=True, hover=True, size='md')])
                                      ])
                ], justify='between'),

                html.Div(style=styles.TEXT,
                         children=['* ', underline_it('Note'), ': the estimates above are standardized and scaled with \
                         respect to time, so they represent a', bold_it('yearly'), 'change.']),

                # Pop variable descriptives
                dbc.Offcanvas(style={'width': styles.OFFCANVAS_WIDTH},
                              id=ids.POP1_RICLPM,
                              children=[dcc.Graph()],
                              title='',
                              is_open=False,
                              placement='end')
    ])
])
