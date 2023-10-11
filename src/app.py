#!/usr/bin/env python
# coding: utf-8

from dash import Dash, html, dcc, callback, Output, Input, State, ctx, no_update, dash_table, get_asset_url
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto

import pandas as pd

from backend_funcs import get_label, desc_plot, model_structure, all_desc, age_desc, plot_overview, read_res1, best_fit1, make_plot1, make_net1, stylenet1, make_table1, make_net2, make_netstyle2, make_tab2, read_res3, timemarks3, make_net3, stylenet3


def badge_it(text, color):
    return dbc.Badge(text, color=color, style={'padding':'4px 5px'})

def bold_it(text):
    return html.Span(f' {text} ', style={'font-weight':'bold'})

def undeline_it(text):
    return html.Span(f'{text}', style={'text-decoration':'underline'})

def wrap_it(n_spaces=1):
    return html.Span( [html.Br()] * n_spaces )

def DEPmarker_checklist(disable_maternal=False):
    opts=[{'label': 'Self-reported',  'value':'sDEP', 'disabled':False},
          {'label': 'Maternal report','value':'mDEP', 'disabled':disable_maternal}]
    return opts
    
def CMRmarker_checklist(reporter='s'):
    # Read all fitted model names 
    #allfiles = [x.split('.')[0] for x in os.listdir(PATH+'mod1')]
    #moth = [s[5:]+'_' for s in allfiles if "mDEP" in s]
    
    sortd = ['FMI','LMI','BMI','waist_circ','android_fatmass','total_fatmass','total_leanmass', 
             'tot_chol','HDL_chol','LDL_chol','insulin','triglyc','CRP']
    
    opts = [{'value':v, 'label':get_label(v) + f' ({v})' if len(v)<4 else get_label(v), 'disabled':False } for v in sortd]

    if reporter=='m':
        sortd_m = ['FMI','LMI','BMI','waist_circ','total_fatmass','total_leanmass']
        for i in opts: 
            if i['value'] not in sortd_m: i['disabled'] = True
    
    return opts


def param_checklist(depname, cmrname, p='lt', best=False):
    pref = '' if p == 'lt' else 'ma'
    cols = ['crimson','green'] if p=='lt' else ['orange', 'lightblue']
    position = 'left' if p=='lt' else 'right'

    if best: 
        val = best_fit1(depname, cmrname, list1=p)
    else: 
        val = ['ltCL_dep','ltCL_cmr','ltAR_dep','ltAR_cmr'] if p == 'lt' else []
    
    return  html.Div(style={'width':'50%','height':'65%','float':position},
                  children=[ dcc.Checklist(id =f'{p}-checklist',
                   options=[{'label': html.Span([badge_it(f'{pref}AR', cols[0]),' depression']),   'value': f'{p}AR_dep'},
                            {'label': html.Span([badge_it(f'{pref}AR', cols[0]),' cardio-metabolic risk']),'value': f'{p}AR_cmr'},
                            {'label': html.Span([badge_it(f'{pref}CL', cols[1]),' depression \u290F cardio-metab.']),'value': f'{p}CL_dep'},
                            {'label': html.Span([badge_it(f'{pref}CL', cols[1]),' cardio-metab. \u290F depression']),'value': f'{p}CL_cmr'}],
                    value = val,                        
                inputStyle={ 'margin-left':'20px','margin-right':'20px'}, labelStyle = {'display': 'block'}) ])



app = Dash(__name__, external_stylesheets=[dbc.themes.LITERA, dbc.icons.BOOTSTRAP], suppress_callback_exceptions=True)
server = app.server

app.layout = dbc.Container([
    # Title
    html.H1('Longitudinal modelling of the co-development of depression and cardio-metabolic risk from childhood to young adulthood',
             style={'textAlign':'center', 'font-weight':'bold'}),
    html.Br(), # space
    # Main body
    dbc.Row([ dbc.Col([ 
        dcc.Tabs(id="tabs", 
                 children=[dcc.Tab(label='Data overview', value='tab-0'),
                           dcc.Tab(label='Cross-lag panel model', value='tab-1'), # style={''}
                           dcc.Tab(label='Cross-lag network analysis', value='tab-2'),
                           dcc.Tab(label='Cross-sectional network analysis', value='tab-3') ], 
                 value='tab-0'),
        html.Div(id='tabs-content') ], # App content
        width={'size': 10, 'offset': 1}), # add left and right margin
    ]) ], fluid=True )

# -----------------------------------------------------------------------------------------------------------------------
@callback( Output('tabs-content', 'children'), Input('tabs', 'value') )

def render_content(tab):
    if tab == 'tab-0': # ================================================================================================
        return html.Div([ html.Br(),
            html.Div(['Overview of the data available from the', html.Span(' ALSPAC ', style={'font-weight':'bold'}), 'cohort.']),
            html.Hr(),
            dcc.Graph(id='overview-fig', figure = plot_overview()),
        ])
        
    elif tab == 'tab-1': # ==============================================================================================
        return  html.Div([ html.Br(),
            html.Div(['Results of the generalized', bold_it('cross-lag panel model'), 
                      'described as model 1 in the paper.',wrap_it(2),'Using the selection pane below, you can decide which depression report (i.e., self or \
                       parental reports) and cardio-metabolic risk (CMR) marker you want to model. You can then inspect the variables included in the \
                       model by clicking on the inspect icon or on the graph nodes directly. Check the table on the right for info on the model fit.',
                       wrap_it(), undeline_it('Note'),': by default, the "classic" cross-lag panel model is presented, but the parameter conbination can be constumized \
                       using the tickboxes on the right (don\'t forget to hit the ', badge_it('Update model', 'grey'),' button to see the changes). Hit the ',
                       badge_it('Best fit', 'silver'),' button to display the best fitting model (i.e. lowest AIC).']),
            html.Hr(),
            # Input 
            dbc.Row([dbc.Col([
                         html.H5(children='Depression score', style={'textAlign':'left'}),
                         dcc.RadioItems(id='dep-selection',
                                        options = DEPmarker_checklist(), value='sDEP', 
                                        inputStyle={'margin-left':'20px','margin-right':'20px'}),
                         html.Br(),
                         html.H5(children='Cardio-metabolic marker', style={'textAlign':'left'}),
                         dcc.Dropdown(id='cmr-selection', 
                                      options=CMRmarker_checklist(),
                                      value='FMI') ], width={'size': 5, 'offset': 1}), 
                      dbc.Col([
                         html.H5(children='Model estimation', style={'textAlign':'left'}),
                         param_checklist('sDEP', 'FMI', p='lt'),
                         param_checklist('sDEP', 'FMI', p='ma'),
                         html.Div( [ dbc.Button('Best fit', id='bestfit-button', color='secondary', n_clicks=0,
                                               style={'font-weight':'bold', 'background-color':'silver','padding':'4px 10px'}), 
                                     dbc.Button('Update model', id='update-button', color='secondary', n_clicks=0,
                                               style={'font-weight':'bold', 'background-color':'grey','padding':'4px 10px',
                                                      'margin-left': '15px'})], 
                                  style={'width':'30%','height':'35%','float':'right'}),
                         html.Div( id='failed-model', style={'color':'red', 'width':'60%','height':'25%','float':'left'}),
                     ], width={'size': 5, 'offset': 1}), 
                    ]),
            html.Hr(),
            # Time plot 
            html.Div([ # html.I(className="bi bi-info-circle-fill me-2", style={'color':'black'}),
                dbc.Accordion([ dbc.AccordionItem([ dcc.Graph(id='time-graph', figure = make_plot1('sDEP','FMI'))],
                                                       title='Inspect the variables included in this model')], start_collapsed=True)]),
            dbc.Row([
                # Network
                dbc.Col([cyto.Cytoscape(id='cyto-graph',
                                    layout={'name': 'preset', 'fit':False},
                                     style={'width': '100%', 'height': '1000px'}, minZoom=1, maxZoom=1, # reduce the range of user zooming 
                                  elements=make_net1('sDEP', 'FMI'), 
                                stylesheet=stylenet1)], width=9), 
                 # Table
                dbc.Col([ html.Br(), html.Div(id='fitm-table', children=[dbc.Table.from_dataframe(df = make_table1('sDEP', 'FMI'), 
                                                                         # style={'text-align':'right','width':'25%'},
                                                                         color='light', striped=True, bordered=True, hover=True, size='lg')]) ],
                        width=3) ]), 
                          
                # Pop variable descriptives
                dbc.Offcanvas(id='pop1', children=[ dcc.Graph()], title='', is_open=False, placement='end', style={'width': 700})
        ])
    
    elif tab == 'tab-2': # ==============================================================================================
        return html.Div([ 
            html.Br(),
            html.Div(['Results of the', bold_it('cross-lag network analyis'), 'performed using the variables listed below.']),
            html.Hr(),
            dbc.Row([
                dbc.Col([ # Temporal network
                    html.H4('Temporal (within-person) network'),   
                    cyto.Cytoscape(id='temp-net', layout={'name':'preset'},
                          style={'width':'30%', 'height':'60%', 'position':'absolute', 'left':150, 'top':350, 'z-index':999},
                          minZoom=1, maxZoom=1, # reduce the range of user zooming 
                          elements = make_net2('t'), 
                          stylesheet= make_netstyle2('t') ) ], width=4),
                dbc.Col([ # Contemporaneous network
                    html.H4('Contemporaneous (within-person) network', style={'textAlign':'center'}),
                    cyto.Cytoscape(id='cont-net', layout={'name':'preset'},
                          style={'width':'30%', 'height':'60%', 'position':'absolute', 'left':850, 'top':350, 'z-index':999},
                          minZoom=1, maxZoom=1, # reduce the range of user zooming 
                          elements = make_net2('c'), 
                          stylesheet= make_netstyle2('c') )], width=4),
                dbc.Col([ # Between person network
                    html.H4('Contemporaneous (between-person) network', style={'textAlign':'right'}),
                    dash_table.DataTable(id='temp-ci-tab', columns=[ {'name': i, 'id': i} for i in make_tab2('t') ],
                                               sort_action='custom', sort_mode='single', sort_by=[], 
                                               fixed_columns={'headers': True, 'data': 1}, # Fix node name column 
                                               # style_as_list_view=True, # Remove vertical lines between columns 
                                               style_header={'fontWeight':'bold'},
                                               style_cell={'fontSize':20, 'font-family':'sans-serif'},
                                               style_cell_conditional=[{'if': {'column_id': 'Node'}, 'width': '250px'}],
                                               style_data={'whiteSpace':'normal', 'height': 'auto','lineHeight':'20px', 
                                                           'minWidth': '100px', 'width': '100px', 'maxWidth': '100px'}, 
                                               style_table={'overflowX': 'auto', 'minWidth': '100%'})
                    # cyto.Cytoscape(id='betw-net', layout={'name':'preset'},
                    #       style={'width':'30%', 'height':'60%', 'position':'absolute', 'left':1550, 'top':350, 'z-index':999},
                    #       minZoom=1, maxZoom=1, # reduce the range of user zooming 
                    #       elements = make_net2('b'), 
                    #       stylesheet= make_netstyle2('b') )
                ], width=4)
            ]),
            # dbc.Row([
            #     dbc.Col([dash_table.DataTable(id='temp-ci-tab', columns=[ {'name': i, 'id': i} for i in make_tab2('t') ],
            #                                    sort_action='custom', sort_mode='single', sort_by=[], 
            #                                    fixed_columns={'headers': True, 'data': 1}, # Fix node name column 
            #                                    # style_as_list_view=True, # Remove vertical lines between columns 
            #                                    style_header={'fontWeight':'bold'},
            #                                    style_cell={'fontSize':20, 'font-family':'sans-serif'},
            #                                    style_cell_conditional=[{'if': {'column_id': 'Node'}, 'width': '250px'}],
            #                                    style_data={'whiteSpace':'normal', 'height': 'auto','lineHeight':'20px', 
            #                                                'minWidth': '100px', 'width': '100px', 'maxWidth': '100px'}, 
            #                                    style_table={'overflowX': 'auto', 'minWidth': '100%'})])
            # ])
             
        ])
     
    elif tab == 'tab-3': # ==============================================================================================
        return html.Div([
            html.Br(), 
            html.Span(['Results of each', bold_it('cross-sectional network model'), 'conducted as follow-up analyses. \
            Please select the timepoint of interest from the slides to visualize the network structure and correspoding centrality indices.']),
            html.Br(), html.Hr(),
            # Slider
            html.Span(children=[ html.Div('Select a timepoint:'), html.Br(style={'line-height':'5'}) ]),
            dcc.Slider(id='cros-net-slider', min=9.7, max=24.2, step=None, value=9.7, marks=timemarks3, included=False ),
            dbc.Row([
                # Network
                dbc.Col([cyto.Cytoscape(id='cros-net', layout={'name':'preset'},
                          style={'width':'40%', 'height':'80%', 'position':'absolute', 'left':150, 'top':380, 'z-index':999},
                          minZoom=1, maxZoom=1, # reduce the range of user zooming 
                          elements = make_net3(9.7)[0], 
                          stylesheet= stylenet3)], width=6), 
                # Table
                dbc.Col([ html.Br(),
                          dash_table.DataTable(id='ci-table', columns=[ {'name': i, 'id': i} for i in make_net3(9.7)[1].columns ],
                                               sort_action='custom', sort_mode='single', sort_by=[], 
                                               fixed_columns={'headers': True, 'data': 1}, # Fix node name column 
                                               # style_as_list_view=True, # Remove vertical lines between columns 
                                               style_header={'fontWeight':'bold'},
                                               style_cell={'fontSize':20, 'font-family':'sans-serif'},
                                               style_cell_conditional=[{'if': {'column_id': 'Node'}, 'width': '250px'}],
                                               style_data={'whiteSpace':'normal', 'height': 'auto','lineHeight':'20px', 
                                                           # 'minWidth': '100%',},
                                                           'minWidth': '100px', 'width': '100px', 'maxWidth': '100px'}, 
                                               style_table={'overflowX': 'auto', 'minWidth': '100%'})
                           ], width={'size': 5, 'offset': 1}),
                 # Pop variable descriptives
                 dbc.Offcanvas(id='pop3', children=[ dcc.Graph()],
                                title='Lab name', is_open=False, placement='end', style={'width': 700}), # backdrop
            ]) 
            
        ])

# Control variable selection (if mother reports are selected, only some CMR are available) 
@callback(
    Output('cmr-selection', 'options'),
    Output('dep-selection', 'options'),
    Input('dep-selection', 'value'),
    Input('cmr-selection', 'value')
)
def update_dropdown_options(dep_selection, cmr_selection):
    
    if cmr_selection not in ['FMI','LMI','BMI','waist_circ','total_fatmass','total_leanmass']:
        return no_update, DEPmarker_checklist(disable_maternal=True)
        
    if dep_selection=='sDEP':
        return CMRmarker_checklist('s'), DEPmarker_checklist()
        
    return CMRmarker_checklist('m'), DEPmarker_checklist()


# Control variable selection, updating the scatterplot and the ticks on the model selection pane
@callback(
    Output('time-graph', 'figure'),
    Input('dep-selection', 'value'),
    Input('cmr-selection', 'value')
)
def update_time_plot(dep_selection, cmr_selection):
    # Prevent update if the selection is not allowed (i.e., combination of markers not available) 
    if dep_selection=='mDEP' and cmr_selection not in ['FMI','LMI','BMI','waist_circ','total_fatmass','total_leanmass']:
        return no_update
    # Default output: only long term paramenters
    return make_plot1(dep_selection, cmr_selection)
    

# Based on variable and model selection display graph (or notify if model didn't converge)
@callback(
    Output('cyto-graph', 'elements'),
    Output('fitm-table', 'children'),
    Output('failed-model','children'),
    Output('lt-checklist', 'value'),
    Output('ma-checklist', 'value'),
    
    Input('dep-selection', 'value'),
    Input('cmr-selection', 'value'),
    Input('update-button', 'n_clicks'),
    Input('bestfit-button', 'n_clicks'),
    State('lt-checklist', 'value'),
    State('ma-checklist', 'value') # prevent_initial_call=True
)

def update_graph(dep_selection, cmr_selection, click_manual, click_best, lt_checklist, ma_checklist):
    
    # Prevent update if the selection is not allowed (i.e., combination of markers not available) 
    if dep_selection=='mDEP' and cmr_selection not in ['FMI','LMI','BMI','waist_circ','total_fatmass','total_leanmass']:
        return no_update, no_update, no_update, no_update, no_update
    
    if ctx.triggered_id == 'update-button':
        
        checked = lt_checklist + ma_checklist
        series = pd.Series([ 1 if x in checked else 0 for x in model_structure.index], index=model_structure.index)
        retrieve_model = model_structure.columns[ model_structure.eq(series, axis=0).all() ]
        
        if len(retrieve_model) == 0: 
            return no_update, no_update, 'Sorry, I did not estimate this model.', no_update, no_update
            
        model_name = retrieve_model[0]
        
        update_graph = make_net1(dep_selection, cmr_selection, which_model=model_name)
        
        if update_graph == 'fail': # prevents any single output updating
             return no_update, no_update, 'Sorry, this model did not converge.', no_update, no_update

        tab_df = make_table1(dep_selection, cmr_selection, which_model=model_name)
        update_tab = dbc.Table.from_dataframe(df = tab_df, color='light', striped=True, bordered=True, hover=True, size='lg')
        
        return update_graph, update_tab, None, no_update, no_update 
    
    elif ctx.triggered_id == 'bestfit-button':
        
        update_graph = make_net1(dep_selection, cmr_selection, which_model='best')
        
        tab_df = make_table1(dep_selection, cmr_selection, which_model='best')
        update_tab = dbc.Table.from_dataframe(df = tab_df, color='light', striped=True, bordered=True, hover=True, size='lg')
        
        return update_graph, update_tab, None, best_fit1(dep_selection, cmr_selection, 'lt'), best_fit1(dep_selection, cmr_selection, 'ma')

    
    update_graph = make_net1(dep_selection, cmr_selection)

    if update_graph == 'fail': # prevents output updating
        return no_update, no_update, 'Sorry, this model did not converge.', no_update, no_update
    
    tab_df = make_table1(dep_selection, cmr_selection)
    update_tab = dbc.Table.from_dataframe(df = tab_df, color='light', striped=True, bordered=True, hover=True, size='lg')

    return make_net1(dep_selection, cmr_selection), update_tab, None, ['ltCL_dep','ltCL_cmr','ltAR_dep','ltAR_cmr'], []

# Dispaly info upon tapping on node
@callback(
    # Output('display-labels', 'children'),
    Output('pop1', 'is_open'),
    Output('pop1', 'title'),
    Output('pop1', 'children'),
    Input('cyto-graph', 'tapNodeData'),
    [State('pop1', 'is_open')],
) 
def displayTapNodeData1(node, is_open):
    if node:
        titl = node['label']
        
        if titl=='Eta':
            return is_open, no_update, no_update
            
        plot = dcc.Graph(figure=desc_plot(titl))
        
        return not is_open, titl, plot
        
    return is_open, no_update, no_update

# TAB 3: display cross-sectional network based on selected timpoint ------------------------------------------------
@callback(
    Output('cros-net', 'elements'),
    Output('ci-table', 'data'),
    Input('cros-net-slider', 'value'),
    Input('ci-table', 'sort_by')
)
def update_crosnet(timepoint, sort_by):
    
    if ctx.triggered_id == 'ci-table':
        
        dff = make_net3(timepoint)[1].sort_values( sort_by[0]['column_id'],
                                                  ascending=sort_by[0]['direction'] == 'desc', inplace=False )
        
        return no_update, dff.to_dict('records')
    
    return make_net3(timepoint)[0], make_net3(timepoint)[1].to_dict('records')

# Dispaly info upon tapping on node
@callback(
    # Output('display-labels', 'children'),
    Output('pop3', 'is_open'),
    Output('pop3', 'title'),
    Output('pop3', 'children'),
    Input('cros-net', 'tapNodeData'),
    [State('pop3', 'is_open')],
) 
def displayTapNodeData3(node, is_open):
    if node:
        titl = node['label']
        plot = dcc.Graph(figure=desc_plot(node['id']))
        return not is_open, titl, plot
        
    return is_open, no_update, no_update

if __name__ == '__main__':
    app.run(debug=True, jupyter_mode="external")

