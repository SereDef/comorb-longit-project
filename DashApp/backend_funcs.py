
import os, pyreadr, json, re, textwrap
import pandas as pd
import numpy as np

from plotly.subplots import make_subplots
import plotly.graph_objects as go

# ---- General set-ups ----------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# Set the path to results, to be fixed at deployment
PATH = './assets/results/'

# Opening Descriptives JSON file
with open('./assets/descrip.json') as jf:
    desc = json.load(jf)

# other descriptives file
all_desc = pd.read_csv('./assets/descdf.csv', index_col=0)

# Matrix with parameter combinations (8x57) used to fit the models in Rscript 1.
model_structure = pd.read_csv('./assets/model_structure.csv').set_index('Unnamed: 0')

labels = {'DEP01':'Felt miserable or unhappy',
          'DEP02':'Didn\'t enjoy anything at all',
          'DEP03':'Felt so tired they just sat around and did nothing',
          'DEP04':'Was very restless',
          'DEP05':'Felt they were no good any more',
          'DEP06':'Cried a lot',
          'DEP07':'Found it hard to think properly or concentrate',
          'DEP08':'Hated themselves',
          'DEP09':'Felt they were a bad person',
          'DEP10':'Felt lonely',
          'DEP11':'Thought nobody really loved them',
          'DEP12':'Thought they would never be as good as other people',
          'DEP13':'Felt they did everything wrong',
          'DEP_score':'Total depression score',
          'height':'Height',
          'weight':'Weight',
          'BMI':'Body mass index',
          'waist_circ':'Waist circumference',
          'waist_hip_ratio':'Waist/hip ratio',
          'total_fatmass':'Total fat mass',
          'total_leanmass':'Total lean mass',
          'trunk_fatmass':'Trunk fat mass',
          'android_fatmass':'Android fat mass',
          'FMI':'Fat mass index',
          'LMI':'Lean mass index',
          'TFI':'Trunk fat mass index',
          'liver_fat':'Liver fat',
          'SBP':'Systolic blood pressure',
          'DBP':'Diastolic blood pressure',
          'PWV':'Pulse wave velocity',
          'IMT':'Intima-media thickness',
          'heart_rate':'Heart rate',
          'LVM':'Left ventricular mass',
          'RWT':'Relative wall thickness',
          'FS':'Fractional shortening',
          'tot_chol':'Total cholesterol',
          'HDL_chol':'HDL-cholesterol',
          'LDL_chol':'LDL-cholesterol',
          'insulin':'Insulin',
          'triglyc':'Triglycerides',
          'glucose':'Glucose',
          'CRP':'C-reactive protein',
          'IL_6':'Interleaukin-6',
          'alcohol':'Alcohol consumption',
          'canabis':'Cannabis consumption',
          'smoking':'Smoking (tobacco)'
         }

def get_label(var): 
    # Remove year from var name if present
    vs = var.split('_')
    if bool(re.match('[0-9]', vs[-1])): 
        var_name = '_'.join(vs[:-1])
    else: var_name = '_'.join(vs)
    # Remove self or maternal report from depression item 
    if 'DEP' in var_name: var_name = var_name[1:]
    # get label
    lab = labels[var_name]
    return lab


def desc_plot(var):
    # Read in descriptives from file 
    d = desc[var]
    sum, count, dens = (pd.DataFrame.from_dict(x) for x in d)

    # Initialize figure 
    fig = go.Figure()

    if bool(re.match(r'.DEP[0-9]+|alcohol|canabis|smoking',var)):
        
        likert = ['Not true','Sometimes','True']
        colors = ['lightblue','salmon','crimson']
        
        x_lab = '<br>'.join(textwrap.wrap(get_label(var), width=20)) # Wrap label text 

        # Stacked histogram 
        for i, c in enumerate(count['count'].iloc[:-1]):
            
            fig.add_trace(go.Bar(x=[0], y=[c], name=likert[i], marker_color=colors[i], opacity=0.6, 
                                 marker_line_width=1.5, marker_line_color='white',
                                 customdata =  [count['prop_noNA'].iloc[i]*100],
                         hovertemplate = """<br> n = <b>%{y}</b> (%{customdata:.1f}%) <br><extra></extra>"""))
            
            fig.update_layout(barmode='stack', autosize=False, width=400, height=500, margin=dict(l=20, r=20, t=20, b=20),
                              yaxis=dict(title_text='Count'),
                              xaxis=dict(ticktext=['',x_lab,''], tickvals=[-0.5, 0, 0.5], tickfont=dict(size=15)) ) 
        return fig

    else: 
        x_lab = get_label(var)
        
        # Distribution plot 
        fig.add_trace(go.Scatter(x=dens['x'], y=dens['dens'], fill='tozeroy', mode='lines', line_color='salmon'))
        # add mean and median 
        fig.add_vline(x=sum['mean'][0], line_width=1.5, line_color='crimson', name='meanvalue')
        fig.add_vline(x=sum['50%'][0], line_width=1, line_dash='dash', line_color='crimson', name='medianvalue')
        # Adjust x-axis 
        fig.update_xaxes(title_text=x_lab, range=[round(sum['min']),round(sum['max'])], 
                        ticks='outside', showline=True, linecolor='black', gridcolor='lightgrey')
        fig.update_layout(autosize=False, width=600, height=300, margin=dict(l=20, r=20, t=20, b=20),
                          yaxis=dict(showticklabels=False, showline=True, linecolor='black'),   
                         )
        
        return fig

# ---- Tab 0: Data overview -----------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

age_desc = all_desc[list(all_desc.columns[all_desc.columns.str.contains('age')])]

def plot_overview(width=1900, height=1000, exclude_cmr = ['BMI','FMI','LMI','TFI','alcohol','canabis','smoking']):
    '''Input: dimentions (optional).
       Creates an interactive scatterplot figure with time of measurement on the x-axis and all measures on the y-axis. 
       Markers indicate the median age at measurement [+/- 95% CI in age range]. These values + number of observations are also displayed when hovering over the points.
    '''
    fig = make_subplots(2,1, vertical_spacing=0.05, row_heights=[13*2, 25*2])

    colors = {'sDEP':'71, 107, 237', 'mDEP':'148, 103, 189','CMR':'214, 39, 40'} # blue, purple and red 
    legend = {'sDEP':'Depression\nself-reports', 'mDEP':'Depression\nmaternal reports','CMR':'Cardio-metabolic\nrisk markers'}
    
    def scat(vartype):
        # Get ages of measurement from age_desc 
        df = age_desc[list(age_desc.columns[age_desc.columns.str.contains(vartype)])]
        
        if vartype=='CMR':
            # Extract variable names from labels (only keeping those of interest)
            cmrkeys = [i for i in list(labels.keys())[14:] if i not in exclude_cmr]
            
            # Specify their label (including manually combining some measures into one row) 
            cmrlab = {k:labels[k] for k in cmrkeys if k in labels}
            cmrlab['weight'] = 'Weight / BMI'
            cmrlab['total_fatmass'] = 'Total fat mass / FMI'
            cmrlab['total_leanmass'] = 'Total lean mass / LMI'
            
            names = list(cmrlab.values())
        else:
            names = list(labels.values())[0:13] # depression items 1-13

        # Create grid of values and ranges
        x, y = [i.flatten() for i in np.meshgrid(df.loc['50%'], names)]
        x_min,_ = np.meshgrid(df.loc['5%'], names)
        x_max,_ = np.meshgrid(df.loc['95%'], names)

        # The CMR grid will need to contain some "NA" (and so do the counts)
        if vartype=='CMR': 
            ages = [i[-1][:-1] for i in df.columns.str.split('_')] # get age without 'y'
    
            cmrmat = pd.DataFrame(columns=ages, index=cmrlab.keys())
            counts = cmrmat.copy()
            
            for age in ages:
                # what measures are available at that age
                t = all_desc[list(all_desc.columns[all_desc.columns.str.contains(age) & ~(all_desc.columns.str.contains('DEP|age'))])]
                meas = [x.split('_'+age)[0] for x in t.columns if x.split('_'+age)[0] not in exclude_cmr]
                nobs = [  t.loc['count', x] for x in t.columns if x.split('_'+age)[0] not in exclude_cmr]
                
                cmrmat.loc[meas, age] = [cmrlab[m] for m in meas]
                counts.loc[meas, age] = [int(c) for c in nobs]
                
            y = list(cmrmat.stack(future_stack=True))
            counts = list(counts.stack(future_stack=True))
        
        else: 
            counts = []
            for i in [vartype[0] + s for s in list(labels.keys())[0:13] ]:
                counts.extend( [int(c) for c in list(all_desc.loc['count', list(all_desc.columns[all_desc.columns.str.contains(i)])])] )
            
        p = go.Scatter(x = x, y = y, mode='markers',
                       marker = dict(size = 10, symbol = 'square', color=f'rgb({colors[vartype]})', opacity = .9),
                       error_x = dict(type='data', symmetric=False,
                                      array= x_max.flatten() - x,
                                      arrayminus= x - x_min.flatten(),
                                      color=f'rgba({colors[vartype]}, 0.3)', # use this for opacity 
                                      thickness=10, width=0),
                       name=legend[vartype],
                       text=[f'{round(i,1)}-{round(j,1)}' for i,j in zip(x_min.flatten(), x_max.flatten() )], customdata = [f'{i}' for i in counts],
                       hovertemplate = """ <b>%{y}</b> <br> %{x:.1f} [%{text}] years <br> <i>N =</i> %{customdata} <br><extra></extra>""" )
        return p
    
    fig.add_trace(scat('sDEP'), 1,1)
    fig.add_trace(scat('mDEP'), 1,1)
    fig.add_trace(scat('CMR'), 2,1)
    
    cmrgroups = {'Anthropometry':[4,'lightblue'], 
                 'Fat distribution':[5,'yellow'], 
                 'Arteries':[4,'orange'], 
                 'Heart':[4,'red'], 
                 'Blood-based    <br>metabolic    <br>markers':[6,'brown'], 
                 'Inflammation':[2,'blue']}
    
    cmrgrouplabels = []; y0 = -0.5
    
    for g in cmrgroups:
        y1 = y0 + cmrgroups[g][0]
        cmrgrouplabels.append(dict(x0 = -0.125, x1 = -0.13, xref='paper',
                              y0 = y0, y1 = y1-0.25, yref='y2', 
                              type='rect', fillcolor=cmrgroups[g][1], opacity=.3, line_width=0,
                              label=dict(text=f'<b>{g}    </b>', textposition='top right', padding=2)
        ))
        y0 = y1 + 0.05
        
    
    fig.update_xaxes(range=[9,26], mirror=True, ticks='outside',linecolor='black', gridcolor='lightgrey', 
                     tickmode = 'linear',tick0 = 9, dtick = 1)
    
    fig.update_yaxes( mirror=True, ticks='inside',linecolor='black') # autorange='reversed',
    
    fig.update_layout(autosize=False, width=width, height=height,
                      yaxis1 = dict(range=(12.7, -0.7)), yaxis2 = dict(range=(24.7,-0.7),), xaxis2= dict(title='Child age (years)'),
                      shapes=cmrgrouplabels, plot_bgcolor='white', margin=dict(l=20, r=20, t=20, b=20), 
                      legend=dict(orientation='h', entrywidth=300, font=dict(size=16),
                                  yanchor='bottom', y=1.02, xanchor='left', x=0) )
    
    return fig

# ---- Tab 1: generalized cross-lagged panel model ------------------------------------------------
# -------------------------------------------------------------------------------------------------

def read_res1(depname, cmrname, path=PATH+'mod1/'):
    '''Input: names of the depression report (sDEP = self or mDEP = parental reports) and cardio-metabolic risk (CMR) marker. 
       Open the .RData file created by Rscript 1.GCLPM (one for dep-cmr marker pair). This contains the following elements:
       - summ: a summary dataframe (with information about which marker is used and the timepoints included + mean ranges and number of observations)
       - fit_meas: fit measures for every parameter combination, when the model converged.
       - estimates: (unstandardized) estimates (+ robust SE, pvalues and CIs) for every parameter combination, when the model converged.
       - failed: list of models that did not converge, with corresponding error or warning message.
       Use: summ, fitm, esti, fail = read_res1('sDEP','FMI') -OR- summ = read_res1('sDEP','BMI')[0]
    '''
    res = pyreadr.read_r(f'{path}{depname}_{cmrname}.RData')
    
    summ = res['dat_summ']
    fitm = res['fit_meas'].T
    esti = res['estimates'].set_index('rep(f, nrow(es))')
    esti['sign'] = [0 if l <= 0 <= u else 1 for l, u in zip(esti['ci.lower'],esti['ci.upper'])] # based on CI because pvalue sometimes NA
    fail = list(res['failed'].index) # TODO: report warning message ...
    
    return(summ, fitm, esti, fail)


def best_fit1(depname, cmrname, list1=None):
    '''Input: names of the depression report (sDEP = self or mDEP = parental reports) and cardio-metabolic risk (CMR) marker.
       By default, returns a dataframe with the model name (indicating the excluded parameters) and structure (0,1...).  
       When list1 is provided, returns the list of long-term (list1='lt') or short-term (list1='ma') parameters estimated in the model.
    '''    
    fitm = read_res1(depname, cmrname)[1]
    
    # Best fitting model (lowest AIC)
    mod = fitm.index[fitm.aic == fitm.aic.min()][0]
    if list1: # Return list of parameters estimated in the model 
        return list(model_structure.index[(model_structure[mod] > 0) & (model_structure.index.str.contains(list1))])
    else: # Return a dataframe with its name and model structure
        return pd.DataFrame(model_structure[mod])


def make_plot1(depname, cmrname):
    '''Input: names of the depression report (sDEP = self or mDEP = parental reports) and cardio-metabolic risk (CMR) marker.
       Creates an interactive scatterplot figure with time of measurement on the x-axis and values of dep and cmr markers on the (double) y-axis. 
       Markers indicate the median [IQR] of each measure. Median values and timepoint are also displayed when hovering over the points.
    '''
    # load summary dataframe
    summ = read_res1(depname,cmrname)[0]
    
    # extract timepoints
    t_dep = [ float(x.split('_')[-1][:-1]) for x in summ.columns[:summ.shape[1]//2] ]
    t_cmr = [ float(x.split('_')[-1][:-1]) for x in summ.columns[summ.shape[1]//2:] ]

    # scatterplot function 
    def scat(t, name, fullname, shortname):
        means = summ.loc['Median', summ.columns.str.contains(name)]
        p = go.Scatter(x = t, y = means, 
                       error_y = dict(type='data', symmetric=False, # visible=True,
                                      array = summ.loc['3rd Qu.', summ.columns.str.contains(name)] - means,
                                      arrayminus = means - summ.loc['1st Qu.', summ.columns.str.contains(name)]),
                       name = fullname, text = [f'{shortname} {n}' for n in range(1,len(t)+1)],
                       marker = dict(size = 10, symbol = 'square',opacity = .8), opacity = .7,
                       hovertemplate = """ <b>%{text}</b> <br> Median: %{y:.2f} <br> Timepoint: %{x} years <br><extra></extra>""")
        return p

    fig = make_subplots(specs=[[{'secondary_y': True}]]) # Specify a double y axis to allow for different scale of depression and CMR
    
    fig.add_trace( scat(t_dep, depname, 'Depression score', 'DEP'), secondary_y=False)
    fig.add_trace( scat(t_cmr, cmrname, cmrname, cmrname), secondary_y=True )

    # Set y-axes
    def yrange(name):
        sub = summ.filter(like=name)
        ymin = sub.min(axis=1)['1st Qu.']; ymax = sub.max(axis=1)['3rd Qu.']
        range = ymax-ymin
        y_max_lower = ymax
        ymin = ymin - (range/10); ymax = ymax + (range/10)
        return [ymin, ymax, y_max_lower]
    
    fig.update_yaxes(title_text='<b>Depression score</b>', secondary_y=False, range=yrange(depname)[:2], 
                     mirror=True, ticks='outside', showline=True, linecolor='black', gridcolor='lightgrey')
    fig.update_yaxes(title_text=f'<b>{cmrname}</b>', secondary_y=True, range=yrange(cmrname)[:2], 
                     mirror=True, ticks='outside', showline=True, linecolor='black', gridcolor='lightgrey')
    # Set x-axis 
    fig.update_xaxes(title_text='Years', mirror=True, ticks='outside', showline=True, linecolor='black', gridcolor='lightgrey')

    # Group "cross-sectional" timepoints using grey background
    crosspoints = []
    for i in range(len(t_dep)):
        xmin = min(t_dep[i], t_cmr[i])-.2; xmax = max(t_dep[i], t_cmr[i])+.2
        # rectangles 
        crosspoints.append( dict(x0 = str(xmin), x1 = str(xmax), y0 = 0, y1 = 1, xref='x', yref='paper', 
                                 type='rect', fillcolor='lightgray', opacity=.3, line_width=0, layer='below') )  
        # text 
        fig.add_trace( go.Scatter(x=[xmin + (xmax-xmin)/2], y=[yrange(depname)[2]], mode='text', text=i+1, 
                                  textposition='top center',
                                  textfont_size=13, textfont_color='dimgray', showlegend=False) )
    # Background
    fig.update_layout(# title = dict(text='Included measures\n', font=dict(size=15), automargin=True, yref='paper'),
                      plot_bgcolor='white', shapes=crosspoints, margin=dict(l=20, r=20, t=20, b=20))
    
    return(fig)


def make_net1(depname, cmrname, which_model = 'maCL_dep-maCL_cmr-maAR_dep-maAR_cmr', width=1000):
    '''Input: names of the depression report (sDEP = self or mDEP = parental reports) and cardio-metabolic risk (CMR) marker; model structure and size of the graph.
       Creates a list of dictionaries to use as input for the elements arg of cytoscape graph. 
       This only creates the core structure, stlyle parameters are defined in stylesheet1.
    '''
    # read data
    summ,_,esti,fail = read_res1(depname, cmrname)
    
    # Best fitting model
    if which_model == 'best': modstr = best_fit1(depname, cmrname).T
    elif which_model in fail: 
        return 'fail' # modstr = best_fit(depname, cmrname).T
    else: modstr = pd.DataFrame(model_structure[which_model]).T
    
    # Extract the estimated paramters from the result files (returns estimates and significance [= '*' or '']
    def extr_est(name=None, which=None, lamb=False, eta_corr=False, imp_corr=False, model_output=esti):
        df = model_output.loc[modstr.index[0]]
        if   lamb:     e,s = df.loc[(df.lhs==f'eta_{name}')&(df.op=='=~'), ['est','sign']].round(2).iloc[which]
        elif eta_corr: e,s = df.loc[(df.lhs=='eta_dep') & (df.rhs=='eta_cmr'),['est','sign']].round(2).iloc[0]
        elif imp_corr: e,s = df.loc[ df.label.str.contains('comv'), ['est','sign']].round(2).iloc[which]
        else:          e,s = df.loc[df.label.str.contains(name)].iloc[::-1][['est','sign']].round(2).iloc[which]
        return(e,s)

    # Ready to draw
    nt = summ.shape[1]//2 # Number of timepoints
    
    # Get timepoints of measurement from summary
    timedic = dict(zip([f'dep{i+1}' if 'DEP' in n else f'cmr{i-summ.shape[1]//2+1}' for i,n in enumerate(summ.columns)],
                       [float(i.split('_')[-1][:-1]) for i in summ.columns]))
    
    pos_top = 30; pos_bot = 550 # Vertical cohordinates (in pixel)
    rshift = 90; obs_pos = 90; imp_pos=190 # Distances from left edge and top (in pixels) 
    vs = ['dep','cmr']
    
    e = [] # Initialize

    e.append({'data': {'id':f'name_dep', 'label':'Depression\nscore'}, 'classes':'notes', 'position':{'x':120,'y':pos_top+obs_pos}})
    e.append({'data': {'id':f'name_cmr', 'label':'\n'.join(textwrap.wrap(get_label(cmrname), width=10))}, 'classes':'notes', 
              'position':{'x':120,'y':pos_bot-obs_pos}})
   
    for eta, pos in enumerate([pos_top, pos_bot]):
        # Eta factors nodes
        e.append({'data': {'id':f'eta{eta}', 'label':'Eta'}, 'classes':'latent', 'position':{'x':width/2 +rshift,'y':pos}})
        
    
    #  Eta factors correlations 
    e.append({'data': {'source':'eta0', 'target':'eta1', 'firstname':'eta_corr', 'label':'%.2f' % extr_est(eta_corr=True)[0] }})
    
    for i in range(1,nt+1): 
        
        e.append({'data': {'source':f'imp_dep{i}', 'target':f'imp_cmr{i}', 'firstname':'imp_corr', 
                           'label':'%.2f' % extr_est(which=i-1, imp_corr=True)[0] }})
        
        for eta, v in enumerate(vs):
            
    # ===== Other nodes
            p = [pos_top+obs_pos, pos_top+imp_pos] if v=='dep' else [pos_bot-obs_pos, pos_bot-imp_pos] # define vertical position
            e.extend([
                # Observed variables
                {'data': {'id':f'{v}{i}','firstname':f'{v.upper()}\n{timedic[v+str(i)]} yrs', 'label':summ.columns[(i-1)+(nt*eta)]},
                 'classes':'observed',
                 'position' : {'x':((width/nt)*i)-(width/nt)/2 +rshift, 'y': p[0]} },
                # Impulses
                {'data': {'id':f'imp_{v}{i}', 'label': f'impulse {i}'}, 
                 'classes':'latent',
                 'position' : {'x':((width/nt)*i)-(width/nt)/2 +rshift, 'y': p[1]} }
            ])
        
    # ===== Edges: lambdas
            e.append({'data': {'source':f'eta{eta}', 'target':f'{v}{i}', 'firstname':'lambda',
                               'label':'%.2f' % extr_est(f'{v}',i-1, lamb=True)[0] }})
            # impulses link
            e.append({'data': {'source':f'imp_{v}{i}', 'target':f'{v}{i}', 'firstname':'imp_link'}})
            
            if i < nt: 
                otherv = abs(eta-1)
                # maAR and AR terms
                if modstr[f'ltAR_{v}'][0]: e.append({'data': {'source':f'{v}{i}', 'target':f'{v}{i+1}',
                                   'weight':'%.2f' % (extr_est(f'^AR_{v}', i-1)[0]*(timedic[f'{v}{i+1}'] - timedic[f'{v}{i}'])), 
                                   'sign': extr_est(f'^AR_{v}', i-1)[1],                       
                                   'label': f'AR{i}', 'firstname':'direct'}})
                if modstr[f'maAR_{v}'][0]: e.append({'data': {'source':f'imp_{v}{i}', 'target':f'{v}{i+1}',
                                   'weight':'%.2f' % (extr_est(f'^maAR_{v}',i-1)[0]/(timedic[f'{v}{i+1}'] - timedic[f'{v}{i}'])), 
                                   'sign': extr_est(f'^maAR_{v}',i-1)[1],                                   
                                   'label': f'maAR{i}', 'firstname':'direct'}})
                # maCL and CL terms
                if modstr[f'ltCL_{v}'][0]: e.append({'data': {'source':f'{v}{i}', 'target':f'{vs[otherv]}{i+1}',
                                   'weight':'%.2f' % (extr_est(f'^CL_{v}', i-1)[0]*(timedic[f'{vs[otherv]}{i+1}'] - timedic[f'{v}{i}'])), 
                                   'sign': extr_est(f'^CL_{v}', i-1)[1],
                                   'label': f'CL{i}', 'firstname':'direct'}})
                if modstr[f'maCL_{v}'][0]: e.append({'data': {'source':f'imp_{v}{i}', 'target':f'{vs[otherv]}{i+1}',
                                   'weight':'%.2f' % (extr_est(f'^maCL_{v}',i-1)[0]/(timedic[f'{vs[otherv]}{i+1}'] - timedic[f'{v}{i}'])), 
                                   'sign': extr_est(f'^maCL_{v}',i-1)[1],
                                   'label': f'maCL{i}', 'firstname':'direct'}})
                      
    return e

# Define the stile of the graph 
stylenet1=[ 
    # Noting down names 
    {'selector':'.notes',
     'style':{'background-color':'white','font-size':18,
              'content':'data(label)','color':'k','font-weight':'bold','text-halign':'left','text-valign':'center','text-wrap':'wrap'}},
    
    # Nodes - shape & color
    {'selector':'.observed',
     'style':{'shape':'rectangle', 'height':40, 'width':80, 'border-width':2,'background-color':'white', 'border-color':'k', 
              'content':'data(firstname)','color':'grey','text-halign':'center','text-valign':'center','text-wrap':'wrap'}},
    
    {'selector':'.latent', 
     'style':{'shape':'round', 'height':20, 'width':20, 'border-width':1,'background-color':'white', 'border-color':'silver'}},
    
    # Edges
    {'selector':'edge[firstname *= "direct"]', # directed paths 
     'style':{'curve-style':'straight', 'target-arrow-shape':'vee', 'width': 3, 'arrow-scale':1.2 }},
    
    {'selector':'edge[firstname *= "imp_link"]', # impulses links
     'style':{'curve-style':'straight', 'target-arrow-shape':'vee', 'width': 1, 'arrow-scale':.8 }},

    # Correlations
    {'selector':'edge[firstname *= "eta_corr"]',
     'style':{'curve-style':'unbundled-bezier','target-arrow-shape':'vee','source-arrow-shape':'vee', 'width': 1, 
              'control-point-distances': [-450,-520,-530,-520,-450],'control-point-weights': [0.01, 0.20, 0.5, 0.80, 0.99],
              'label':'data(label)','font-size':15, 'text-background-color':'silver', 'text-background-opacity':.7 }},
    
    {'selector':'edge[firstname *= "imp_corr"]', 
     'style':{'curve-style':'unbundled-bezier','target-arrow-shape':'vee','source-arrow-shape':'vee','width': 1, 
              'label':'data(label)','font-size':15, 'text-background-color':'silver', 'text-background-opacity':.7 }},
    
    # Lambdas
    {'selector':'edge[firstname *= "lambda"]', 
     'style':{'curve-style':'straight', 'target-arrow-shape':'vee', 'width': 1, 'arrow-scale':.8,
              'label':'data(label)','font-size':15, 'text-background-color':'silver', 'text-background-opacity':.7 }},

    # Dashed lines 
    # {'selector':'edge[weight < 0.01][weight > -0.01]', 'style':{'line-style':'dashed'}},
    # {'selector':'edge[sign < 1]', 'style':{'line-style':'dashed'}},
]
# Set the color of each edge type and the distance between sorce and label displaying its weight (to avoid overlapping) 
d = {'AR':['red',40],'maAR':['orange',70],'CL':['green',220],'maCL':['lightblue',40]}

for c in d.keys():
    
    stylenet1.extend([{'selector':f'[label *= "{c}"][sign > 0]', 
                      'style': {'line-color':d[c][0], 'target-arrow-color':d[c][0],
                                'source-label':'data(weight)', 'source-text-offset': d[c][1],'font-size':20, 'font-weight':'bold',
                                'text-background-color':d[c][0], 'text-background-opacity':.5 }},
                      {'selector':f'[label *= "{c}"][sign < 1]', 
                      'style': {'line-style':'dashed', 'width': 1.5,
                                'line-color':d[c][0], 'target-arrow-color':d[c][0],
                                'source-label':'data(weight)', 'source-text-offset': d[c][1],'font-size':18,
                                'text-background-color':d[c][0], 'text-background-opacity':.5 }}
                     ])


def make_table1(depname, cmrname, which_model='maCL_dep-maCL_cmr-maAR_dep-maAR_cmr'):
    '''Input: names of the depression report (sDEP = self or mDEP = parental reports) and cardio-metabolic risk (CMR) marker; model structure.
       Extracts the fit measures for the specified model and stores into a table to be displayed next to the graph. 
    '''
    fitm = read_res1(depname, cmrname)[1]
    
    if which_model == 'best':
        which_model = fitm.index[fitm.aic == fitm.aic.min()][0]# Best fitting model (lowest AIC)
    
    dt = pd.DataFrame(fitm.loc[which_model, ['npar', 'df', 'chisq', 'pvalue','cfi', 'tli','rmsea','srmr','aic', 'bic', ]])
    dt = dt.rename(columns={ which_model:' '}).round(3)
    dt.insert(loc=0, column='Fit measures', value=['Number of parameters','Degrees of freedom','\u03C7\u00b2','P-value (\u03C7\u00b2)',
                                                   'CFI','TLI','RMSEA','SRMR','AIC','BIC'])

    def format_row(x, form='{0:.2f}'):
        try: return form.format(x)
        except: return x

    dt.loc[['npar','df']] = dt.loc[['npar','df']].applymap(format_row, form='{0:.0f}')
    dt.loc[['aic','bic']] = dt.loc[['aic','bic']].applymap(format_row, form='{0:.1f}')
    dt.loc[['chisq','cfi','tli','rmsea','srmr']] = dt.loc[['chisq','cfi','tli','rmsea','srmr']].applymap(format_row, form='{0:.2f}')
    dt.loc[['pvalue']] = dt.loc[['pvalue']].applymap(format_row, form='{0:.3f}')
    
    return(dt)

# ---- Tab 2: cross-lagged panel network model ---------------------------------------------------
# -------------------------------------------------------------------------------------------------

def read_res2(which_net, path=PATH+'mod2/'):
    '''Input: temporal/contemp/between person network input. 
       Open the .RData file created by Rscript 2.CLNPM. This contains the following elements:
       - fit: fit measures + number of observations the network is based on.
       - layout: the spring graphical disposition computed by `qgraph`
       - t_net/c_net/b_net: dataframe with all edge weights
       - t_cent/c_cent/b_cent:: 95% confidence intervals for those weights
       Use: fit, lay, wm, ci = read_res2(t')
    '''
    res = pyreadr.read_r(f'{path}unpruned.RData')
    
    fit = res['fit'] # [[]]
    
    # layout computed by qgraph (spring algorithm) 
    lay = res['layout'].rename(columns={0:'x_og',1:'y_og'})
    # rescale to pixels 
    lay['x'] = np.interp(lay['x_og'], (lay['x_og'].min(), lay['x_og'].max()), (0, 500))
    lay['y'] = np.interp(lay['y_og'], (lay['y_og'].min(), lay['y_og'].max()), (0, 500))
    # add class 
    lay['class'] = ['dep' if t else 'cmr' for t in lay.index.str.contains('DEP')]

    nw = res[f'{which_net}_net']
    nw['from'] = nw['from'].map({ x : lay.index[x-1] for x in nw['from'] })
    nw['to'] = nw['to'].map({ x : lay.index[x-1] for x in nw['to'] })
    nw['dir'] = ['neg' if x<0 else 'pos' for x in nw.weight]
    
    ci = res[f'{which_net}_cent']
    
    return fit, lay, nw, ci

def make_tab2(which_net):
    '''Input: network type (temporal, contemporaneous or between person).
       Creates the table with centrality indices tab.
    '''
    ci = read_res2(which_net)[3] # read in data 

    # Centrality indices tab
    ci_tab = ci.round(2)
    ci_tab.insert(1, 'Node', [get_label(name) for name in ci_tab.index]) # Replace node name with its label 

    return ci_tab
    
def make_net2(which_net):
    '''Input: network type (temporal, contemporaneous or between person). 
       Creates a network structure.
    '''
    fit,lay,nw,_ = read_res2(which_net) # read in data 
    
    # already trimmed estimates
    nodes = [{'data': {'id':node, 'label': '\n'.join(textwrap.wrap(get_label(node), width=20))}, 'classes':lay.loc[node,'class'], 
             'position':{'x':lay.loc[node,'x'], 'y':lay.loc[node,'y'] }} 
       for node in lay.index ]
    
    edges = [{'data': {'source':a, 'target':b, 'weight':round(abs(w)*10,3), 'width':round(abs(w)*10,2)}, 'classes':c} 
       for a,b,w,c in nw.itertuples(index=False)]
    
    network = nodes+edges

    return network 

def make_netstyle2(which_net):
    
    sn2 = [ {'selector': 'node', 'style': {'label': 'data(label)', 'text-wrap':'wrap'} },
           # Edge opacty and width
           {'selector': 'edge', 'style': {'opacity': 'data(weight)', 'width': 'data(width)'}},
           # Color nodes by group
           {'selector': '.dep', 'style': {'background-color': 'lightblue'} },
           {'selector': '.cmr', 'style': {'background-color': 'pink'} },
           # Color edges by positive/negative weights
           {'selector': '.neg', 'style': {'line-color': 'red', 'target-arrow-color':'red'} },
           {'selector': '.pos', 'style': {'line-color': 'blue','target-arrow-color':'blue'} } ]
    
    if which_net=='t':
        sn2 = sn2 + [{'selector': 'edge', 'style':{'curve-style':'straight', 'target-arrow-shape':'vee','arrow-scale':.8 }}]
    
    return sn2


# ---- Tab 3: cross-sectional network models ------------------------------------------------------
# -------------------------------------------------------------------------------------------------

def read_res3(time, path=PATH+'mod3/'):
    '''Input: timepoint of interest. 
       Open the .RData file created by Rscript 2.CLNPM (one for each timepoint). This contains the following elements:
       - wm: dataframe with all edge weights
       - ci: 95% confidence intervals for those weights
       - fit: fit measures + number of observations the network is based on.
       - layout: the spring graphical disposition computed by `qgraph`
       Use: wm, ci, fit, lay = read_res3('9.8y-9.8y')
    '''
    res = pyreadr.read_r(f'{path}crosnet_{time}y.RData')
    # weight matrix
    wm = res['wm']; wm['link'] = wm.index; wm[['a','b']] = wm.link.str.split(' ', expand = True)
    wm = wm.loc[wm.a!=wm.b, ] # remove links to between an edge and itself
    wm = wm.reset_index()[['a','b','V1']].rename(columns={'a':'node1','b':'node2','V1':'weight'})
    wm['dir'] = ['neg' if x<0 else 'pos' for x in wm.weight]
    # centrality indices
    ci = res['ci']; ci['class'] = ['dep' if t else 'cmr' for t in ci.node.str.contains('DEP')]
    # fit measures and number of observations 
    fit = res['fit'].T.round(3)
    # layout computed by qgraph (spring algorithm) 
    lay = res['layout'].rename(columns={0:'x_og',1:'y_og'})
    # rescale to pixels 
    lay['x'] = np.interp(lay['x_og'], (lay['x_og'].min(), lay['x_og'].max()), (0, 700))
    lay['y'] = np.interp(lay['y_og'], (lay['y_og'].min(), lay['y_og'].max()), (0, 700))
    
    return wm, ci, fit, lay


times = [file.split('_')[1][:-7] for file in os.listdir(PATH+'mod3') if file.endswith('.RData')]
times.sort() # just for readability, not necessary 

# create marks dict to use in slider, also provides a map for make_net3
timemarks3 = dict()
for t in times:
    where = round(np.mean([float(n) for n in t.split('-')]),1)
    timemarks3[where] = {'label': f'\n{t} years', 'style': {'transform':'rotate(45deg)', 'whitespace':'nowrap'} } # 'color': '#f50'


def make_net3(timepoint):
    '''Input: timepoint of interest. 
       Creates a network structure and the table with centrality indices.
    '''
    lab = timemarks3[timepoint]['label'][1:-6] # get label and remove '\n' in the beginning and ' years' at the end
    
    wm,ci,_,lay = read_res3(lab) # read in data 
    
    # tim estimates?
    wm_trim = wm.loc[abs(wm.weight)>0.01,].reset_index(drop=True)
    
    nodes = [{'data': {'id':node, 'label': '\n'.join(textwrap.wrap(get_label(node), width=20))}, 'classes':group, 
             'position':{'x':lay.loc[node,'x'], 'y':lay.loc[node,'y'] }} 
       for node,group in ci[['node','class']].itertuples(index=False) ]
    
    edges = [{'data': {'source':a, 'target':b, 'weight':w, 'width':round(abs(w)*20,2)}, 'classes':c} 
       for a,b,w,c in wm_trim.itertuples(index=False)]
    
    network = nodes+edges

    # Centrality indices tab
    ci_tab = ci.round(2).drop(columns='class')
    ci_tab.insert(1, 'Node', [get_label(name) for name in ci_tab['node']]) # Replace node name with its label 
    ci_tab = ci_tab.drop(columns='node')

    return network, ci_tab


stylenet3 = [ {'selector': 'node', 'style': {'label': 'data(label)', 'text-wrap':'wrap'} },
           # Edge opacty and width
           {'selector': 'edge', 'style': {'opacity': 'data(weight)', 'width': 'data(width)'}},
           # Color nodes by group
           {'selector': '.dep', 'style': {'background-color': 'lightblue'} },
           {'selector': '.cmr', 'style': {'background-color': 'pink'} },
           # Color edges by positive/negative weights
           {'selector': '.neg', 'style': {'line-color': 'red'} },
           {'selector': '.pos', 'style': {'line-color': 'blue'} } ]

