# Load dependencies
invisible(lapply(c('psychonetrics','qgraph','dplyr','tidyr','psych','gdata'), require, 
                 character.only = TRUE));

data <- readRDS('../mats/raw_data.rds')

# which(sapply(data, is.character)) 

# Select variables from dataframe
sel <- function(var, times=NULL){
  subs = names(data)[grep(paste(var, collapse='|'), names(data))]
  if (!is.null(times)) { subs = subs[grep(paste(paste0('_',times), collapse='|'), subs)] }
  return(subs)
}

# ==============================================================================
# ========================= CROSS-LAG NETWORK MODEL ============================
# ==============================================================================

# Adjust odd timepoint BMI 10.7 -> 10.6 
data$BMI_10.6y <- data$BMI_10.7y

# Define waves of measurement (tot = 13 timepoints, ~1 year apart)
var_times = c(# '9.6y',
              '9.8y',
              '10.6y', 
              '11.8y',
              '12.8y',
              # '13.1y',
              '13.8y',
              '15.4y', # note: cmr only
              # '16y','17y',
              '16.6y',
              '17.8y',
              '18.7y', # note: dep only
              '21.9y', # note: dep only
              '22.9y', # note: dep only
              '23.8y', # note: dep only
              '24.5y') # note: cmr only

# Define variables of interest (27 nodes = 13 dep + 14 cmr)
var_names = c(paste0('sDEP', sprintf('%02d', 1:13)), # depression items
              'BMI', # collinear: 'height', 'weight'
              'waist_circ', # collinear: 'waist_hip_ratio'
              'FMI', # 'total_fatmass',
              'android_fatmass', # 'total_leanmass','LMI','trunk_fatmass','TFI',
              # 'liver_fat', # only one observation
              'DBP',
              'SBP',
              'PWV',
              'IMT', # note: measured twice
              # 'heart_rate', # unsure repeated measures are comparable
              # 'LVM', 'RWT', 'FS', # note: measured twice
              'HDL_chol',
              'LDL_chol', # 'tot_chol'
              'insulin',
              'triglyc',
              'glucose',
              'CRP' # 'IL_6' # only one timepoint
              # 'alcohol', 'canabis', 'smoking' # note: measured twice
              )

# Initialize design matrix (variable by times)
m <- matrix(ncol=length(var_times), nrow=length(var_names), dimnames=list(var_names, var_times))

# Fill matrix with variable names (when observed), else leave NA. 
for (v in row.names(m)){
  obs <- sel(v)
  for (t in colnames(m)) {
    time_present <- grepl(t, obs)
    if (any(time_present)) {  m[v,t] <- obs[which(time_present)] }
  }
}

# ----------------- pre-processing: min-max normalization -----------------------

normalize <- function(x, ...) {
  return( (x - min(x,...)) / (max(x,...) - min(x,...)) )
}

normdata <- sapply(data[,-c(1,2)], normalize, na.rm = TRUE)

# panelgvar() specifies a graphical vector-autoregression (GVAR) model on panel data
# When using only observed variables in the network, the panel data model takes the form of a random intercept cross-lagged panel model, 
# except the contemporaneous and between-subjects structures are modeled as networks. 
# This allows for the fixed-effects decomposition into temporal, contemporaneous, and between-subjects networks. 
# NOTE: when full-information maximum likelihood (FIML) estimation is used, the panel data model is a multi-level GVAR model 
# (as implemented in mlVAR) with only random intercepts (no random network parameters). 
# In the GVAR, temporal dependencies are modeled via a regression on the previous measurement occasion, 
# which leads to a matrix of regression coefficients that can also be used to draw a directed network model 
# often termed the temporal network because it encodes predictive effects over time. 
# The remaining variances and covariances (i.e., the covariance structure after controlling for the previous 
# measurement occasion) can be modeled as a GGM, which is also termed the contemporaneous network. 
# When time series of multiple subjects are available, a third GGM can be formed on the between-subject effects 
# (relationships between stable means)—also termed the between-subject network.

start_time <- Sys.time()
umod <- psychonetrics::panelgvar(normdata, vars = m, 
                                 estimator ='FIML', 
                                 verbose = TRUE) %>% 
  psychonetrics::runmodel() # run unpruned 
Sys.time() - start_time

save(umod, file='../results/mod2/unpruned_norm.RData')

ms = c('beta', 
  'omega_zeta_within', 'delta_zeta_within','sigma_epsilon_within',
  'omega_zeta_between','delta_zeta_between','sigma_epsilon_between')

pdf(paste0('CIplot_',format(Sys.Date(), '%y%m%d'),'.pdf'))
psychonetrics::CIplot(umod, ms) 
dev.off()

# ==============================================================================
load('../results/mod2/unpruned_norm.RData')
# ==============================================================================

save_info <- function(modobj = umod, modname='unpruned_norm_info', thresh=0.01, adjecency=FALSE) {
  # Get fit measures
  fit <- as.data.frame(modobj@fitmeasures)
  
  # Get networks 
  named_matrix <- function(type, model=modobj, rename=var_names){
    # type gets values = 'PDC', 'omega_zeta_within','omega_zeta_between'
    mat <- psychonetrics::getmatrix(model, type) # extract an estimated matrix from model
    rownames(mat) <- colnames(mat) <- rename
    return(mat)
  }
  
  # Get the three networks
  t_mat <- named_matrix('PDC')  # Temporal 
  c_mat <- named_matrix('omega_zeta_within')  # Contemporaneous
  b_mat <- named_matrix('omega_zeta_between') # Between
  
  layout <- qgraph::averageLayout(t_mat, c_mat, b_mat, layout = "spring")
  row.names(layout) <- var_names
  
  # Created plots and extract each matrix weight 
  plotnets <- function(network, title, filename, groups=c(rep('dep',13), rep('cmr',14))) {
    
    g <- qgraph::qgraph(network, title=title, threshold=thresh,
                        layout = layout, labels = var_names, 
                        theme = "colorblind", repulsion = 0.7, maximum = 1, 
                        groups = groups, color = c('#FFEAF5','#E1EDFF'),
                        vsize = 8, label.scale.equal=T, vsize2=7, shape ='ellipse',
                        border.width=0.5, label.cex = 1.1, legend = FALSE,
                        filetype='png', filename=file.path('../results/mod2/',filename))
    
    edges = as.data.frame(g[['Edgelist']][c('from','to','weight')])
    
    return(edges)
  }
  
  t_net = plotnets(t_mat, 'Temporal network', 't_plot')
  b_net = plotnets(b_mat, 'Between-person network', 'b_plot')
  c_net = plotnets(c_mat, 'Contemporaneous network', 'c_plot')
  
  descript <- function(matrix, type='crossNet') {
    if (type=='crossNet') { vector <- as.vector(gdata::lowerTriangle(matrix)) } else { vector <- as.vector(matrix) }
    estimates <- vector[vector !=0]
    descriptives <- round(psych::describe(estimates),3)
    return(t(descriptives))
  }  
  
  # Parameter descriptives
  desc <- data.frame(descript(t_mat, type='temporal'), descript(c_mat), descript(b_mat))
  names(desc) <- c('temp','cont','betw')
  
  print(t(desc))
  
  # Centrality indices
  t_cent <- qgraph::centrality_auto(t_mat)$node.centrality
  c_cent <- qgraph::centrality_auto(c_mat)$node.centrality
  b_cent <- qgraph::centrality_auto(b_mat)$node.centrality
  
  if (adjecency) {
    # pruned adjacency matrices
    t_adja <- 1*(psychonetrics::getmatrix(modobj, 'beta')!=0)
    b_adja <- 1*(psychonetrics::getmatrix(modobj, 'omega_zeta_between')!=0)
    c_adja <- 1*(psychonetrics::getmatrix(modobj, 'omega_zeta_within')!=0)
    # save to csv...?
  }
  
  # Variance covariance matrix for the input data
  var_cov <- cov(data[,sel(var_names, var_times)], use = 'pairwise.complete.obs')
  
  save(fit, layout, t_net, c_net, b_net, t_cent, c_cent, b_cent, var_cov,
      file = paste0('../results/mod2/',modname,'.RData'))
}

save_info()

# ==============================================================================
# The end :) 
# ==============================================================================