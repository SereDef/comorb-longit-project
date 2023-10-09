library(psychonetrics)
library(qgraph)
library(dplyr)
library(tidyr)
library(psych)
library(gdata)

data <- readRDS('../mats/raw_data.rds')

# which(sapply(data, is.character)) 
# Select variables from dataframe
sel <- function(var, times=NULL){
  subs = names(data)[grep(paste(var, collapse='|'), names(data))]
  if (!is.null(times)) { subs = subs[grep(paste(paste0('_',times), collapse='|'), subs)] }
  return(subs)
}
# sel(c(paste0('sDEP',0:2)), 23.8)

# ==============================================================================
# ========================= CROSS-LAG NETWORK MODEL ============================
# ==============================================================================

# Helper -----------------------------------------------------------------------
# Select subset of timepoints
# sub <- data[, sel(c('9.8y','10.6y','10.7y','11.8y', # T1
#                     '17.8y','15.4y',                # T2
#                     '23.8y','24.5y'))]              # T3
# # Quickly check the available timepoints and counts per variable 
# check <- function(var) {
#   message(var, '\nCOUNT:')
#   print(sapply(sub[,grep(var, names(sub))], function(x) sum(!is.na(x)) ))
#   message('SUMMARY')
#   print(summary(sub[,grep(var, names(sub))]))
# }

# Order alphabetically (i.e. per variable, per timepoint)
# sort(gsub('_9.8', '_09.8', colnames(sub)))
# check('waist_circ')

# Initiate dataframe -----------------------------------------------------------
df = data.frame('sex' = data$sex)

add <- function(vars, times, output=df) {
  # Convert times to string and add "y"
  times <- paste0(as.character(times), 'y')
  # Get vercor of names in oder 
  get_names <- unlist(lapply(vars, function(x) paste(x, times, sep='_')))
  # append to output dataframe 
  output[, get_names] <- data[,get_names]
  return(output)
}

# Add depression scores 
# df = add(paste0('sDEP', sprintf('%02d', 1:13) ), 
#          c(10.6, 17.8, 23.8))
# # Add cardio-metabolic variables 
# df = add(c('BMI','FMI','LMI'), 
#          c(11.8, 17.8, 24.5)) # excluding: 9.8, (10.7) and 15.4
# df = add(c('waist_circ'), 
#          c( 9.8, 15.4, 24.5)) # excluding 10.6 and 11.8
# df = add(c('DBP','SBP','PWV'), 
#          c(10.6, 17.8, 24.5)) # excluding: 9.8, 10.7 and 15.4
# df = add(c('HDL_chol','LDL_chol', 'insulin', 'triglyc','CRP'), # 'tot_chol'
#          c( 9.8, 17.8, 24.5)) # excluding: 15.4

# 10.6 -> 10.7
data$BMI_10.6y <- data$BMI_10.7y

var_times = c(# '9.6y',
          '9.8y',
          '10.6y', 
          '11.8y',
          '12.8y',
          # '13.1y',
          '13.8y',
          '15.4y', # cmr only
          # '16y','17y',
          '16.6y',
          '17.8y',
          '18.7y', # dep only
          '21.9y', # dep only
          '22.9y', # dep only
          '23.8y', # dep only
          '24.5y') # cmr only

var_names = c(paste0('sDEP', sprintf('%02d', 1:13)),
         # 'height', 'weight'
         'BMI',
         'waist_circ',
         # 'waist_hip_ratio', # ralated to waist_circ
         # 'total_fatmass', # index to height
         # 'total_leanmass', # index to height
         'FMI',
         # 'LMI',
         # 'trunk_fatmass','TFI',
         'android_fatmass',
         # 'liver_fat', # only one observation
         'DBP',
         'SBP',
         'PWV',
         'IMT', # 2 times
         # 'heart_rate', # unsure measure compares
         # 'LVM', # 2 times
         # 'RWT', # 2 times
         # 'FS', # 2 times
         # 'tot_chol',
         'HDL_chol',
         'LDL_chol',
         'insulin',
         'triglyc',
         'glucose',
         'CRP')
         # 'IL_6', # only one timepoint
         # 'alcohol', # 2 times
         # 'canabis', # 2 times
         # 'smoking') # 2 times

m <- matrix(ncol=length(var_times), nrow=length(var_names), dimnames=list(var_names, var_times))

for (v in row.names(m)){
  obs <- sel(v)
  for (t in colnames(m)) {
    time_present <- grepl(t, obs)
    if (any(time_present)) {
      m[v,t] <- obs[which(time_present)]
    }
  }
}

des = m # [c(paste0('sDEP0', 1:3),'BMI','android_fatmass','HDL_chol'),]

start_time <- Sys.time()
umod <- psychonetrics::panelgvar(data, vars = des, # lambda=
                                 # missing ='pairwise', 
                                 estimator ='FIML', 
                                 # optimizer = 'nlminb',
                                 verbose = TRUE) %>% 
  psychonetrics::runmodel() # run unpruned 
Sys.time() - start_time

save(umod, file='../results/mod2/unpruned.RData')

ms = c('beta', 
  'omega_zeta_within', 'delta_zeta_within','sigma_epsilon_within',
  'omega_zeta_between','delta_zeta_between','sigma_epsilon_between')

pdf('CIplot_2200923.pdf')
CIplot(umod, ms) 
dev.off()

# ==============================================================================

# removed: 'alcohol','canabis','smoking','android_fatmass','FS','glucose','heart_rate','IMT','LVM','RWT', # only at 18 and 24
# 'IL_6', # only at 10
# 'liver_fat', # only at 24 (earlier binary variable)
# 'waist_hip_ratio', # only 10/12 and 24
# 'CMR_age', 'DEP_age', 'DEP_score','TFI','height',
# 'fatmass','leanmass', 'weight', # using vars indexed to height

# ==============================================================================

names = c('felt\nmiserable\nor unhappy',
          'did not\nenjoy\nanything',
          'so tired\njust sat\naround',
          'was very\nrestless',
          'felt they\nwere no good\nanymore',
          'cried\na lot',
          'hard to\nthink or\nconcentrate',
          'hated\nthemselves',
          'felt they\nwere a\nbad person',
          'felt\nlonely',
          'nobody\nreally loved\nthem',
          'never\nas good as\nothers',
          'felt did\neverything\nwrong',
          'Body mass\nindex (BMI)',
          'Waist\ncircumference',
          'Fat mass\nindex (FMI)', # 'Lean mass\nindex (LMI)',
          'Android fat mass',
          'DBP',
          'SBP',
          'Pulse wave\nvelocity',
          'IMT',
          'HDL\ncholesterol',
          'LDL\ncholesterol',
          'Insulin',
          'Triglycerides',
          'Glucose',
          'C-reactive\nprotein')

load('/Users/Serena/Desktop/panel_network/results/mod2/unpruned_220923.RData')

# Get fit measures
ufit <- umod@fitmeasures

write.csv(ufit)

named_matrix <- function(model, type, rename=var_names){
  # type gets values = 'PDC' (namedTemporal), 'omega_zeta_within' (namedContemporaneous),
  # 'omega_zeta_between' (namedBetween)
  mat <- psychonetrics::getmatrix(model, type) # extract an estimated matrix from model
  rownames(mat) <- colnames(mat) <- rename
  return(mat)
}

# Get three networks
t_mat <- named_matrix(umod, 'PDC') # Temporal 
c_mat <- named_matrix(umod, 'omega_zeta_within') # Contemporaneous
b_mat <- named_matrix(umod, 'omega_zeta_between') # Between

layout <- qgraph::averageLayout(t_mat, c_mat, b_mat, layout = "spring")
# write.csv(layout, '../results/mod2/layout.csv')

plotnets <- function(network, title, filename, groups=c(rep('dep',13), rep('cmr',14))) {
  
  g <- qgraph::qgraph(network, title=title,
              labels=names, 
              threshold=0.01,
              theme = "colorblind", 
              repulsion = 0.7, 
              maximum = 1, 
              layout= 'spring',# layout, 
              groups = groups, 
              color = c('#FFEAF5','#E1EDFF'),
              vsize = 8, 
              label.scale.equal=T,
              vsize2=7,
              shape ='ellipse',
              border.width=0.5,
              label.cex = 1.1, legend = FALSE,
              filetype='png', 
              filename=file.path('/Users/Serena/Desktop/panel_network/results/mod2/',filename))
  
  return(g)
}

t = plotnets(t_mat, 'Temporal network', 't_plot')
b = plotnets(b_mat, 'Between-person network', 'b_plot')
c = plotnets(c_mat, 'Contemporaneous network', 'c_plot')

descript <- function(matrix, type='crossNet') {
  if (type=='crossNet') { vector <- as.vector(gdata::lowerTriangle(matrix)) } else { vector <- as.vector(matrix) }
  estimates <- vector[vector !=0]
  descriptives <- round(psych::describe(estimates),3)
  return(t(descriptives))
}  

# Parameter descriptives
desc <- data.frame(descript(t_mat, type='temporal'), descript(c_mat), descript(b_mat))
names(desc) <- c('temp','cont','betw')

desc = t(desc)

# Centrality indices df
t_cent <- qgraph::centrality_auto(t_mat)$node.centrality
c_cent <- qgraph::centrality_auto(c_mat)$node.centrality
b_cent <- qgraph::centrality_auto(b_mat)$node.centrality


# Pruned adjacency matrices
# t_adjacency <- 1*(psychonetrics::getmatrix(pmod, 'beta')!=0)
# b_adjacency <- 1*(psychonetrics::getmatrix(pmod, 'omega_zeta_between')!=0)
# c_adjacency <- 1*(psychonetrics::getmatrix(pmod, 'omega_zeta_within')!=0)


try = round(c_mat, 2)


try2 = 1*(c_mat!=1)
# ---------------------------------------------------------------------------
# Remove sex for now
# d = df[,-1]
# 
# nT = 3 # number of timepoints 
# nS = ncol(d) # number of symptoms/markers
# # Design matrix with nT columns and nS rows
# des <- matrix(as.vector(colnames(d)), nrow=nS/nT, byrow = TRUE) 

# Unpruned model, get fit
# panelgvar() specifies a graphical vector-autoregression (GVAR) model on panel data
# When using only observed variables in the network, the panel data model takes the form of a random intercept cross-lagged panel model, 
# except the contemporaneous and between-subjects structures are modeled as networks. 
# This allows for the fixed-effects decomposition into temporal, contemporaneous, and between-subjects networks. 
# NOTE: when full-information maximum likelihood (FIML) estimation is used, the panel data model is a multi-level GVAR model 
# (as implemented in mlVAR) with only random intercepts (no random network parameters). 
# start_time <- Sys.time()
# umod <- psychonetrics::panelgvar(d, vars = des, # lambda=
#                                  # missing ='pairwise', 
#                                  estimator ='FIML', 
#                                  # optimizer = 'nlminb',
#                                  verbose = TRUE) %>% 
#   psychonetrics::runmodel() # run unpruned 
# Sys.time() - start_time
# 
# save(umod, file='../results/mod2/unpruned.RData')
# Warning messages:
#   1: In doTryCatch(return(expr), name, parentenv, handler) :
#   restarting interrupted promise evaluation
# 2: In psychonetrics::runmodel(.) :
#   Information matrix or implied variance-covariance matrix was not positive semi-definite. This can happen because the model is not identified, or because the optimizer encountered problems. Try running the model with a different optimizer using setoptimizer(...).
# 3: In psychonetrics::runmodel(.) :
#   Model might not have converged properly: mean(abs(gradient)) > 1.
unprunedfit <- umod@fitmeasures


# pruned model, get fit, parameters, and matrices
pruned <- umod %>% psychonetrics::prune(recursive=T) # (recursively) remove parameters that are not significant and refit the model
prunedfit <- pruned@fitmeasures
prunedparameters <- pruned %>% psychonetrics::parameters()
# ==============================================================================
# In the GVAR, temporal dependencies are modeled via a regression on the previous measurement occasion, 
# which leads to a matrix of regression coefficients that can also be used to draw a directed network model 
# often termed the temporal network because it encodes predictive effects over time. 
# The remaining variances and covariances (i.e., the covariance structure after controlling for the previous 
# measurement occasion) can be modeled as a GGM, which is also termed the contemporaneous network. 
# When time series of multiple subjects are available, a third GGM can be formed on the between-subject effects 
# (relationships between stable means)—also termed the between-subject network.


# Get three networks
prunedtemporal <- namedMatrix(pruned, names, 'PDC')
prunedcontemporaneous <- namedMatrix(pruned, names, 'omega_zeta_within')
prunedbetween <- namedMatrix(pruned, names, 'omega_zeta_between')

layout <- qgraph::averageLayout(prunedtemporal, prunedcontemporaneous, prunedbetween, 
                                layout = "spring")
write.csv(layout, '../results/mod2/layout.csv')

groups <- c(rep('dep',13), rep('cmr',12))

plotnets <- function(network, title, filename) {
  
  g <- qgraph(network, title=title,
              labels=names, 
              threshold=0.01,
             theme = "colorblind", 
             repulsion = 0.7, 
             maximum = 1, 
             layout= 'spring',# layout, 
             groups = groups, 
             color = c('#FFEAF5','#E1EDFF'),
             vsize = 8, 
             label.scale.equal=T,
             vsize2=7,
             shape ='ellipse',
             border.width=0.5,
             label.cex = 1.1, legend = FALSE,
             filetype='png', 
             filename=file.path('./assets',filename))
  
  return(g)
}

t = plotnets(prunedtemporal, 'Temporal network', 'tempPlot')
b = plotnets(prunedbetween, 'Between-person network', 'betweenPlot')
c = plotnets(prunedcontemporaneous, 'Contemporaneous network', 'contempPlot')

descript <- function(matrix, type='crossNet') {
  if (type=='crossNet') { vector <- as.vector(gdata::lowerTriangle(matrix)) } else { vector <- as.vector(matrix) }
  estimates <- vector[vector !=0]
  descriptives <- psych::describe(estimates)
  return(descriptives)
}  

prunedtemp.desc <- descript(prunedtemporal, type='temporal')
prunedcont.desc <- descript(prunedcontemporaneous)
prunedbet.desc  <- descript(prunedbetween)

# pruned centrality
PDCcentp <- qgraph::centrality_auto(prunedtemporal)
PDCcentP <- PDCcentp$node.centrality
contCentp <- qgraph::centrality_auto(psychonetrics::getmatrix(pruned, 'omega_zeta_within'))
contCentP <- contCentp$node.centrality
betCentp <- qgraph::centrality_auto(psychonetrics::getmatrix(pruned, 'omega_zeta_between'))
betCentP <- betCentp$node.centrality

# pruned adjacency matrices
prune.adjacencyT <- 1*(psychonetrics::getmatrix(pruned, 'beta')!=0)
prune.adjacencyB <- 1*(psychonetrics::getmatrix(pruned, 'omega_zeta_between')!=0)
prune.adjacencyC <- 1*(psychonetrics::getmatrix(pruned, 'omega_zeta_within')!=0)

# pruned invariance
# genderComp <- df %>% tidyr::drop_na(sex)

# Pconfigural <-  psychonetrics::panelgvar(genderComp, vars = des, 
#                                          estimator = 'FIML', missing = 'pairwise', 
#                                          groups= 'sex',
#                                          beta = prune.adjacencyT, 
#                                          omega_zeta_within = prune.adjacencyC, 
#                                          omega_zeta_between = prune.adjacencyB) %>% 
#   psychonetrics::runmodel()
# prunedconfig.fit <- Pconfigural@fitmeasures
# 
# Pconstrained <- Pconfigural %>% psychonetrics::groupequal('beta') %>% 
#   psychonetrics::groupequal('omega_zeta_within') %>% 
#   psychonetrics::groupequal("omega_zeta_between") %>% 
#   psychonetrics::runmodel()
# prunedconstrained.fit <- Pconstrained@fitmeasures
# 
# prunedcomparison <- psychonetrics::compare(Pconfigural = Pconfigural, 
#                                            Pconstrained = Pconstrained)
# 
# aicDif <- Pconfigural@fitmeasures$aic.ll-Pconstrained@fitmeasures$aic.ll
# bicDif <- Pconfigural@fitmeasures$bic-Pconstrained@fitmeasures$bic
# prunedICs <- cbind(aicDif, bicDif)
# 
# #stepup model, get fit, parameters and matrices
# stepup <- pruned %>% psychonetrics::stepup()
# stepupfit <- stepup@fitmeasures
# stepupparameters <- stepup %>% parameters()
# 
# stepuptemporal <- namedTemporal(stepup, names)
# stepupcontemporaneous <- namedContemporaneous(stepup, names)
# stepupbetween <- namedBetween(stepup, names)
# stepuptemp.desc <- temporalDescript(stepuptemporal )
# stepupcont.desc <- crossNetDescript(stepupcontemporaneous )
# stepupbet.desc <- crossNetDescript(stepupbetween )
# 
# #stepup centrality
# PDCcents<- centrality_auto(stepuptemporal )
# PDCcentS <- PDCcents$node.centrality
# contCents<- centrality_auto(getmatrix(stepup, "omega_zeta_within"))
# contCentS <- contCents$node.centrality
# betCents<- centrality_auto(getmatrix(stepup, "omega_zeta_between"))
# betCentS <- betCents$node.centrality
# 
# #stepup adjacency matrices
# stepupadjacencyT <- 1*(getmatrix(stepup, "beta")!=0)
# stepupadjacencyB<- 1*(getmatrix(stepup, "omega_zeta_between")!=0)
# stepupadjacencyC<- 1*(getmatrix(stepup, "omega_zeta_within")!=0)
# 
# # pruned invariance
# Sconfigural<- panelgvar(genderComp, vars= des, estimator= "FIML", missing= "pairwise", groups= "Gender",
#                         beta= stepupadjacencyT, omega_zeta_within= stepupadjacencyC, omega_zeta_between= stepupadjacencyB)%>% runmodel
# stepupconfig.fit <- Sconfigural@fitmeasures
# 
# Sconstrained<- Sconfigural%>% groupequal("beta")%>% groupequal("omega_zeta_within")%>% groupequal("omega_zeta_between")%>% runmodel()
# stepupconstrained.fit <- Sconstrained@fitmeasures
# 
# stepupcomparison<- psychonetrics::compare(Sconfigural= Sconfigural, Sconstrained= Sconstrained)
# 
# SaicDif<- Sconfigural@fitmeasures$aic.ll-Sconstrained@fitmeasures$aic.ll
# SbicDif<- Sconfigural@fitmeasures$bic-Sconstrained@fitmeasures$bic
# stepupICs <- cbind(SaicDif, SbicDif)


# ==============================================================================
# =========== Cross-sectional network analysis (single timepoint) ==============
# ==============================================================================

library(qgraph)

run_net <- function(times, mean_times='',
                    remove = c('score','age','weight','total_fatmass','total_leanmass',
                               'trunk_fatmass','TMI','TFI','height'), 
                    verbose = TRUE) { 
  # Subset columns 
  d <- data[,sel(times)]
  d = d[, -which( names(d) %in% sel(remove) )] # e.g., remove total depression score and age at measurement
  
  # Select cases with at least one observation 
  d_obs = d[ apply(d, 1, function(x) { sum(is.na(x))<ncol(d) } ),]; n_obs <- nrow(d_obs)
  
  if (verbose) { 
    n_depitems <- sum(names(d) %in% sel('DEP'))
    message(n_depitems, ' depression items and ', ncol(d)-n_depitems, ' CMR markers.')
    cat(n_obs, 'observations.\n') 
    cat(names(d), sep='\n') }
  
  # Correlation matrix
  mat <- cor(d_obs, use = 'pairwise.complete.obs', method = 'spearman')
  
  message('Fitting network...')
  
  # Fit the network
  gra <- qgraph::qgraph(mat, 
                        labels=names(d_obs),
                        graph = 'glasso', # with EPIC model selection
                        sampleSize = n_obs,
                        # estimator ='FIML', 
                        layout = 'spring', # node placement
                        tuning = 0, # EBIC hyperparameter (gamma): [0-0.5] higher=fewer connections
                        lambda.min.ratio = 0.01 # minimal lambda ratio used in EBICglasso, defaults to 0.01.
  ) # For more control on the estimation procedure, use EBICglasso() [= qgraph(..., graph="glasso")]
  
  # Weight matrix
  wm <- getWmat(gra) # estimated weights matrix (can be obtained directly using EBICglasso)
  wm <- as.data.frame(matrix(wm, dimnames=list(t(outer(colnames(wm), rownames(wm), FUN=paste)), NULL)))
  
  # Centrality indices
  ci <- centralityTable(gra)[,3:5]
  ci <- reshape(ci, idvar = 'node', timevar = 'measure', direction = 'wide')
  names(ci) = gsub('value.', '', names(ci))
  
  # Fit measures 
  fit <- as.data.frame(ggmFit(gra, covMat=mat, sampleSize=n_obs)$fitMeasures)
  fit['N_obs'] <- n_obs
  
  # Layout
  layout <- gra$layout
  row.names(layout) <- gra$graphAttributes$Nodes$names
  
  # Save output
  if (mean_times=='') { temp = paste(gsub('y','',times), collapse='-') } else { temp = mean_times } 
  save(wm, ci, fit, layout, file = paste0('../results/mod3/crosnet_',temp,'y.RData'))
  
  return(gra)
}

cn9.7 <- run_net(times = c('9.6y','9.8y'))
cn10.65 <- run_net(times = c('10.6y','10.7y'))
cn11.75 <- run_net(times = c('11.7y','11.8y'))

cn12.8 <- run_net(times = c('12.8y'))
cn13.1 <- run_net(times = c('13.1y'))
cn13.8 <- run_net(times = c('13.8y'))
cn15.4 <- run_net(times = c('15.4y')) # note: CMR only 

cn16.6 <- run_net(times = c('16.6y','17y'))
cn16.7 <- run_net(times = c('16.7y','17y'))

cn17.8 <- run_net(times = c('17.8y')) # note: also lifestyle 

cn18.7 <- run_net(times = c('18.7y')) # note: depression only 
cn21.9 <- run_net(times = c('21.9y')) # note: depression only 
cn22.9 <- run_net(times = c('22.9y')) # note: depression only 

cn24.1 <- run_net(times = c('23.8y','24.5y'))

# averageLayout(cn24.1, cn22.9)

# cor_auto detects ordinal variables (with up to 7 unique integer values) and uses
# lavaan to estimate polychoric, polyserial and Pearson correlations.
# PROBLEM: tot_chol_9.8y, LDL_chol_9.8y
# mat <- qgraph::cor_auto(d)

# cp <- centralityPlot(gra) # centrality indices

# ==============================================================================
# library(mgm) # Estimating Time-Varying Mixed Graphical Models in High-Dimensional Data
# 
# d <- data[,sel(c('9.6y','9.8y'))]
# d = d[, -which(names(d)%in%sel(c('score','age')))]
# 
# dc <- d[complete.cases(d), ] # doesnt take missign values 
# 
# vtype <- c(rep('c',13), rep('g',18))
# vlevl <- c(rep( 3, 13), rep( 1 ,18))
# 
# set.seed(1)
# fit_mgm <- mgm::mgm(data = dc,
#                    type = vtype,
#                    levels = vlevl,
#                    k = 2,
#                    lambdaSel = "CV")
# 
# qgraph::qgraph(fit_mgm$pairwise$wadj,
#        edge.color = fit_mgm$pairwise$edgecolor,
#        layout = 'spring',
#        nodeNames = names(d),
#        legend = TRUE)
# fit_mgm$pairwise # = weighted adjacency matrix and the signs (if defined) of the parameters
# fit_mgm$interactions # list of all recovered interactions (cliques) and a list of parameters associated with all cliques; 
# fit_mgm$intercepts # all estimated thresholds/intercepts
# fit_mgm$nodemodels # list with the p glmnet objects from which all above results are computed. 
# We inspect the weigthed adjacency matrix stored in fit_mgm$pairwise$wadj

# ==============================================================================
# ==============================================================================
