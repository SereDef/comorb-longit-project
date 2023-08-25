# Load dependency 
library(lavaan)
library(tidySEM)
library(parallel)

# Read in data
data <- readRDS('../Mats/raw_data.rds')

# ------------------------------------------------------------------------------
# -------------------------- Set-up and functions ------------------------------
# ------------------------------------------------------------------------------

# parallel::detectCores()

## Construct matrix determining all possible model structures -------------------
m = t(expand.grid(lapply(numeric(4), function(x) c(1, 0)))) # all possible combinations
m1 = matrix(1, 4, 12) # matrix of 4x12 ones to append later
# remove combinations not used and reorder
m0 = m[, colSums(m)!=1]; m0 = m0[,order(colSums(-m0))] # dim = 4x12
# Construct the full matrix and rename rows and columns 
mat = data.frame( cbind( rbind(m0,m1), rbind(m1[,-1], m0[,-1])), 
                  row.names = c('maCL_dep','maCL_cmr','maAR_dep','maAR_cmr',
                                'ltCL_dep','ltCL_cmr','ltAR_dep','ltAR_cmr'))
names(mat) <- c('full_st','no_maCL_dep','no_maCL_cmr','no_maAR_dep','no_maAR_cmr',
                'no_maCL','no_ma_dep','no_ma_CLcmr_ARdep','no_ma_CLdep_ARcmr','no_ma_cmr','no_maAR','no_ma',
                'no_ltCL_dep','no_ltCL_cmr','no_ltAR_dep','no_ltAR_cmr','no_ltCL','no_lt_dep',
                'no_lt_CLcmr_ARdep','no_lt_CLdep_ARcmr','no_lt_cmr','no_ltAR' ,'no_lt') 
rm(m0,m1)
# ------------------------------------------------------------------------------

## Select subsets of data (i.e., which variables and timepoints)
sel <- function(var, times=NULL) { 
  subs = names(data)[grep(paste(var, collapse='|'), names(data))]
  if (!is.null(times)) { subs = subs[grep(paste(paste0('_',times), collapse='|'), subs)] }
  return(subs)
}

## Create dataset and extract metadata (used in run_model)
make_df <- function(dep_name, dep_times, cmr_name, cmr_times, verbose=T) {
  
  dep = data[,sel(dep_name, dep_times)]; cmr = data[,sel(cmr_name, cmr_times)]
  
  if(ncol(dep)!=ncol(cmr)) { message('Error: unequal number of timepoints') }
  
  # Display characteristics of selection 
  temp_distance <- function(times) { tds <- c()
  for (i in 1:length(times)-1) { tds = c(tds, times[i+1]-times[i]) }
  return(tds)
  }
  
  df <- data.frame('id'= data$IDC)
  for (i in 1:ncol(dep)) { df[paste0('dep',i)] <- dep[,i]; df[paste0('cmr',i)] <- cmr[,i] }
  
  df = df[,-1]
  
  # select sample (at least one observation)
  samp = df[rowSums(is.na(df)) != ncol(df), ]
  
  if (verbose) {
    # Display temporal structure
    cat('Distance between constructs:', abs(dep_times-cmr_times), sep='\t')
    cat('\nDepression, temporal gap:', temp_distance(dep_times), sep='\t')
    cat('\nCMR marker, temporal gap:', temp_distance(cmr_times), sep='\t')
    
    # Display correlation
    message('\nCorrelations:')
    c = cor(df, use='pairwise.complete.obs', method='spearman')
    print(round(c,2))
    
    # Display sample size
    message('Sample size: ', nrow(samp))
  }
  
  output <- list("data" = samp, 
             "dep_time" = temp_distance(dep_times), 
             "cmr_time" = temp_distance(cmr_times) )
  
  return(output)
}

## Specify the formula for the model (used in run_sem)
formula <- function(dep_temp, cmr_temp, stationarity=T, mod='full_st') {
  # How many occasions 
  n_ocs = length(dep_temp)+1
  # Define model structure according to mat matrix
  for (r in row.names(mat)) {
    if (mat[mod][r,]==0) { assign(tolower(r),'0') } else { assign(gsub('lt','', tolower(r)), gsub('lt','', r)) }
  }
  # if (fix_lambda1) { l1='' } else { l1 = 'NA*'}
  # first loading is freely estimated and last one if fixed to 1 
  unit_effects = paste0('# Unit effects\neta_dep =~ NA*', paste0('dep', 1:(n_ocs-1), collapse=' + '), paste0(' + 1*dep',n_ocs), '\n',
                                        'eta_cmr =~ NA*', paste0('cmr', 1:(n_ocs-1), collapse=' + '), paste0(' + 1*cmr',n_ocs), '\n') 
  impulses = '# Impulses\n'
  for (i in 1:n_ocs) { impulses = paste0(impulses, 'u_dep',i,' =~ dep',i,'\ndep',i,' ~~ 0*dep',i,'\n') }
  for (i in 1:n_ocs) { impulses = paste0(impulses, 'u_cmr',i,' =~ cmr',i,'\ncmr',i,' ~~ 0*cmr',i,'\n') }
  
  regressions = '# Regressions\n'
  
  for (i in n_ocs:2) { 
    # assign different name to each parameter unless it was set to 0
    if (ar_dep!=0) { ar_depi = paste0(ar_dep,i-1) } else { ar_depi = ar_dep }
    if (cl_dep!=0) { cl_depi = paste0(cl_dep,i-1) } else { cl_depi = cl_dep }
    if (maar_dep!=0) { maar_depi = paste0(maar_dep,i-1) } else { maar_depi = maar_dep }
    if (macl_dep!=0) { macl_depi = paste0(macl_dep,i-1) } else { macl_depi = macl_dep }
    # Add depression regressions
    regressions = paste0(regressions, 
    'dep',i,' ~ ',ar_depi,'*dep',i-1,' + ',cl_depi,'*cmr',i-1,' + ',maar_depi,'*u_dep',i-1,' + ',macl_depi,'*u_cmr',i-1,'\n') }
  
  for (i in n_ocs:2) {
    # assign different name to each parameter unless it was set to 0
    if (ar_cmr!=0) { ar_cmri = paste0(ar_cmr,i-1) } else { ar_cmri = ar_cmr }
    if (cl_cmr!=0) { cl_cmri = paste0(cl_cmr,i-1) } else { cl_cmri = cl_cmr }
    if (maar_cmr!=0) { maar_cmri = paste0(maar_cmr,i-1) } else { maar_cmri = maar_cmr }
    if (macl_cmr!=0) { macl_cmri = paste0(macl_cmr,i-1) } else { macl_cmri = macl_cmr }
    # Add cardiometabolic regressions
    regressions = paste0(regressions, 
    'cmr',i,' ~ ',cl_cmri,'*dep',i-1,' + ',ar_cmri,'*cmr',i-1,' + ',macl_cmri,'*u_dep',i-1,' + ',maar_cmri ,'*u_cmr',i-1,'\n') }
  
  comevement = '# Comovement\n'
  for (i in 1:n_ocs) { comevement = paste0(comevement,'u_dep',i,' ~~ comv',i,'*u_cmr',i,'\n')}
  
  restrictions = '# Restrictions\n'
  for (i in 1:n_ocs) { 
    if (i < n_ocs) { ar = paste0(' + 0*u_dep', seq(i+1, n_ocs), collapse='') } else { ar = ''}
    restrictions = paste0(restrictions, 'u_dep',i,' ~~ 0*eta_dep + 0*eta_cmr',ar, ' + ',
                          paste0('0*u_cmr',c(1:n_ocs)[-i], collapse=' + '), '\n') }
  for (i in 1:n_ocs) {  
    if (i < n_ocs) { ar = paste0(' + 0*u_cmr', seq(i+1, n_ocs), collapse='') } else { ar = ''}
    restrictions = paste0(restrictions, 'u_cmr',i,' ~~ 0*eta_dep + 0*eta_cmr', ar, '\n') }
  
  if (stationarity) {
    constraints = '# Constraints\n'
    for (i in 2:(n_ocs-1)) {
      if (ar_dep!=0) { ar_dep_con = paste0('AR_dep',i,' == (AR_dep1)/',dep_temp[i]/dep_temp[1],'\n') } else { ar_dep_con=''}
      if (cl_dep!=0) { cl_dep_con = paste0('CL_dep',i,' == (CL_dep1)/',dep_temp[i]/dep_temp[1],'\n') } else { cl_dep_con=''}
      if (ar_cmr!=0) { ar_cmr_con = paste0('AR_cmr',i,' == (AR_cmr1)/',cmr_temp[i]/cmr_temp[1],'\n') } else { ar_cmr_con=''}
      if (cl_cmr!=0) { cl_cmr_con = paste0('CL_cmr',i,' == (CL_cmr1)/',cmr_temp[i]/cmr_temp[1],'\n') } else { cl_cmr_con=''}
      constraints = paste0(constraints, ar_dep_con, ar_cmr_con, cl_dep_con, cl_cmr_con) # maar_dep_con, macl_dep_con, maar_cmr_con, macl_cmr_con)
    }
  } else { constraints = '' }
  
  f = paste0(unit_effects, impulses, regressions, comevement, restrictions, constraints)
  return(f)
}

## Fit a single SEM model (used in run_model)
run_sem <- function(mod, stationarity, data, dep_temp, cmr_temp, bootn=1000) {
  m = lavaan::sem(formula(dep_temp, cmr_temp, stationarity, mod), data, se ='bootstrap', bootstrap = bootn)
  cat('.')
  return(m)
}

## Save model output as a graph in .png format (used in run_model)
save_semgraph <- function(model, n_ocs, name) {
  # TODO: cannot make it work with do.call...
  # nodes = c('eta_dep',rep('',n_ocs-1), paste0('dep',1:n_ocs), paste0('u_dep',1:n_ocs), 
  #           paste0('u_cmr',1:n_ocs), paste0('cmr',1:n_ocs), 'eta_cmr', rep('',n_ocs-1))
  # Specify the node structure
  if (n_ocs == 6) {
    lay = tidySEM::get_layout('eta_dep', '', '', '', '', '', 
                              'dep1', 'dep2', 'dep3', 'dep4', 'dep5', 'dep6', 
                              'u_dep1', 'u_dep2', 'u_dep3', 'u_dep4', 'u_dep5', 'u_dep6', 
                              'u_cmr1', 'u_cmr2', 'u_cmr3', 'u_cmr4', 'u_cmr5', 'u_cmr6', 
                              'cmr1', 'cmr2', 'cmr3', 'cmr4', 'cmr5', 'cmr6', 
                              'eta_cmr', '', '', '', '', '', rows=6)
  
    } else if (n_ocs == 5) {
    lay = tidySEM::get_layout('eta_dep', '', '', '', '',
                              'dep1', 'dep2', 'dep3', 'dep4', 'dep5',
                              'u_dep1', 'u_dep2', 'u_dep3', 'u_dep4', 'u_dep5',
                              'u_cmr1', 'u_cmr2', 'u_cmr3', 'u_cmr4', 'u_cmr5',
                              'cmr1', 'cmr2', 'cmr3', 'cmr4', 'cmr5', 
                              'eta_cmr', '', '', '', '', rows=6)
  
    } else if (n_ocs == 4) {
    lay = tidySEM::get_layout('eta_dep', '', '', '', 
                              'dep1', 'dep2', 'dep3', 'dep4', 
                              'u_dep1', 'u_dep2', 'u_dep3', 'u_dep4',
                              'u_cmr1', 'u_cmr2', 'u_cmr3', 'u_cmr4', 
                              'cmr1', 'cmr2', 'cmr3', 'cmr4',
                              'eta_cmr', '', '', '', rows=6)
  } else { cat('Too few occasions!\n') }
  
  # Remove edges estimated or set to 0.00 for readability
  e = get_edges(model); e$show[e$label=='0.00'] <- F
  # Create graph 
  p = tidySEM::graph_sem(model = model, layout = lay, edges=e)
  # Save image
  ggplot2::ggsave(paste0('../results/mod1_graphs/',name,'.png'), p, device = "png", width= 20, height = 15)
}

## Run the analyses! 
run_model <- function(dep_name, dep_times, cmr_name, cmr_times, verbose=F){
  # Make dataframe
  d <- make_df(dep_name, dep_times, cmr_name, cmr_times, verbose=verbose)
  if (verbose) { print(summary(d$data)) }
  
  # save variable names and summary for dashboard output 
  dep_summ <- do.call(cbind, lapply(data[,sel(dep_name, dep_times)], summary))
  cmr_summ <- do.call(cbind, lapply(data[,sel(cmr_name, cmr_times)], summary))
  dat_summ <- cbind(dep_summ, cmr_summ) 
  dat_summ <- rbind(dat_summ, 'N_obs'= nrow(data) - dat_summ["NA's",])
  
  cat('\n---------------------------------------\nFitting the model...')
  
  stat_seq <- c(T, rep(F, 2)) # ncol(mat)-1)) # Stationarity is true in the full models but false in all others
  
  system.time( # Run this bitch in parallel #parallel::mc
    fits <- mapply('run_sem', mod=names(mat)[1:3], stationarity=stat_seq, 
                 MoreArgs=list(data=d$data, dep_temp=d$dep_time, cmr_temp=d$cmr_time), 
                 SIMPLIFY=F) #, mc.cores=12)
  )
  cat('\nDone!\n---------------------------------------\n')
  
  # Initialize output dataframes
  fit_meas = data.frame(matrix(nrow = 46, ncol = 0)); estimates = data.frame(); failed = data.frame()
  
  for (f in names(fits)){
    
    if(lavInspect(fits[[f]], "converged")){
      # Inspect 
      # print(lavaan::summary(fits[[f]]))
      # Print out the layout 
      save_semgraph(fits[[f]], n_ocs=ncol(d$data)/2, name=paste(substr(dep_name,1,4),cmr_name,f, sep='_'))
      # Fit measures
      fm <- as.data.frame(lavaan::fitmeasures(fits[[f]])); names(fm) <- f
      fit_meas <- cbind(fit_meas, fm)
      # cat('\nFit measures:', fit_meas)
      # Parameters
      # stad_est = lavaan::standardizedSolution(mfit, type="std.all")[,'est.std'] # same but standardized
      es <- lavaan::parameterEstimates(fits[[f]])
      estimates <- rbind(estimates, cbind(rep(f, nrow(es)), es))
      
    } else { 
      message(f, ' - Model did not converge!\n') 
      problem <- t(data.frame(names(warnings()))); row.names(problem) <- rep(f, nrow(problem))
      failed <- rbind(failed, problem)               
      } 
  }
  save(dat_summ, fit_meas, estimates, failed, file = paste0('../results/',substr(dep_name,1,4),'_',cmr_name,'.RData'))
  
  # return(fits)
}
# ==============================================================================
# -------------------------------- Analysis ------------------------------------
# Self-reported depression =====================================================
# – BMI: 10 – 13 – 14 – 16 - 18 - 24 # (with mean BMI 15.5 – 18)?
run_model('sDEP_score', c(10.6, 12.8, 13.8, 16.6, 17.8, 23.8),
          'BMI',        c(10.7, 12.8, 13.8, 16,   17.8, 24.5)) # there is 17 too

# - total fat mass/ FMI : 10/11 – 12/13 – 14 – 15½ / 16½ - 18 – 24
run_model('sDEP_score', c(10.6, 12.8, 13.8, 16.6, 17.8, 23.8),
          'FMI',        c( 9.8, 11.8, 13.8, 15.4, 17.8, 24.5))

# Take it from here
run_model('sDEP_score', c(10.6, 12.8, 13.8, 16.6, 17.8, 23.8),
       'total_fatmass', c( 9.8, 11.8, 13.8, 15.4, 17.8, 24.5))

# – waist circumference: 10 – 13 – 16 – 14/25
run_model('sDEP_score', c(10.6, 12.8, 16.6, 23.8),
          'waist_circ', c(10.6, 12.8, 15.4, 24.5))


# Mother reported ==============================================================
# – BMI: 10 – 12 – 13 – 16  # (with mean BMI 15.5 – 18) 
run_model('mDEP_score', c(9.6, 11.7, 13.1, 16.7),
          'BMI',        c(9.8, 11.8, 12.8, 16)) # there is 17 too but lower correlation (0.2 vs. 0.4)

# – waist circumference: 10 – 12 – 13 – 16
run_model('mDEP_score', c(9.6, 11.7, 13.1, 16.7),
          'waist_circ', c(9.8, 11.8, 12.8, 15.4))

# – total fat mass/ FMI : 10 – 12 – 13.5 – 17.5 (or with mean fm 15.5 – 18) 
run_model('mDEP_score', c(9.6, 11.7, 13.1, 16.7),
          'FMI',        c(9.8, 11.8, 13.8, 17.8))

run_model('mDEP_score', c(9.6, 11.7, 13.1, 16.7),
       'total_fatmass', c(9.8, 11.8, 13.8, 17.8))

# – total lean mass/ LMI : 10 – 12 – 13.5 – 17.5 (or with mean lm 15.5 – 18)
run_model('mDEP_score', c(9.6, 11.7, 13.1, 16.7),
          'LMI',        c(9.8, 11.8, 13.8, 17.8))

run_model('mDEP_score', c(9.6, 11.7, 13.1, 16.7),
      'total_leanmass', c(9.8, 11.8, 13.8, 17.8))
# ==============================================================================

#mfit <- lavaan::sem(formula(d$dep_time, d$cmr_time), d$data)

# mfit2 <- lavaan::sem(formula(d$dep_time, d$cmr_time, 
#                              macl_dep='0',
#                              maar_dep='0' ), d$data)
# 
# , macl_dep=0, macl_cmr=0, maar_dep=0, maar_cmr=0
# ==============================================================================
# Inspect 
#lavaan::summary(mfit)

# round(lavaan::fitmeasures(mfit)[c('chisq','pvalue','cfi','rmsea','aic','bic')],3)
# 
# # Extractor functions
# lavaan::coef(mfit) # estimated parameters
# nest = lavaan::parameterEstimates(mfit) # data.frame containing all the model parameters 
# # similar but standardized parameter estimates
# sest = lavaan::standardizedSolution(mfit, type="std.all")

# f <- function(est){ 
#   print(t(data.frame(d$dep_time, d$cmr_time)))
#   cat('Unstandardized\n'); print(nest[grep(est, nest$label),])
#   cat('\nStandardized\n'); sest[grep(est, sest$label),] 
#   }
# f('comv')
# f('^AR_dep')
# f('^AR_cmr')
# f('maAR_dep') # dep short term AR
# f('maAR_cmr') # cmr short term AR
# f('^CL_dep') # cmr -> dep CL
# f('^CL_cmr') # dep -> cmr CL
# f('maCL_dep') # cmr -> dep short term CL
# f('maCL_cmr') # dep -> cmr short term CL
# 
# 
# # lavaan::fitted(mfit) # model-implied (fitted) covariance matrix
# lavaan::resid(mfit) # (unstandardized) residuals of a fitted model = difference between the observed and implied covariance matrix
# lavaan::lavResiduals(mfit) # gives more extensive information about the residuals (both raw and standardized + summary statistics).
# 
# lavaan::summary(mfit, standardized=T)
# 
# lavaan::inspect(mfit)
# 
# lavaan::modindices(mfit, sort = TRUE, maximum.number = 10)
# 
# lavaan::summary(mfit, fit.measures = TRUE)
# 
# # Self reported depression and BMI =============================================
# # Self depression score  - BMI: 11 – 13 – 14 – 17/15.5 – 18 – 24
# sel('BMI')
# 
# d <- make_df('sDEP_score', c(10.6, 12.8, 13.8, 16.6, 17.8), #, 23.8
#              'BMI',        c(10.7, 12.8, 13.8, 16.0, 17.8)) #, 24.5
# 
# mfit <- lavaan::sem(formula(5), d)
# 
# summary(mfit)
# # The model chi-square is highly significant, suggesting poor global model fit
# 
# varTable(mfit)
# # Parameter estimation can be hampered when the variances of variables in the model 
# # differ substantially (orders of magnitude).
# 
# # These look atrocious: CFI is < .95 (and much less than even .9), 
# # RMSEA is much greater than the .08 level that we would consider just ‘okay.’
# 
# round(inspect(mfit, what="cor.all")[1:12,1:12],2)
# # How does this compare to the observed correlations?
# round(lavCor(mfit),2)
# 
# round(resid(mfit, "cor")$cov,1)
# 
# modificationindices(mfit, minimum.value = 20)
# 
# sest = standardizedSolution(mfit, type="std.all")
# 
# sest[sest$label=='a',] # dep AR
# sest[sest$label=='g',] # cmr AR
# 
# sest[sest$label=='c',] # dep short term AR
# sest[sest$label=='i',] # cmr short term AR
# 
# sest[sest$label=='b',] # cmr -> dep CL
# sest[sest$label=='f',] # dep -> cmr CL
# 
# sest[sest$label=='d',] # cmr -> dep short term CL
# sest[sest$label=='h',] # dep -> cmr short term CL
# 
# semPlot::semPaths(mfit,what='std', nCharNodes=6, sizeMan=10,
#                   edge.label.cex=1.25, curvePivot = TRUE, fade=FALSE)
# 
# # Mother reported depression and fat mass 10-13 years
# m_totfat3 <- data.frame('dep1'= data[,'mDEP_score_9.6y'], 'fat1'= data[,'total_fatmass_9.8y'],
#                         'dep2'= data[,'mDEP_score_11.7y'],'fat2'= data[,'total_fatmass_11.8y'],
#                         'dep3'= data[,'mDEP_score_13.1y'],'fat3'= data[,'total_fatmass_12.8y'])# sDEP_score_12.8y
# # Mother reported depression and fat mass 10-17 years (including average fat between 15.4 and 17.8 measures)
# m_totfat4 <- cbind(m_totfat3, data.frame('dep4'= data[,'mDEP_score_16.7y']))
# m_totfat4$fat4 <- ifelse(rowSums(is.na(data[,sel('total_fatmass',c(15,17))]))<2, 
#                                  rowMeans(data[,sel('total_fatmass',c(15,17))], na.rm=T), NA)
# 
# # Self reported depression and fat mass 12-24 years
# m_totfat <- data[,c('sDEP_score_12.8y','total_fatmass_11.8y', # mDEP_score_11.7y
#                     'sDEP_score_17.8y','total_fatmass_17.8y',
#                     'sDEP_score_23.8y','total_fatmass_24.5y')]
# 
# names(m_totfat) <- c('dep1','fat1','dep2','fat2','dep3','fat3')
# 
# 
# samp = m_totfat[rowSums(is.na(m_totfat)) != ncol(m_totfat), ]
# 
# # mf <- cov(samp, use='pairwise.complete.obs')
# 
# 
# 
# mfit <- lavaan::sem(formula(4), m_totfat4)
# mfit
# summary(mfit)
# inspect(mfit)
# summary(mfit, fit.measures = TRUE)
# 
# mfit <- lavaan::sem(mod, m_totfat)
# 
# summary(mfit)
# 
# mfit <- lavaan::sem(mod, 
#                     sample.cov = mf, 
#                     sample.nobs = nrow(samp), 
#                     sample.mean = colMeans(samp,na.rm=T) )
# summary(mfit)
# 
# 
# library(lavaan)
# 
# # read in summary data:
# lower <- '
# 1.0
# .942 1.0
# .930  .948 1.0
# .914  .899  .892 1.0
# .922  .898  .888  .963 1.0
# .887  .887  .889  .887  .923 1.0
# .811  .778  .764  .680  .724  .710 1.0
# .828  .824  .803  .751  .756  .752  .95 1.0
# .824  .807  .791  .718  .742  .743  .939  .969 1.0
# .832  .807  .788  .732  .763  .753  .952  .963  .984 1.0
# .838  .809  .791  .744  .784  .782  .932  .963  .973  .978 1.0
# .831  .803  .788  .734  .775  .779  .914  .953  .959  .971  .988 1.0 '
# 
# corr <- lavaan::getCov(lower)
# 
# sds<-c(1.0904,1.0959,1.09499,1.07238,1.12339,1.08305,
#        .8700575,.91433,.927362,.92195445,.9011104,.897775)
# means<-c(5.261,5.414,5.394,5.425,5.424,5.426,
#          7.698,7.661,7.728,7.727,7.750,7.751 )
# 
# # produce the covariance matrix from the summary data:
# swb.cov <- lavaan::cor2cov(corr,sds,names=c("y1","y2","y3","y4","y5","y6",
#                                             "x1","x2","x3","x4","x5","x6"))
# 
# # specify the model:
# swb.mod <- '
# #unit effects
# eta_x =~ x6 + x5 + x4 + x3 + x2 + x1
# eta_y =~ y6 + y5 + y4 + y3 + y2 + y1
# 
# #impulses
# u_x1 =~ x1
# x1 ~~ 0*x1
# u_x2 =~ x2
# x2 ~~ 0*x2
# u_x3 =~ x3
# x3 ~~ 0*x3
# u_x4 =~ x4
# x4 ~~ 0*x4
# u_x5 =~ x5
# x5 ~~ 0*x5
# u_x6 =~ x6
# x6 ~~ 0*x6
# u_y1 =~ y1
# y1 ~~ 0*y1
# u_y2 =~ y2
# y2 ~~ 0*y2
# u_y3 =~ y3
# y3 ~~ 0*y3
# u_y4 =~ y4
# y4 ~~ 0*y4
# u_y5 =~ y5
# y5 ~~ 0*y5
# u_y6 =~ y6
# y6 ~~ 0*y6
# 
# #regressions
# x6 ~ a*x5 + b*y5 + c*u_x5 + d*u_y5 
# x5 ~ a*x4 + b*y4 + c*u_x4 + d*u_y4 
# x4 ~ a*x3 + b*y3 + c*u_x3 + d*u_y3 
# x3 ~ a*x2 + b*y2 + c*u_x2 + d*u_y2 
# x2 ~ a*x1 + b*y1 + c*u_x1 + d*u_y1 
# y6 ~ f*x5 + g*y5 + h*u_x5 + i*u_y5 
# y5 ~ f*x4 + g*y4 + h*u_x4 + i*u_y4 
# y4 ~ f*x3 + g*y3 + h*u_x3 + i*u_y3 
# y3 ~ f*x2 + g*y2 + h*u_x2 + i*u_y2 
# y2 ~ f*x1 + g*y1 + h*u_x1 + i*u_y1 
# 
# #co-movements
# u_x1 ~~ u_y1
# u_x2 ~~ u_y2
# u_x3 ~~ u_y3
# u_x4 ~~ u_y4
# u_x5 ~~ u_y5
# u_x6 ~~ u_y6
# 
# #restrictions
# u_x1 ~~ 0*eta_x + 0*eta_y + 0*u_x2 + 0*u_x3 + 0*u_x4 + 0*u_x5 + 0*u_x6 +
#         0*u_y2 + 0*u_y3 + 0*u_y4 + 0*u_y5 + 0*u_y6 
# u_x2 ~~ 0*eta_x + 0*eta_y +  0*u_x3 + 0*u_x4 + 0*u_x5 + 0*u_x6 +
#         0*u_y1 + 0*u_y3 + 0*u_y4 + 0*u_y5 + 0*u_y6
# u_x3 ~~ 0*eta_x + 0*eta_y +  0*u_x4 + 0*u_x5 + 0*u_x6 +
#         0*u_y1 + 0*u_y2 + 0*u_y4 + 0*u_y5 + 0*u_y6
# u_x4 ~~ 0*eta_x + 0*eta_y + 0*u_x5 + 0*u_x6 +
#         0*u_y1 + 0*u_y2 + 0*u_y3 + 0*u_y5 + 0*u_y6
# u_x5 ~~ 0*eta_x + 0*eta_y + 0*u_x6 +
#         0*u_y1 + 0*u_y2 + 0*u_y3 + 0*u_y4 + 0*u_y6
# u_x6 ~~ 0*eta_x + 0*eta_y + 0*u_y1 + 0*u_y2 + 0*u_y3 + 0*u_y4 + 0*u_y5
# u_y1 ~~ 0*eta_x + 0*eta_y + 0*u_y2 + 0*u_y3 + 0*u_y4 + 0*u_y5 + 0*u_y6 
# u_y2 ~~ 0*eta_x + 0*eta_y + 0*u_y3 + 0*u_y4 + 0*u_y5 + 0*u_y6
# u_y3 ~~ 0*eta_x + 0*eta_y + 0*u_y4 + 0*u_y5 + 0*u_y6
# u_y4 ~~ 0*eta_x + 0*eta_y + 0*u_y5 + 0*u_y6
# u_y5 ~~ 0*eta_x + 0*eta_y + 0*u_y6 
# u_y6 ~~ 0*eta_x + 0*eta_y '
# 
# #fit the model and show output:
# swb.fit <- lavaan::sem(swb.mod, sample.cov = swb.cov, sample.nobs = 135, 
#                        sample.mean = means)
# summary(swb.fit)
