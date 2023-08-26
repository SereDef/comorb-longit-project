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
rm(m0,m1); write.csv(mat, '../mats/model_structure.csv')
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
formula <- function(dep_temp, cmr_temp, mod='full_st') { # stationarity = T # just set it to true when full model is eastimated
  # How many occasions 
  n_ocs = length(dep_temp)+1
  # Define model structure according to mat matrix
  for (r in row.names(mat)) {
    if (mat[mod][r,]==0) { assign(gsub('lt','', tolower(r)),'0') } else { assign(gsub('lt','', tolower(r)), gsub('lt','', r)) }
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
  
  if (mod=='full_st') {
    constraints = '# Constraints\n'
    for (i in 2:(n_ocs-1)) {
      # if (ar_dep!=0) { ar_dep_con = paste0('AR_dep',i,' == (AR_dep1)/',dep_temp[i]/dep_temp[1],'\n') } else { ar_dep_con=''}
      # if (cl_dep!=0) { cl_dep_con = paste0('CL_dep',i,' == (CL_dep1)/',dep_temp[i]/dep_temp[1],'\n') } else { cl_dep_con=''}
      # if (ar_cmr!=0) { ar_cmr_con = paste0('AR_cmr',i,' == (AR_cmr1)/',cmr_temp[i]/cmr_temp[1],'\n') } else { ar_cmr_con=''}
      cl_cmr_con = paste0('CL_cmr',i,' == (CL_cmr1)/',cmr_temp[i]/cmr_temp[1],'\n')
      ar_dep_con = paste0('AR_dep',i,' == (AR_dep1)/',dep_temp[i]/dep_temp[1],'\n')
      cl_dep_con = paste0('CL_dep',i,' == (CL_dep1)/',dep_temp[i]/dep_temp[1],'\n')
      ar_cmr_con = paste0('AR_cmr',i,' == (AR_cmr1)/',cmr_temp[i]/cmr_temp[1],'\n')
      cl_cmr_con = paste0('CL_cmr',i,' == (CL_cmr1)/',cmr_temp[i]/cmr_temp[1],'\n')
      constraints = paste0(constraints, ar_dep_con, ar_cmr_con, cl_dep_con, cl_cmr_con) # maar_dep_con, macl_dep_con, maar_cmr_con, macl_cmr_con)
    }
  } else { constraints = '' }
  
  f = paste0(unit_effects, impulses, regressions, comevement, restrictions, constraints)
  return(f)
}

## Save model output as a graph in .png format (used in run_sem)
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

## Fit a single SEM model (used in run_model)
run_sem <- function(mod, data, dep_temp, cmr_temp, dep_name, cmr_name) {
  # Run model
  m <- lavaan::sem(formula(dep_temp, cmr_temp, mod), data, missing='fiml', se ='robust',
                   verbose=T) #, bootstrap = 1000)
  
  if (lavInspect(m, "converged")){
    # Print out the layout 
    save_semgraph(m, n_ocs=ncol(data)/2, name=paste(substr(dep_name,1,4),cmr_name,mod, sep='_'))
  } else {
    # Return warning message ( overwrites the model object )
    m <- paste(gsub('lavaan WARNING:', '', names(warnings())), collapse='\n') }
  
  return(m)
}

## Run the analyses! 
run_model <- function(dep_name, dep_times, cmr_name, cmr_times, which_models=names(mat), verbose=F){
  # Make dataframe
  d <- make_df(dep_name, dep_times, cmr_name, cmr_times, verbose=verbose)
  if (verbose) { print(summary(d$data)) }
  
  # save variable names and summary for dashboard output 
  dep_summ <- do.call(cbind, lapply(data[,sel(dep_name, dep_times)], summary))
  cmr_summ <- do.call(cbind, lapply(data[,sel(cmr_name, cmr_times)], summary))
  dat_summ <- cbind(dep_summ, cmr_summ) 
  dat_summ <- rbind(dat_summ, 'N_obs'= nrow(data) - dat_summ["NA's",])
  
  cat('\n---------------------------------------\nFitting the models...')
  
  # stat_seq <- c(T, rep(F, 2)) # ncol(mat)-1)) # Stationary is true in the full models but false in all others
  
  # Run this bitch (in parallel)
  start <- Sys.time(); cat(' started at: ', start)
  fits <- parallel::mcmapply('run_sem', mod=which_models,
                 MoreArgs=list(data=d$data, dep_temp=d$dep_time, cmr_temp=d$cmr_time, dep_name=dep_name, cmr_name=cmr_name), 
                 SIMPLIFY=F) # , mc.cores=12)
  end <- Sys.time()
  
  cat('\nDone! Runtime: ',round(end - start, 2),' hours\n---------------------------------------\n')
  
  # Initialize output dataframes
  fit_meas = data.frame(matrix(nrow = 55, ncol = 0)); estimates = data.frame(); failed = data.frame()
  
  for (f in names(fits)){
    if( is.character(fits[[f]]) ) {
      message(f) # , '\n\t', fits[[f]]) 
      problem <- data.frame(fits[[f]], row.names=f); names(problem) <- 'problem'
      failed <- rbind(failed, problem)           
    } else { 
      # Inspect # print(lavaan::summary(fits[[f]]))
      # Fit measures
      fm <- as.data.frame(lavaan::fitmeasures(fits[[f]])); names(fm) <- f
      fit_meas <- cbind(fit_meas, fm) # cat('\nFit measures:', fit_meas)
      # Parameters
      # stad_est = lavaan::standardizedSolution(mfit, type="std.all")[,'est.std'] # same but standardized
      es <- lavaan::parameterEstimates(fits[[f]])
      estimates <- rbind(estimates, cbind(rep(f, nrow(es)), es))    
    } 
  }
  save(dat_summ, fit_meas, estimates, failed, file = paste0('../results/mod1/',substr(dep_name,1,4),'_',cmr_name,'.RData'))
  
  # return(fits)
}

# ==============================================================================
# -------------------------------- Analysis ------------------------------------
# ==============================================================================

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
