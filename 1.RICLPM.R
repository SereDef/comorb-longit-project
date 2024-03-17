# ==============================================================================
# =========== 1. RANDOM-INTERCEPT CROSS-LAG PANEL MODEL (RI-CLPM) ==============
# ==============================================================================

args = commandArgs(trailingOnly = TRUE) # provide access to a copy of the command line arguments supplied

if (length(args) == 0) {
  stop("Supply output folder name!")
} else {
  out_folder <- args[1] # e.g. 'mod1_ri_pstat' # parameters stationary
}

# Set parameters
dir.create(out_folder)

stationarity <- !(grepl('pfree', out_folder)) # model regression coefs as stationary over time 

# The type and scale of time constraints
time_constraints=c('*','years')
# Speak or quite
verbose=FALSE

# Load dependencies
invisible(lapply(c('lavaan','foreach'), require, character.only = TRUE));
# Note: I also tried parallel and pbapply for parallel processing but foreach worked best

# Read in data
data <- readRDS('./raw_data.rds')

# ==============================================================================
# -------------------------- Set-up and functions ------------------------------
# ==============================================================================

## Select subsets of data (i.e., which variables and timepoints)
sel <- function(var, times=NULL) { 
  subs = names(data)[grep(paste(var, collapse='|'), names(data))]
  if (!is.null(times)) { subs = subs[grep(paste(paste0('_',times), collapse='|'), subs)] }
  return(subs)
}

## Create dataset and extract metadata (used in run_model)
make_df <- function(name1, times1, name2, times2, 
                    group=NULL, # stratify sample?
                    var1='dep',
                    var2='cmr', 
                    return_time_scale='years',
                    transform=TRUE,
                    normalize=TRUE, 
                    verbose=TRUE) {
  
  # Select the two sets of measurement
  d1 = data[,sel(name1, times1)]
  d2 = data[,sel(name2, times2)]
  
  # Check input 
  if(ncol(d1)!=ncol(d2)) { message('Error: unequal number of timepoints') }
  
  # save variable names and summary for dashboard output 
  dat_summ <- cbind(do.call(cbind, lapply(d1, summary)),
                    do.call(cbind, lapply(d2, summary)))
  dat_summ <- rbind(dat_summ, 'n_obs'= nrow(data) - dat_summ["NA's",])
  
  # Rename variables
  names(d1) <- paste0(var1,1:ncol(d1))
  names(d2) <- paste0(var2,1:ncol(d2))
  
  df <- as.data.frame(cbind(d1,d2))
  
  # Add strata variable(s)
  if (!is.null(group)) df[,group] <- data[,group]
  
  # Subset sample (at least one observation per variable)
  samp = df[(rowSums(!is.na(d1)) > 0) & (rowSums(!is.na(d2)) > 0), ]
  message('Removing ', nrow(df)-nrow(samp),' rows with insufficient data.')
  
  # Data summary after subsetting
  dat_summ_samp <- do.call(cbind, lapply(samp, summary))
  dat_summ_samp <- rbind(dat_summ_samp, 'n_obs'= nrow(samp) - dat_summ_samp["NA's",])
  
  dat_summ <- round(cbind(dat_summ, dat_summ_samp),2)
  
  # Correlations
  corrs = round(cor(samp[,c(names(d1), names(d2))], use='pairwise.complete.obs', method='pearson'), 2)
  
  if (transform) {
    # Square root for depression and log for CMR
    samp[,names(d1)] <- sapply(samp[,names(d1)], sqrt) 
    samp[,names(d2)] <- sapply(samp[,names(d2)], log1p) # note: avoid -Inf values in case of 0s
  }
  
  # Normalize data 
  if (normalize) { 
    minmax_norm <- function(x, ...) { return( (x - min(x,...)) / (max(x,...) - min(x,...)) )}
    samp[,c(names(d1), names(d2))] <- sapply(samp[,c(names(d1), names(d2))], minmax_norm, na.rm = TRUE) 
  }
  
  if (verbose) {
    # Display characteristics of selection 
    temp_distance <- function(times) { 
      tds <- c()
      for (i in 1:length(times)-1) tds = c(tds, times[i+1]-times[i])
      return(tds)
    }
    # Display temporal structure
    cat('Distance between constructs:', abs(times1-times2), sep='\t')
    cat('\nDepression, temporal gap:', temp_distance(times1), sep='\t')
    cat('\nCMR marker, temporal gap:', temp_distance(times2), sep='\t')
    
    # Display correlation
    message('\nCorrelations:')
    print(corrs)
    
    # Display sample size
    message('\nSample size: ', nrow(samp))
  }
  
  # Finally, add age at measurement
  if (return_time_scale=='years') { 
    scale <- 1; fdigit <- 1
  } else if (return_time_scale=='months') { 
    scale <- 12; fdigit <- 0
  }
  age_vars1 <- data[,sel(paste0(toupper(var1),'_age'), times1)] * scale
  age_vars2 <- data[,sel(paste0(toupper(var2),'_age'), times2)] * scale 
  
  med_ages1 <- unname(round(sapply(age_vars1, median, na.rm=TRUE), fdigit))
  med_ages2 <- unname(round(sapply(age_vars2, median, na.rm=TRUE), fdigit))
  
  return(list('sample'=samp, 'corrs'=corrs, 'summary'=dat_summ,
              'med_ages1'=med_ages1, 'med_ages2'=med_ages2))
}

## Specify the formula for the model (used in run_sem)
riclpm_formula <- function(var1='dep', var2='cmr', 
                           n_ocs=NULL, meas_time=list(c(), c()), 
                           rel='*', # specify temporal equality using power ('^') or linear ('*') relations 
                           strata=c(),
                           stationarity=TRUE,
                           verbose=TRUE) {
  
  # Check input ----------------------------------------------------------------
  if (is.null(n_ocs) & (length(meas_time[[1]])==0 | length(meas_time[[1]])==0)) { 
    stop('Provide total number of occasions or time points of measurement!') }
  
  if (length(meas_time[[1]]) != length(meas_time[[2]])) {
    stop('You need the same number of time points for each of the two variables.') }
  
  if (!is.null(n_ocs) & !length(meas_time[[1]]) %in% c(0, n_ocs)) { # contraddictory info 
    stop('Number of occations and measurement times provided do not agree.') }
  
  # ----------------------------------------------------------------------------
  # How many occasions 
  if (length(meas_time[[1]]) > 0) { 
    temp_var1 <- meas_time[[1]]
    temp_var2 <- meas_time[[2]] 
  } else { 
    temp_var1 <- temp_var2 = 1:n_ocs 
  }
  
  if (is.null(n_ocs)) n_ocs <- length(temp_var1)
  
  # Create between components (random intercepts) ------------------------------
  ri_l = '1*' # Set random intercept loadings to 1
  random_intercepts <- paste0('# Random intercepts\n',
     'ri_',var1,' =~ ', paste0(ri_l,var1,1:n_ocs, collapse=' + '), '\n',
     'ri_',var2,' =~ ', paste0(ri_l,var2,1:n_ocs, collapse=' + '), '\n')
  
  # Create within-components (or within-person centered variables) -------------
  # Note: factor loadings are set to 1 in non-stationary model and freely estimated in stationary models
  if (stationarity) { l = 'NA*' } else { l = '1*' }
  impulses <- paste0('# Impulses\n',
     paste0('w_',var1,1:n_ocs,' =~ ',l,var1,1:n_ocs, collapse='\n'), '\n',
     paste0('w_',var2,1:n_ocs,' =~ ',l,var2,1:n_ocs, collapse='\n'), '\n')
  
  # Define estimate names, in main and grouped analyses 
  ename <- function(etype, var, n, strata){
    
    if (var != '') var <- paste0('_',var)
    
    if (length(strata) < 1) {
      name <- paste0(etype, var, n)
    } else {
      name <- paste0('c(', paste0(etype, var, n, '_',strata, collapse=', '), ')')
    }
    return(name)
  }
  
  # Estimate lagged effects between within-person centered variables -----------
  regressions <- '# Regressions\n'
  
  for (i in 1:(n_ocs-1)) {    # iterator

    ar_1 <- ename('AR', var1, i, strata=strata)
    ar_2 <- ename('AR', var2, i, strata=strata)
    cl_1 <- ename('CL', var1, i, strata=strata)
    cl_2 <- ename('CL', var2, i, strata=strata)
      
    regressions <- paste0(regressions, 
       'w_',var1,i+1,' ~ ',ar_1,' * w_',var1,i,' + ',cl_1,' * w_',var2,i,'\n',
       'w_',var2,i+1,' ~ ',ar_2,' * w_',var2,i,' + ',cl_2,' * w_',var1,i,'\n') 
  }
  
  # Estimate covariance between within-person centered variables and RIs
  rcovs <- sapply(2:n_ocs, function(i) ename('rcov','',i,strata=strata))
  covariances <- paste0('# Covariances\n',
     # Estimate covariance (or correlation) between within-components at first wave
     'w_',var1,'1 ~~ ', ename('cor','',1,strata=strata), ' * w_',var2,'1\n',
     # Estimate residual covariances (between residuals of within-person vars, i.e., innovations)
     paste0('w_',var1,2:n_ocs,' ~~ ',rcovs,' * w_',var2,2:n_ocs, collapse='\n'), '\n', 
     # Estimate covariance of the random intercepts
     'ri_',var1,' ~~ ', ename('covRI','','',strata=strata), ' * ri_',var2,'\n')
  
  # Set variances of within-components at first wave to 1 in standardized model
  if (stationarity) { w1 = '1*' } else { w1 = '' }
  
  # Estimate (residual) variances of within-person centered variables and RIs
  rvars1 <- sapply(2:n_ocs, function(i) ename('rvar', var1, i, strata=strata))
  rvars2 <- sapply(2:n_ocs, function(i) ename('rvar', var2, i, strata=strata))
  variances = paste0('# Variances\n',
     # Estimate variances of within-components at first wave (or set it to 1)
     paste0('w_',var1,'1 ~~ ',w1,'w_',var1,'1\n',
            'w_',var2,'1 ~~ ',w1,'w_',var2,'1\n'),
     # Estimate (and label) residual variances of within-person centered variables
     paste0('w_',var1,2:n_ocs,' ~~ ',rvars1,'*w_',var1,2:n_ocs, collapse='\n'),'\n', 
     paste0('w_',var2,2:n_ocs,' ~~ ',rvars2,'*w_',var2,2:n_ocs, collapse='\n'),'\n', 
     # Estimate variance of random intercepts
     'ri_',var1,' ~~ ri_',var1,'\nri_',var2,' ~~ ri_',var2, '\n')
  
  # Finally impose constraints if any
  constraints = ar_var1_con = ar_var2_con = cl_var1_con = cl_var2_con = cor_con = rvar_con = ''
  
  if (stationarity) { # Do not constrain if the models are stratified
    constraints = '\n# Constraints\n'
      
    # paste0('# Constrain grand means over time\n',
    # paste0(var1,1:n_ocs, collapse=' + '), ' ~ mean_',var1,'*1\n',
    # paste0(var2,1:n_ocs, collapse=' + '), ' ~ mean_',var2,'*1\n')
    
    # Weight for unequal temporal gaps between measurements
    i = 1:(n_ocs-2) # new iterator
    # set precision of temporal difference measure (NOTE: floats don't work with power relations)
    if (rel=='*') { prec=1 } else if (rel=='^') { prec=0 } 
    
    ar_var1_con <- paste0('AR_',var1,i,  rel,round(temp_var1[i+1]-temp_var1[i], prec),
                      ' == AR_',var1,i+1,rel,round(temp_var1[i+2]-temp_var1[i+1],prec), collapse='\n')
    ar_var2_con <- paste0('AR_',var2,i,  rel,round(temp_var2[i+1]-temp_var2[i], prec), 
                      ' == AR_',var2,i+1,rel,round(temp_var2[i+2]-temp_var2[i+1], prec), collapse='\n')
    cl_var1_con <- paste0('CL_',var1,i,  rel,round(temp_var1[i+1]-temp_var2[i], prec), 
                      ' == CL_',var1,i+1,rel,round(temp_var1[i+2]-temp_var2[i+1], prec), collapse='\n')
    cl_var2_con <- paste0('CL_',var2,i,  rel,round(temp_var2[i+1]-temp_var1[i], prec), 
                      ' == CL_',var2,i+1,rel,round(temp_var2[i+2]-temp_var1[i+1], prec), collapse='\n')
    
    # Allow for standardized estimates, even when constrained 
    i = 1:(n_ocs-1) # iterator
      
    # Compute correlations between the within-components themselves at each wave 
    cor_con <- paste0('cor', 2:n_ocs,
       ' := AR_',var1,i,'*CL_',var2,i,' + ',
           'CL_',var1,i,'*AR_',var2,i,' + ',
           'AR_',var1,i,'*AR_',var2,i,'*cor',i,' + ',
           'CL_',var1,i,'*CL_',var2,i,'*cor',i,' + rcov',2:n_ocs, collapse='\n')
    
    # Constrain residual variances of within-components such that the total variance of each 
    # within-component equals 1
    rvar_con <- paste0(
      'rvar_',var1,2:n_ocs,
      ' == 1 - (AR_',var1,i,'*AR_',var1,i,' + ',
               'CL_',var1,i,'*CL_',var1,i,' + 2*AR_',var1,i,'*CL_',var1,i,'*cor',i,')\n',
      'rvar_',var2,2:n_ocs,
      ' == 1 - (AR_',var2,i,'*AR_',var2,i,' + ',
               'CL_',var2,i,'*CL_',var2,i,' + 2*AR_',var2,i,'*CL_',var2,i,'*cor',i,')', collapse = '\n')
  }
  
  constraints = paste0(constraints, 
                       ar_var1_con,'\n', ar_var2_con,'\n', cl_var1_con,'\n', cl_var2_con,'\n',
                       cor_con,'\n', rvar_con)
  
  # Paste everything together 
  f = paste0(random_intercepts, impulses, regressions, covariances, variances, 
             constraints)
  
  if (verbose) { message('MODEL SYNTAX:'); cat(f) }
  
  return(f)
}

## Fit a single SEM model (used in run_all_models)
run_lavaan <- function(param, group=NULL, transform=TRUE, normalize=TRUE) {
  
  if (!is.null(group)) { subset = paste0('_',group) } else { subset = '' }
  
  cat('\n------------------------------------------------------\n', 
      param[[1]],' ~ ',param[[3]], '   ', subset,
      '\n------------------------------------------------------\n')
  
  dc <- make_df(param[[1]], param[[2]], param[[3]], param[[4]], 
                group = group,
                transform=transform,
                normalize=normalize, 
                return_time_scale=time_constraints[2],
                verbose=verbose)

  times <- dc[c('med_ages1','med_ages2')]
  print(times)
  
  group_levels <- levels(dc$sample[,group])
  print(group_levels)
  
  set.seed(310896)
  
  start <- Sys.time(); cat(' started at: ', as.character(start), '\n\n') 
  # Run model
  m <- lavaan::lavaan(riclpm_formula(meas_time = times, 
                                     rel=time_constraints[1],
                                     stationarity = stationarity,
                                     strata=group_levels,
                                     verbose=verbose),
                      data = dc$sample, 
                      se = 'robust',
                      missing = 'fiml', 
                      fixed.x = FALSE,
                      group = group, 
                      meanstructure = TRUE,
                      int.ov.free = TRUE)
  end <- Sys.time()
  
  cat('\nDone! Runtime: ', difftime(end, start, units = 'mins'),
      'mins\n------------------------------------------------------\n')
  
  if (!lavaan::lavInspect(m, 'converged')){
    # Return warning message ( overwrites the model object )
    m <- paste(gsub('lavaan WARNING:', '', names(warnings())), collapse='\n')
    cat(m,'\n')
  } else {
    print(summary(m))
    
    # Fit measures
    fit_measures <- as.data.frame(lavaan::fitmeasures(m))
    # Parameters
    estimates_un <- lavaan::parameterEstimates(m) # unstandardized estimates
    standardized <- lavaan::standardizedSolution(m, type="std.all")[,c('est.std','ci.lower','ci.upper','pvalue')] # same but standardized

    correrations <- dc$corrs
    # save variable names and summary for dashboard output 
    data_summary <- dc$summary

    save(fit_measures, estimates_un, standardized, correrations, data_summary, times,
         file = paste0(out_folder,'/',substr(param[[1]],1,4),'_',param[[3]],subset,'.RData'))
  }
  
  return(m)
}

# ==============================================================================
# -------------------------------- Analysis ------------------------------------
# ==============================================================================

models <- list(
  # Self-reported depression ===================================================
  # - total fat mass/ FMI : 10/11 – 12/13 – 14 – 15½ / 16½ - 18 – 24
  list('sDEP_score', c(10.6, 12.8, 13.8, 16.6, 17.8, 23.8),
       'FMI',        c( 9.8, 11.8, 13.8, 15.4, 17.8, 24.5)),
  
  # – BMI: 10 – 13 – 14 – 16 - 18 - 24 # (with mean BMI 15.5 – 18)?
  list('sDEP_score', c(10.6, 12.8, 13.8, 16.6, 17.8, 23.8),
       'BMI',        c(10.7, 12.8, 13.8, 16,   17.8, 24.5)), # there is 17 too
  
  # - total lean mass/ LMI : 10/11 – 12/13 – 14 – 15½ / 16½ - 18 – 24
  list('sDEP_score', c(10.6, 12.8, 13.8, 16.6, 17.8, 23.8),
       'LMI',        c( 9.8, 11.8, 13.8, 15.4, 17.8, 24.5)),
  
  list('sDEP_score', c(10.6, 12.8, 13.8, 16.6, 17.8, 23.8),
       'total_fatmass', c( 9.8, 11.8, 13.8, 15.4, 17.8, 24.5)),
  
  list('sDEP_score', c(10.6, 12.8, 13.8, 16.6, 17.8, 23.8),
       'total_leanmass', c( 9.8, 11.8, 13.8, 15.4, 17.8, 24.5)),
  
  # # – waist circumference: 10 – 13 – 16 – 14/25
  list('sDEP_score', c(10.6, 12.8, 16.6, 23.8),
       'waist_circ', c(10.6, 12.8, 15.4, 24.5)),
  
  # Android fat mass 
  list('sDEP_score', c(13.8, 16.6, 17.8, 23.8),
       'android_fatmass', c(13.8, 15.4, 17.8, 24.5)),
  
  # Metabolic markers 
  list('sDEP_score', c(10.6, 16.6, 17.8, 23.8),
       'tot_chol',   c( 9.8, 15.4, 17.8, 24.5)),
  
  list('sDEP_score', c(10.6, 16.6, 17.8, 23.8),
       'HDL_chol',   c( 9.8, 15.4, 17.8, 24.5)),
  
  list('sDEP_score', c(10.6, 16.6, 17.8, 23.8),
       'LDL_chol',   c( 9.8, 15.4, 17.8, 24.5)),
  
  list('sDEP_score', c(10.6, 16.6, 17.8, 23.8),
       'insulin',    c( 9.8, 15.4, 17.8, 24.5)),
  
  list('sDEP_score', c(10.6, 16.6, 17.8, 23.8),
       'triglyc',    c( 9.8, 15.4, 17.8, 24.5)),
  
  list('sDEP_score', c(10.6, 16.6, 17.8, 23.8),
       'CRP',        c( 9.8, 15.4, 17.8, 24.5)),
  
  # Blood pressure, PWV and heart rate and glucose have three occasions...
  
  # Mother reported ============================================================
  # – total fat mass/ FMI : 10 – 12 – 13.5 – 17.5 (or with mean fm 15.5 – 18)
  list('mDEP_score', c(9.6, 11.7, 13.1, 16.7),
       'FMI',        c(9.8, 11.8, 13.8, 17.8)),
  
  # – BMI: 10 – 12 – 13 – 16  # (with mean BMI 15.5 – 18) 
  list('mDEP_score', c(9.6, 11.7, 13.1, 16.7),
       'BMI',        c(9.8, 11.8, 12.8, 16)), # there is 17 too but lower correlation (0.2 vs. 0.4)
  
  # – total lean mass/ LMI : 10 – 12 – 13.5 – 17.5 (or with mean lm 15.5 – 18)
  list('mDEP_score', c(9.6, 11.7, 13.1, 16.7),
       'LMI',        c(9.8, 11.8, 13.8, 17.8)),
  
  list('mDEP_score', c(9.6, 11.7, 13.1, 16.7),
       'total_fatmass', c(9.8, 11.8, 13.8, 17.8)),
  
  list('mDEP_score', c(9.6, 11.7, 13.1, 16.7),
       'total_leanmass', c(9.8, 11.8, 13.8, 17.8)),
  
  # – waist circumference: 10 – 12 – 13 – 16
  list('mDEP_score', c(9.6, 11.7, 13.1, 16.7),
       'waist_circ', c(9.8, 11.8, 12.8, 15.4))
)

# Run each depression / CMR marker combination in parallel =====================
foreach(i=1:length(models)) %dopar% {
  
  params = models[[i]]
  
  mtot = run_lavaan(param=params)
  
  # Stratify by sex
  if (!stationarity) { 
    msex = run_lavaan(param=params, group='sex') 
    lavTestLRT(mtot, msex)
  }
}

# ==============================================================================
# DESCRIBE
# summary(df)
# 
# df = df[!is.na(df$sex) & df$sex=='Female',]
# 
# for (v in names(df)){
#   if (startsWith(v, 'dep')) {
#     cat(round(median(df[,v], na.rm=TRUE), 1), ' (',
#         round(min(df[,v], na.rm=TRUE), 1),'–',
#         round(max(df[,v], na.rm=TRUE), 1),
#         ') [', round(sum(is.na(df[,v]))/nrow(df)*100),
#         '%]\n', sep='')
#   }
# }
# for (v in names(df)){
#   if (startsWith(v, 'cmr')) {
#     cat(round(median(df[,v], na.rm=TRUE), 1), ' (',
#         round(min(df[,v], na.rm=TRUE), 1),'–',
#         round(max(df[,v], na.rm=TRUE), 1),
#         ') [', round(sum(is.na(df[,v]))/nrow(df)*100),
#         '%]\n', sep='')
#   }
# }
# 
# t(stats::addmargins(t(table(df$sex, useNA = 'ifany')),
#                     margin=1,
#                     FUN=list(list(Percent = function(x){ paste0(round((x/nrow(df))*100),'%') } ))))
# 
# t(stats::addmargins(t(table(df$ethnicity, useNA = 'ifany')),
#                     margin=1,
#                     FUN=list(list(Percent = function(x){ paste0(round((x/nrow(df))*100),'%') } ))))
# 
# t(stats::addmargins(t(table(df$m_edu, useNA = 'ifany')),
#                     margin=1,
#                     FUN=list(list(Percent = function(x){ paste0(round((x/nrow(df))*100),'%') } ))))
