# ==============================================================================
# =========== 1. RANDOM-INTERCEPT CROSS-LAG PANEL MODEL (RI-CLPM) ==============
# ==============================================================================

# args = commandArgs(trailingOnly = TRUE) # provide access to a copy of the command line arguments supplied
# 
# if (length(args) == 0) {
#   stop("Supply output folder name!")
# } else {
#   out_folder <- args[1] # e.g. 'mod1_ri_pstat' # parameters stationary
# }

# Set parameters
dir.create(out_folder)

stationarity <- !(grepl('pfree', out_folder)) # model long term effects stationary over time 

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
make_df <- function(name1, times1, name2, times2, strata=NULL,
                    var1='dep',var2='cmr', normalize=TRUE, verbose=TRUE) {
  
  # Select the two sets of measurement
  d1 = data[,sel(name1, times1)]
  d2 = data[,sel(name2, times2)]
  # Check input 
  if(ncol(d1)!=ncol(d2)) { message('Error: unequal number of timepoints') }
  
  # Rename variables
  names(d1) <- paste0(var1,1:ncol(d1))
  names(d2) <- paste0(var2,1:ncol(d2))
  
  df <- cbind(d1,d2)
  
  # Normalize data 
  if (normalize) { 
    minmax_norm <- function(x, ...) { return( (x - min(x,...)) / (max(x,...) - min(x,...)) )}
    df <- sapply(df, minmax_norm, na.rm = TRUE) 
  }
  
  # Stratify
  if (!is.null(strata)) { 
    df <- df[data[strata[1]]==strata[2],]
    message('Stratified sample: ',nrow(df))
  }
  # Subset sample (at least one observation)
  if (any(rowSums(is.na(df)) != ncol(df))) message('Removing empty rows.')
  samp = df[rowSums(is.na(df)) != ncol(df), ]

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
    c = round(cor(samp, use='pairwise.complete.obs', method='spearman'),2)
    print(c)
    
    # Display sample size
    message('Sample size: ', nrow(samp))
  }
  
  # Finally add age in moths for integer operators
  age_vars1 <- data[,sel(paste0(toupper(var1),'_age'), times1)] * 12 # back to months
  age_vars2 <- data[,sel(paste0(toupper(var2),'_age'), times2)] * 12 # back to months
  
  med_ages1 <- unname(round(sapply(age_vars1, median, na.rm=TRUE)))
  med_ages2 <- unname(round(sapply(age_vars2, median, na.rm=TRUE)))
  
  return(list(samp,c, med_ages1, med_ages2))
}

## Specify the formula for the model (used in run_sem)
riclpm_formula <- function(var1='dep', var2='cmr', 
                           n_ocs=NULL, meas_time=list(c(), c()),
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
    temp_var1 = meas_time[[1]]
    temp_var2 = meas_time[[2]] 
  } else { 
    temp_var1 = temp_var2 = 1:n_ocs 
  }
  
  if (is.null(n_ocs)) { n_ocs = length(temp_var1) }
  
  # Create between components (random intercepts)
  random_intercepts = paste0('# Random intercepts\n',
                             'ri_',var1,' =~ ', paste0('1*',var1,1:n_ocs, collapse=' + '), '\n',
                             'ri_',var2,' =~ ', paste0('1*',var2,1:n_ocs, collapse=' + '), '\n')
  
  # Create within-components (or within-person centered variables)
  # Note: factor loadings are set to 1 in non-stationary model and freely estimated in stationary models
  if (stationarity) { l = 'NA*' } else { l = '1*' }
  
  impulses = paste0('# Impulses\n',
                    paste0('w_',var1,1:n_ocs,' =~ ',l,var1,1:n_ocs, collapse='\n'), '\n',
                    paste0('w_',var2,1:n_ocs,' =~ ',l,var2,1:n_ocs, collapse='\n'), '\n')
  
  # Estimate lagged effects between within-person centered variables
  regressions = '# Regressions\n'
  for (i in 2:n_ocs) { 
    # if (n_strata > 1) { ... ar_1 <- paste0('c(', paste0(rep(paste0('AR_',var1,i-1),strata), collapse=', '), ')') 
    ar_1 <- paste0('AR_',var1,i-1)
    ar_2 <- paste0('AR_',var2,i-1)
    cl_1 <- paste0('CL_',var1,i-1)
    cl_2 <- paste0('CL_',var2,i-1)
    
    regressions = paste0(regressions, 
      'w_',var1,i,' ~ ',ar_1,' * w_',var1,i-1,' + ',cl_1,' * w_',var2,i-1,'\n',
      'w_',var2,i,' ~ ',ar_2,' * w_',var2,i-1,' + ',cl_2,' * w_',var1,i-1,'\n') 
  }
  
  # Estimate covariance between within-person centered variables 
  covariances = paste0('# Covariances\n',
                       # Estimate correlation between within-components at first wave
                       'w_',var1,'1 ~~ cor1 * w_',var2,'1\n',
                       # Estimate residual covariances (between residuals of within-person vars, i.e., innovations)
                       paste0('w_',var1,2:n_ocs,' ~~ rcov',2:n_ocs,' * w_',var2,2:n_ocs, collapse='\n'), '\n', 
                       # Estimate covariance of the random intercepts
                       'ri_',var1,' ~~ covRI * ri_',var2,'\n')
  
  if (stationarity) { w1 = '1*' } else { w1 = '' }
  
  variances = paste0('# Variances\n',
                     # Set variances of within-components at first wave to 1
                     paste0('w_',var1,'1 ~~ ',w1,'w_',var1,'1\n',
                            'w_',var2,'1 ~~ ',w1,'w_',var2,'1\n'),
                     # Estimate (and label) residual variances of within-person centered variables
                     paste0('w_',var1,2:n_ocs,' ~~ rvar_',var1,2:n_ocs,'*w_',var1,2:n_ocs, collapse='\n'),'\n', 
                     paste0('w_',var2,2:n_ocs,' ~~ rvar_',var2,2:n_ocs,'*w_',var2,2:n_ocs, collapse='\n'),'\n', 
                     # Estimate variance of random intercepts
                     'ri_',var1,' ~~ ri_',var1,'\nri_',var2,' ~~ ri_',var2, '\n')
  
  constraints = ar_var1_con = ar_var2_con = cl_var1_con = cl_var2_con = cor_con = rvar_con = ''
  
  if (stationarity) {
    constraints = '# Constraints\n'
      
    # paste0('# Constrain grand means over time\n',
    # paste0(var1,1:n_ocs, collapse=' + '), ' ~ mean_',var1,'*1\n',
    # paste0(var2,1:n_ocs, collapse=' + '), ' ~ mean_',var2,'*1\n')
    
    # Weight for unequal temporal gaps 
    rel='^' # '*'
    i = 1:(n_ocs-2)
    # round to nearest integer?
    prec=2
    ar_var1_con <- paste0('AR_',var1,i,rel,round(temp_var1[i+1]-temp_var1[i], prec),
                          ' == AR_',var1,i+1,rel,round(temp_var1[i+2]-temp_var1[i+1], prec),
                          collapse='\n')
    ar_var2_con <- paste0('AR_',var2,i,rel,round(temp_var2[i+1]-temp_var2[i], prec), 
                          ' == AR_',var2,i+1,rel,round(temp_var2[i+2]-temp_var2[i+1], prec),
                          collapse='\n')
    cl_var1_con <- paste0('CL_',var1,i,rel,round(temp_var1[i+1]-temp_var2[i], prec), 
                          ' == CL_',var1,i+1,rel,round(temp_var1[i+2]-temp_var2[i+1], prec),
                          collapse='\n')
    cl_var2_con <- paste0('CL_',var2,i,rel,round(temp_var2[i+1]-temp_var1[i], prec), 
                          ' == CL_',var2,i+1,rel,round(temp_var2[i+2]-temp_var1[i+1], prec),
                          collapse='\n')
    
    # Compute correlations between the within-components themselves at each wave, 
    i = 1:(n_ocs-1) # iterator
    
    cor_con <- paste0('cor',2:n_ocs,
                      ' := AR_',var1,i,'*CL_',var2,i,' + ',
                          'CL_',var1,i,'*AR_',var2,i,' + ',
                 'AR_',var1,i,'*AR_',var2,i,'*cor',i,' + ',
                 'CL_',var1,i,'*CL_',var2,i,'*cor',i,' + rcov',2:n_ocs, collapse='\n')
    
    # Constrain residual variances of within-components such that the total variance of each 
    # within-component equals 1
    rvar_con <- paste0('rvar_',var1,2:n_ocs,' == 1 - (AR_',var1,i,'*AR_',var1,i,' + ',
                       'CL_',var1,i,'*CL_',var1,i,' + 2*AR_',var1,i,'*CL_',var1,i,'*cor',i,')\n',
                       'rvar_',var2,2:n_ocs,' == 1 - (AR_',var2,i,'*AR_',var2,i,' + ',
                       'CL_',var2,i,'*CL_',var2,i,' + 2*AR_',var2,i,'*CL_',var2,i,'*cor',i,')',
                       collapse = '\n')
  }
  
  constraints = paste0(constraints, 
                       ar_var1_con,'\n', ar_var2_con,'\n', cl_var1_con,'\n', cl_var2_con, 
                       '\n', cor_con,'\n',rvar_con)
  
  f = paste0(random_intercepts, impulses, regressions, covariances, variances, 
             constraints)
  
  if (verbose) { cat(f) }
  
  return(f)
}

## Fit a single SEM model (used in run_all_models)
run_lavaan <- function(param, group=NULL, normalize=TRUE) {
  
  if (!is.null(group)) { subset = paste(c('',group), collapse='_') } else { subset = '' }
  
  cat('\n------------------------------------------------------\n', 
      param[[1]],' ~ ',param[[3]], '   ', subset,
      '\n------------------------------------------------------\nFitting the models...')
  
  dc <- make_df(param[[1]], param[[2]], param[[3]], param[[4]], strata = group, 
               normalize=normalize, verbose=TRUE)
  d <- dc[[1]] # Dataframe
  c <- dc[[2]] # Correlation matrix
  
  set.seed(310896)
  
  start <- Sys.time(); cat(' started at: ', as.character(start), '\n\n') 
  # Run model
  m <- lavaan::lavaan(riclpm_formula(meas_time = dc[3:4], 
                                     stationarity = stationarity),
                      data = d,
                      se = 'robust',
                      missing = 'fiml', # group = group, 
                      meanstructure = TRUE,
                      int.ov.free = TRUE
                      )
  end <- Sys.time()
  
  cat('\nDone! Runtime: ', difftime(end, start, units = 'mins'),
      'mins\n------------------------------------------------------\n')
  
  # if (!lavaan::lavInspect(m, 'converged')){
  #   # Return warning message ( overwrites the model object )
  #   m <- paste(gsub('lavaan WARNING:', '', names(warnings())), collapse='\n') 
  #   cat(m,'\n')
  #   
  # } else { 
  #   # Fit measures
  #   fm <- as.data.frame(lavaan::fitmeasures(m))
  #   # Parameters
  #   stad_esti <- lavaan::standardizedSolution(m, type="std.all")[,'est.std'] # same but standardized
  #   estimates <- lavaan::parameterEstimates(m)
  #   
  #   save(c, fm, estimates, stad_esti, file = paste0(out_folder,'/',
  #                                     substr(param[[1]],1,4),'_',param[[3]],subset,'.RData'))
  # }
  
  return(m)
}

# ==============================================================================
# -------------------------------- Analysis ------------------------------------
# ==============================================================================

models <- list(
  # Self-reported depression ===================================================
  # – BMI: 10 – 13 – 14 – 16 - 18 - 24 # (with mean BMI 15.5 – 18)?
  list('sDEP_score', c(10.6, 12.8, 13.8, 16.6, 17.8, 23.8),
       'BMI',        c(10.7, 12.8, 13.8, 16,   17.8, 24.5)), # there is 17 too
  
  # - total fat mass/ FMI : 10/11 – 12/13 – 14 – 15½ / 16½ - 18 – 24
  list('sDEP_score', c(10.6, 12.8, 13.8, 16.6, 17.8, 23.8),
       'FMI',        c( 9.8, 11.8, 13.8, 15.4, 17.8, 24.5)),
  
  list('sDEP_score', c(10.6, 12.8, 13.8, 16.6, 17.8, 23.8),
       'total_fatmass', c( 9.8, 11.8, 13.8, 15.4, 17.8, 24.5)),
  
  # - total lean mass/ LMI : 10/11 – 12/13 – 14 – 15½ / 16½ - 18 – 24
  list('sDEP_score', c(10.6, 12.8, 13.8, 16.6, 17.8, 23.8),
       'LMI',        c( 9.8, 11.8, 13.8, 15.4, 17.8, 24.5)),
  
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
  
  # Mother reported ==============================================================
  # – BMI: 10 – 12 – 13 – 16  # (with mean BMI 15.5 – 18) 
  list('mDEP_score', c(9.6, 11.7, 13.1, 16.7),
       'BMI',        c(9.8, 11.8, 12.8, 16)), # there is 17 too but lower correlation (0.2 vs. 0.4)
  
  # – total fat mass/ FMI : 10 – 12 – 13.5 – 17.5 (or with mean fm 15.5 – 18)
  list('mDEP_score', c(9.6, 11.7, 13.1, 16.7),
       'FMI',        c(9.8, 11.8, 13.8, 17.8)),
  
  list('mDEP_score', c(9.6, 11.7, 13.1, 16.7),
       'total_fatmass', c(9.8, 11.8, 13.8, 17.8)),
  
  # – total lean mass/ LMI : 10 – 12 – 13.5 – 17.5 (or with mean lm 15.5 – 18)
  list('mDEP_score', c(9.6, 11.7, 13.1, 16.7),
       'LMI',        c(9.8, 11.8, 13.8, 17.8)),
  
  list('mDEP_score', c(9.6, 11.7, 13.1, 16.7),
       'total_leanmass', c(9.8, 11.8, 13.8, 17.8)),
  
  # – waist circumference: 10 – 12 – 13 – 16
  list('mDEP_score', c(9.6, 11.7, 13.1, 16.7),
       'waist_circ', c(9.8, 11.8, 12.8, 15.4))
)

# Run each depression / CMR marker combination in parallel =====================
foreach(i=1:length(models)) %dopar% {
  
  params = models[[i]]
  run_lavaan(param=params)
  
  # Stratify by sex
  run_lavaan(param=params, group= c('sex','1'))
  run_lavaan(param=params, group= c('sex','2'))
}

m = run_lavaan(param=models[[2]])

e = stad_esti[stad_esti$label!='',]

# AR_dep1^27== AR_dep2^12
# AR_dep2^12== AR_dep3^33
# AR_dep3^33== AR_dep4^14
# AR_dep4^14== AR_dep5^73
# AR_cmr1^23== AR_cmr2^25
# AR_cmr2^25== AR_cmr3^19
# AR_cmr3^19== AR_cmr4^28
# AR_cmr4^28== AR_cmr5^81

# CL_dep1^36== CL_dep2^25
# CL_dep2^25== CL_dep3^33
# CL_dep3^33== CL_dep4^28
# CL_dep4^28== CL_dep5^73

# CL_cmr1^14== CL_cmr2^12
# CL_cmr2^12== CL_cmr3^19
# CL_cmr3^19== CL_cmr4^14
# CL_cmr4^14== CL_cmr5^81
