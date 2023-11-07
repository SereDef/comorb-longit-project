# ==============================================================================
# ============ 1. RANDOM-INTERCEPT CROSS-LAG PANEL MODEL (GCLPM) ===============
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

stationarity = !(grepl('pfree', out_folder)) # model long term effects stationary over time 

# Load dependencies
invisible(lapply(c('lavaan','tidySEM','foreach'), require, character.only = TRUE));
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
                    var1='dep',var2='cmr', verbose=TRUE) {
  
  d1 = data[,sel(name1, times1)]; d2 = data[,sel(name2, times2)]
  
  if(ncol(d1)!=ncol(d2)) { message('Error: unequal number of timepoints') }
  
  # Display characteristics of selection 
  temp_distance <- function(times) { tds <- c()
  for (i in 1:length(times)-1) { tds = c(tds, times[i+1]-times[i]) }
  return(tds)
  }
  
  names(d1) <- paste0(var1,1:ncol(d1))
  names(d2) <- paste0(var2,1:ncol(d2))
  
  df <- cbind(d1,d2)
  
  # select sample (at least one observation)
  if (any(rowSums(is.na(df)) != ncol(df))) { message ('Removing empty rows.')}
  samp = df[rowSums(is.na(df)) != ncol(df), ]
  
  if (verbose) {
    # Display temporal structure
    cat('Distance between constructs:', abs(times1-times2), sep='\t')
    cat('\nDepression, temporal gap:', temp_distance(times1), sep='\t')
    cat('\nCMR marker, temporal gap:', temp_distance(times2), sep='\t')
    
    # Display correlation
    message('\nCorrelations:')
    c = cor(df, use='pairwise.complete.obs', method='spearman')
    print(round(c,2))
    
    # Display sample size
    message('Sample size: ', nrow(samp))
  }
  
  return(samp)
}

## Specify the formula for the model (used in run_sem)
riclpm_formula <- function(var1='dep', var2='cmr', n_ocs=NULL, verbose=TRUE,
                           meas_time=list(c(), c()), stationarity=TRUE) {
  
  # Input check ----------------------------------------------------------------
  if (is.null(n_ocs) & length(meas_time[[1]])==0) { # user did not provide number or timepoints
    stop('Provide number of occasions or time points of measurement!') }
  
  if (length(meas_time[[1]]) != length(meas_time[[2]])) { # contraddictory info 
    stop('You need the same number of time points for each of the two variables.') }
  
  if (!is.null(n_ocs) & !length(meas_time[[1]]) %in% c(0, n_ocs)) { # contraddictory info 
    stop('Number of occations and measurement times provided do not agree.') }
  
  # ----------------------------------------------------------------------------
  # How many occasions 
  if (length(meas_time[[1]])>0) { 
    temp_var1 = meas_time[[1]]
    temp_var2 = meas_time[[2]] 
  } else { temp_var1 = temp_var2 = 1:n_ocs }
  
  if (is.null(n_ocs)) { n_ocs = length(temp_var1) }
  
  # Create between components (random intercepts)
  random_intercepts = paste0('# Random intercepts\n',
                             'ri_',var1,' =~ ', paste0('1*',var1,1:n_ocs, collapse=' + '), '\n',
                             'ri_',var2,' =~ ', paste0('1*',var2,1:n_ocs, collapse=' + '), '\n')
  
  # Create within-person centered variables
  impulses = paste0('# Impulses\n',
                    paste0('w_',var1,1:n_ocs,' =~ 1*',var1,1:n_ocs, collapse='\n'), '\n',
                    paste0('w_',var2,1:n_ocs,' =~ 1*',var2,1:n_ocs, collapse='\n'), '\n')
  
  # Estimate lagged effects between within-person centered variables
  regressions = '# Regressions\n'
  for (i in 2:n_ocs) { regressions = paste0(regressions, 
                                            'w_',var1,i,' ~ AR_',var1,i-1,' * w_',var1,i-1,' + CL_',var1,i-1,' * w_',var2,i-1,'\n',
                                            'w_',var2,i,' ~ AR_',var2,i-1,' * w_',var2,i-1,' + CL_',var2,i-1,' * w_',var1,i-1,'\n') 
  }
  
  # Estimate covariance between within-person centered variables at first wave
  # & between residuals of within-person centered variables (i.e., innovations)
  covariances = paste0('# Covariances\n',
                       paste0('w_',var1,1:n_ocs,' ~~ cov',1:n_ocs,' * w_',var2,1:n_ocs, collapse='\n'), '\n', 
                       # Estimate covariance of the random intercepts
                       'ri_',var1,' ~~ covRI * ri_',var2,'\n')
  
  # Estimate (residual) variance of within-person centered variables
  # note: residual variables from the second occasion onward
  variances = paste0('# Variances\n',
                     paste0('w_',var1,1:n_ocs,' ~~ w_',var1,1:n_ocs, collapse='\n'),'\n', 
                     paste0('w_',var2,1:n_ocs,' ~~ w_',var2,1:n_ocs, collapse='\n'),'\n', 
                     # Estimate variance of random intercepts
                     'ri_',var1,' ~~ ri_',var1,'\nri_',var2,' ~~ ri_',var2, '\n')
  
  ar_var1_con = ar_var2_con = cl_var1_con = cl_var2_con = ''
  
  if (stationarity) {
    rel='*' # '*'
    for (i in 1:(n_ocs-1)) {
      for (o in c(1:(n_ocs-1))[-i]) {
        
        ar_var1_con = paste0(ar_var1_con, 'AR_',var1,i,rel,round(temp_var1[i+1]-temp_var1[i],2),
                             ' == AR_',var1,o,rel,round(temp_var1[o+1]-temp_var1[o],2),'\n')
        ar_var2_con = paste0(ar_var2_con, 'AR_',var2,i,rel,round(temp_var2[i+1]-temp_var2[i],2),
                             ' == AR_',var2,o,rel,round(temp_var2[o+1]-temp_var2[o],2),'\n')
        
        cl_var1_con = paste0(cl_var1_con, 'CL_',var1,i,rel,round(temp_var1[i+1]-temp_var2[i],2),
                             ' == CL_',var1,o,rel,round(temp_var1[o+1]-temp_var2[o],2),'\n')
        cl_var2_con = paste0(cl_var2_con, 'CL_',var2,i,rel,round(temp_var2[i+1]-temp_var1[i],2),
                             ' == CL_',var2,o,rel,round(temp_var2[o+1]-temp_var1[o],2),'\n')
      }
    }
  }
  
  constraints = paste0('# Constraints\n', ar_var1_con, ar_var2_con, cl_var1_con, cl_var2_con)
  
  f = paste0(random_intercepts, impulses, regressions, covariances, variances, 
             constraints)
  
  if (verbose) { cat(f) }
  
  return(f)
}

## Fit a single SEM model (used in run_all_models)
run_lavaan <- function(params, normalize=TRUE) {
  
  d <- make_df(param[[1]], param[[2]], param[[3]], param[[4]], verbose=TRUE)
  
  if (normalize) { 
    minmax_norm <- function(x, ...) { return( (x - min(x,...)) / (max(x,...) - min(x,...)) )}
    d <- sapply(d, minmax_norm, na.rm = TRUE) 
  }
  
  cat('\n------------------------------------------------------\n',
      param[[1]],' ~ ',param[[3]])
  cat('\n------------------------------------------------------\nFitting the models...')
  
  set.seed(310896)
  
  start <- Sys.time(); cat(' started at: ', as.character(start), '\n') 
  # Run model
  m <- lavaan::lavaan(riclpm_formula(meas_time=list(param[[2]], param[[4]]), 
                                     stationarity=stationarity), 
                data = d, se = 'robust',
                missing = 'fiml', 
                meanstructure = TRUE, 
                int.ov.free = TRUE
  )
  end <- Sys.time()

  cat('\nDone! Runtime: ',round(end - start, 2),
      ' hours\n------------------------------------------------------\n')
  
  if (!lavInspect(m, "converged")){
    # Return warning message ( overwrites the model object )
    m <- paste(gsub('lavaan WARNING:', '', names(warnings())), collapse='\n') 
    cat(m,'\n')
    
  } else { 
    # Fit measures
    fm <- as.data.frame(lavaan::fitmeasures(m))
    # Parameters
    stad_est <- lavaan::standardizedSolution(m, type="std.all")[,'est.std'] # same but standardized
    estimates <- lavaan::parameterEstimates(m)
    
    save(fm, stad_est, estimates, file = paste0(out_folder,'/',substr(param[[1]],1,4),'_',param[[3]],'.RData'))
  }
  
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
  param = models[[i]]
  run_lavaan(params=param)
}
