# ==============================================================================
# ============================ 1. G-ESTIMATION =================================
# ==============================================================================
# Advantages of g-estimation: 
# 1. no models for the (distribution of) covariates are required, so these can be 
#    either continuous or not, and be either time-invariant or time-varying 
#    (thus possibly affected by exposure). 
# 2. unbiased estimation requires correctly specifying a model for the exposure 
#    or the outcome, but not necessarily both. 

# NOTE: assuming a simple additive SNMM = parameters encoding causal effects do not 
# depend (i.e. are not moderated or modified by) any other variable.

# Time varying confounders: pubertal stage and physical exercise
# Do I believe these are the only ones exposure-confounder feedback over time? 

# ------------------------------------------------------------------------------
# args = commandArgs(trailingOnly = TRUE) # provide access to a copy of the command line arguments supplied
# 
# if (length(args) == 0) {
#   stop("Supply output folder name!")
# } else {
#   out_folder <- args[1] # e.g. 'mod1_ri_pstat' # parameters stationary
# }

# Set parameters
out_folder <- '../tmp_g_estim'
dir.create(out_folder)

# Load dependencies
require(lavaan)
# library(ALSPAC.helpR)

# Critical advantages of using lavaan for g-estimation:
# 1. FIML 
# 2. Bootstrapped SEs
# 3. Fit measures
# 4. Equality constraints on the (regression) model parameters can be used for
#    Sensitivity analysis for unmeasured confounding
# can all be easily implemented during estimation.
# PROBLEM: glm for propensity scores for now set to na.exclude (can use lavaan for fiml)
# PROBLEM 2 calculating the residual scores accounting for NAs...

# Read in data
data <- readRDS('../mats/raw_data.rds')

## Select subsets of data (i.e., which variables and timepoints)
sel <- function(var, times=NULL) { 
  subs = names(data)[grep(paste(var, collapse='|'), names(data))]
  if (!is.null(times)) { subs = subs[grep(paste(paste0('_',times), collapse='|'), subs)] }
  return(subs)
}

# ==============================================================================
g_estimate <- function(tvar, tinv, Y_name, X_name, L_name,
                       normalize=TRUE, prop_score_family = 'gaussian',
                       SE_method='robust', # change to 'bootstrap'
                       verbose=TRUE)
  {
  # Number of observations
  n_ocs <- nrow(tvar)
  
  # Step 0: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Construct dataframe (with renamed and normalized variables) and compute the 
  # propensity scores for the exposure
  
  # NOTE: assumes the order: dep, cmr, confounder (time-varying)
  if (Y_name=='dep') { y=1; x=2 } else { y=2; x=1 }
  
  Y = data[, sel(names(tvar)[y], tvar[,y])]
  X = data[, sel(names(tvar)[x], tvar[,x])]
  L = data[, sel(names(tvar)[3], tvar[,3])]
  C = data[, tinv]
    
  # if(ncol(Y)!=ncol(X) | ncol(Y)!=ncol(L)) { stop('Unequal number of timepoints') }
    
  names(Y) <- paste0(Y_name,1:ncol(Y))
  names(X) <- paste0(X_name,1:ncol(X))
  names(L) <- paste0(L_name,1:ncol(L))
    
  df <- cbind(Y,X,L,C)
    
  # select sample (at least one observation)
  if (verbose & any(rowSums(is.na(df)) != ncol(df))) { 
      message ('Removing ',sum(rowSums(is.na(df)) == ncol(df)),' empty rows.') 
  }
  samp = df[rowSums(is.na(df)) != ncol(df), ]
    
  # Display sample size
  message('Sample size: ', nrow(samp))
    
  if (verbose) {
    # Display correlation
    message('\nCorrelations:')
    cormat = cor(sapply(df, as.numeric), use='pairwise.complete.obs', method='spearman')
    print(round(cormat,2))
  }
  
  if (normalize) { 
    minmax_norm <- function(x, ...) { 
      x = as.numeric(x)
      return( (x - min(x,...)) / (max(x,...) - min(x,...)) )
    }
    samp <- data.frame(sapply(samp, minmax_norm, na.rm = TRUE))
  } else { samp <- data.frame(sapply(samp, as.numeric)) }
  
  
  # Compute propensity scores (PS) for the exposures
  message('\nComputing propensity scores for the exposure...')
  for (x in (ncol(X)-1):1){
    # For now, only one more lag is considered
    if (x > 1) { lag1 <- paste0(X_name,x-1,' + ',L_name,x-1,' + ',Y_name,x-1,' + ') 
    } else { lag1 <- '' }
    
    f <- paste0(X_name,x,' ~ ',lag1,L_name,x,' + ',Y_name,x,' + ',
                paste(tinv,collapse=' + '))
    if (verbose) { cat(f, '\n') }
    
    PSfit <- lavaan::sem(f, data=samp, missing='fiml', fixed.x=FALSE)
    samp[, paste0('PS_',X_name,x)] <- lavaan::lavPredict(PSfit, type='yhat')[, paste0(X_name,x)]
    
    # PSfit <- glm(f, data=samp, family=prop_score_family, na.action=na.exclude)
    # samp[, paste0('PS_',X_name,x)] <- predict(PSfit, type='response', newdata=samp)
    # browser()
  }
  cat('\n')
  if (verbose) { print(summary(samp[, grepl('PS_',names(samp))]))}
  
  # Define model formula to feed to lavaan
  make_formula <- function(lag=1, res_y = '', add_lag=TRUE) {
    # Save me some typing
    Y=Y_name; X=X_name; L=L_name
    
    f <- paste0('# Lag-',lag,' regressions\n')
    for (i in n_ocs:(lag+1)) {
      # For now, only one more lag is considered
      if (add_lag & i > lag+1) { more_lag <- paste0(' + ',X,i-(lag+1),' + ',Y,i-(lag+1),' + ',L,i-(lag+1)) 
      } else { more_lag = ''}
      
      f <- paste0(f, 
                  res_y,Y,i,' ~ psi',i,i-lag,' * ',X,i-lag,' + PS_',X,i-lag,' + ',Y,i-lag,' + ',L,i-lag, 
                  more_lag,' + ',paste0(T_invary, collapse= ' + '),'\n') 
    }
    
    if (verbose) { cat(f) }
    return(f)
  }
  
  # Step 1: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Estimate all lag-1 effects (ψ54, ... ψ21) by fitting a linear regression model
  # for each outcome time point (Y5, ... Y2) [i.e. in the same joint model]
  message('\nFitting lag-1 models...')
  mod1 = make_formula()
  fit1 <- lavaan::sem(mod1, data=samp, se=SE_method, # TODO: specify bootstrap instead?
                      missing='fiml', fixed.x=FALSE)
  # Inspect
  if (verbose) { print(summary(fit1)) } # fit.meas = TRUE
  # fitmeas1 <- as.data.frame(lavaan::fitmeasures(fit1)) # this takes a minute
  
  # Extract the g-estimators of these lag-1 effects (ψ54, ... ψ21)
  Gest1 <- lavaan::parameterEstimates(fit1)
  
  # Step 2: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Remove the lag-1 effect of X from Y: e.g. RES_Y3 = Y3 − ψ32_hat * X2.
  message('Computing lag-1 residualized outcomes...')
  for (i in n_ocs:3) {
    psi <-  Gest1[grepl(paste0('psi',i), Gest1$label), 'est']
    samp[, paste0('RES_',Y_name,i)] <- (samp[, paste0(Y_name,i)] - psi * samp[, paste0(Y_name,i-1)])
  }
  if (verbose) { print(summary(samp[, grepl('RES_',names(samp))])) }
  
  # Step 3: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Fit the lag-2 models for the residualised Ys (including only lag-2 covariates).
  message('\nFitting lag-2 models...')
  mod2 <- make_formula(lag=2, res_y = 'RES_', add_lag=FALSE)
  fit2 <- lavaan::sem(mod2, data=samp, se=SE_method, # TODO: specify bootstrap instead?
                      missing='fiml', fixed.x=FALSE)
  # Inspect
  if (verbose) { print(summary(fit2)) } # fit.meas = TRUE
  # fitmeas2 <- as.data.frame(lavaan::fitmeasures(fit2)) # this takes a minute
  
  # Extract the g-estimators of these lag-2 effects (ψ53, ψ42, ψ31)
  Gest2 <- lavaan::parameterEstimates(fit2)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Bind all g-estimators together and return 
  allGest <- rbind( Gest1[grepl('psi', Gest1$label), ],
                    Gest2[grepl('psi', Gest2$label), ])
  
  return(allGest)
}

# ==============================================================================
# Specify variables
# - total fat mass/ FMI : 10/11 – 12/13 – 13/14 – 15/16 - 17/18
T_vary <- data.frame('sDEP_score' = c(10.6, 12.8, 13.8, 16.6, 17.8), # 23.8
                            'FMI' = c( 9.8, 11.8, 13.8, 15.4, 17.8), # 24.5
                  'phys_activity' = c( 9.6, 11.7, 13.1, 15.3, 17.0)) # 8.2, 10.7, 14.7, 16.0

T_invary <- c('sex','m_edu','p_edu')
# TODO: selecting a set of pre-treatment covariates that suffice to eliminate all 
# non-causal associations

# ==============================================================================
# NOTE: we refer to the exposure as X, the outcome as Y and the time-varying covariate as L.

gs_dep = g_estimate(T_vary, T_invary, Y_name='dep', X_name='cmr', L_name='vpa', normalize = FALSE)

gs_cmr = g_estimate(T_vary, T_invary, Y_name='cmr', X_name='dep', L_name='vpa')


# ==============================================================================

# Standard errors and confidence intervals can be estimated using a nonparametric 
# percentile bootstrap procedure (don't use the s.e. for ψ31 from the fitted model 
# alone, because they do not adequately capture the variability in the estimates ψ32ˆ)

# ------------------------------------------------------------------------------
# A fundamental, empirically untestable assumption for valid causal inference is 
# the absence of unmeasured exposure-outcome confounding. 
# Sensitivity analysis for unmeasured confounding
# Unmeasured confounding can manifest in correlated exposure and outcome residual 
# errors. First, fix the residual correlations at a given value. Estimate the exposure
# effects under this fixed constraint. Repeat both steps using different fixed values 
# of the residual correlations: systematically investigate how the effect estimates 
# change depending on the given strength of unmeasured confounding. 
# Stronger correlations indicate more severe violations of the unconfoundedness assumptions, 
# whereas a zero correlation corresponds to the assumption of no unmeasured confounding. 
# The correlations can be readily parametrized in lavaan as follows.
# A sensitivity analysis then entails fixing the residual correlations to different 
# (non-zero) values by replacing rho by a value between -1 and 1, then carrying out 
# the g-estimation procedure for each fixed value to gauge how different the effect 
# estimates can be
model_lag1 <- '
Y3 ~ psi32*X2 + X1 + L2 + Y2 + L1 + Y1 + C
Y2 ~ psi21*X1 + L1 + Y1 + C
## correlated residuals
Y2 ~~ cov_yx*X1
Y3 ~~ cov_yx*X2
## residual variances
Y2 ~~ var_y2*Y2
Y3 ~~ var_y3*Y3
X1 ~~ var_x1*X1
X2 ~~ var_x2*X2
var_y2 > 0
var_y3 > 0
var_x1 > 0
var_x2 > 0
## fixed correlations
cov_yx/sqrt(var_y2*var_x1) == rho
cov_yx/sqrt(var_y3*var_x2) == rho
'