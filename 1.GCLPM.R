# ==============================================================================
# ============= 1. GENERALIZED CROSS-LAG PANEL MODEL (GCLPM) ===================
# ==============================================================================

args = commandArgs(trailingOnly = TRUE) # provide access to a copy of the command line arguments supplied

if (length(args) == 0) {
  stop("Supply output folder name!")
} else {
  out_folder <- args[1] # e.g. 'mod1_ri_pstat' # parameters stationary
}

# Set parameters
dir.create(out_folder)

fix_lambdas = !(grepl('lfree', out_folder)) # heterogeneity allowed to vary over time 
stationarity = !(grepl('pfree', out_folder)) # model long term effects stationary over time 

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

gclpm_model_matrix <- function(var1='dep', var2='cmr', AR_always = TRUE, 
                               save_csv=FALSE) {
  
  # All possible combinations of the 8 main parameters = 256 potential models 
  m = t(expand.grid(lapply(numeric(8), function(x) c(1, 0))))
  
  # Remove combinations that have less than 2 parameters in the model
  m = m[, colSums(m)>1] # = 247 potential models 
  
  # Order by number of parameters included and rename 
  params = c('ltAR','ltCL','maAR','maCL')
  param_names = c( paste(params, var1, sep='_'), paste(params, var2, sep='_'))
  
  mat = data.frame( m[,order(colSums(-m))], row.names = param_names)
  
  # Generate column names for model structure (indicating params which are excluded)
  nmat = mat 
  for (r in 1:nrow(mat)){ # TODO: can't make it work with apply (can't get row.name)
    nmat[r, ] <- replace(nmat[r, ], which(nmat[r, ]==0), rownames(nmat[r, ]))
  }
  nms = gsub('1-','', gsub('-1', '', lapply(nmat, paste0, collapse='-'))); nms[1]='full'
  names(mat) <- nms
  
  rm(m, nmat, nms, r)
  
  # Remove some comparisons, to speed up processing: 
  # e.g., always estimate long term AR terms
  if (AR_always) { mat <- mat[, -grep('ltAR', names(mat))] } # --> 64 models
  
  if (save_csv) { write.csv(mat, '../mats/model_structure.csv') }
  
  return(mat)
}

mat = gclpm_model_matrix()

# ------------------------------------------------------------------------------

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
gclpm_formula <- function(var1='dep', var2='cmr', 
                          n_ocs=NULL, meas_time=list(c(), c()), 
                          mod='full', 
                          rel='*', # specify temporal equality using power ('^') or linear ('*') relations 
                          strata=c(),
                          stationarity=FALSE,
                          fix_lambdas=FALSE, 
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
  
  # Define model structure according to mat matrix
  mat = gclpm_model_matrix(var1=var1, var2=var2)
  
  for (r in row.names(mat)) {
    varname <- tolower(gsub(var1, 'var1', gsub(var2, 'var2', r)))
    if (mat[mod][r,]==0) { assign(gsub('lt','', varname),'0') 
    } else { assign(gsub('lt','', varname), gsub('lt','', r)) }
  }
  
  # Create between components (unit effects) -----------------------------------
  unit_effects = paste0('# Unit effects\n',
                        'eta_',var1,' =~ ', paste0('l',n_ocs:1,'*',var1,n_ocs:1, collapse=' + '), '\n',
                        'eta_',var2,' =~ ', paste0('l',n_ocs:1,'*',var2,n_ocs:1, collapse=' + '), '\n') 
  
  # Create impulses (or within-person centered variables) ----------------------
  impulses = '# Impulses\n'
  for (i in 2:n_ocs) { impulses = paste0(impulses, 'u_',var1,i,' =~ ',var1,i,'\n',var1,i,' ~~ 0*',var1,i,'\n') }
  for (i in 2:n_ocs) { impulses = paste0(impulses, 'u_',var2,i,' =~ ',var2,i,'\n',var2,i,' ~~ 0*',var2,i,'\n') }
  
  # Estimate lagged effects between within-person centered variables -----------
  regressions <- '# Regressions\n'
  
  for (i in n_ocs:2) { 
    # assign different name to each parameter unless it was set to 0
    if (ar_var1!=0) { ar_var1i = paste0(ar_var1,i-1) } else { ar_var1i = ar_var1 }
    if (cl_var1!=0) { cl_var1i = paste0(cl_var1,i-1) } else { cl_var1i = cl_var1 }
    if (ar_var2!=0) { ar_var2i = paste0(ar_var2,i-1) } else { ar_var2i = ar_var2 }
    if (cl_var2!=0) { cl_var2i = paste0(cl_var2,i-1) } else { cl_var2i = cl_var2 }
    
    if (i > 2) { # No MA terms because no impulses at first occasion
      
      if (maar_var1!=0) { maar_var1i = paste0(maar_var1,i-1) } else { maar_var1i = maar_var1 }
      if (macl_var1!=0) { macl_var1i = paste0(macl_var1,i-1) } else { macl_var1i = macl_var1 }
      if (maar_var2!=0) { maar_var2i = paste0(maar_var2,i-1) } else { maar_var2i = maar_var2 }
      if (macl_var2!=0) { macl_var2i = paste0(macl_var2,i-1) } else { macl_var2i = macl_var2 }
      
      ma_terms_var1 = paste0(' + ',maar_var1i,'*u_',var1,i-1,' + ',macl_var1i,'*u_',var2,i-1) 
      ma_terms_var2 = paste0(' + ',maar_var2i,'*u_',var2,i-1,' + ',macl_var2i,'*u_',var1,i-1) 
      
    } else { ma_terms_var1 = ma_terms_var2 = '' }
    
    # Add depression regressions
    regressions = paste0(regressions, 
                         var1,i,' ~ ',ar_var1i,'*',var1,i-1,' + ',cl_var1i,'*',var2,i-1, ma_terms_var1,'\n',
                         var2,i,' ~ ',ar_var2i,'*',var2,i-1,' + ',cl_var2i,'*',var1,i-1, ma_terms_var2,'\n') 
  }
  
  comevement = '# Comovement\ndep1 ~~ comv1*cmr1\n' # first occasion no impulse
  for (i in 2:n_ocs) { comevement = paste0(comevement,'u_',var1,i,' ~~ comv',i,'*u_',var2,i,'\n')}
  
  restrictions = '# Restrictions\n'
  for (i in 2:n_ocs) { 
    if (i < n_ocs) { ar = paste0(' + 0*u_',var1, seq(i+1, n_ocs), collapse='') } else { ar = ''}
    restrictions = paste0(restrictions, 'u_',var1,i,' ~~ 0*eta_',var1,' + 0*eta_',var2, ar, ' + ',
                          paste0('0*u_',var2,c(2:n_ocs)[-(i-1)], collapse=' + '), '\n') }
  for (i in 2:n_ocs) {  
    if (i < n_ocs) { ar = paste0(' + 0*u_',var2, seq(i+1, n_ocs), collapse='') } else { ar = ''}
    restrictions = paste0(restrictions, 'u_',var2,i,' ~~ 0*eta_',var1,' + 0*eta_',var2, ar, '\n') }
  
  constraints = '# Constraints\n'
  
  if (fix_lambdas) {
    for (i in 1:n_ocs) { constraints = paste0(constraints, 'l',i,' == 1\n') }
  } # only last lambda is set to 1
  
  ar_var1_con = ar_var2_con = cl_var1_con = cl_var2_con = ''
  
  if (stationarity) {
    
    if (rel=='*') { prec=1 } else if (rel=='^') { prec=0 } 
    
    for (i in 1:(n_ocs-1)) {
      for (o in c(1:(n_ocs-1))[-i]) {
        
        ar_var1_con = paste0(ar_var1_con, 'AR_',var1,i,rel,round(temp_var1[i+1]-temp_var1[i],prec),
                             ' == AR_',var1,o,rel,round(temp_var1[o+1]-temp_var1[o],prec),'\n')
        ar_var2_con = paste0(ar_var2_con, 'AR_',var2,i,rel,round(temp_var2[i+1]-temp_var2[i],prec),
                             ' == AR_',var2,o,rel,round(temp_var2[o+1]-temp_var2[o],prec),'\n')
        
        if (cl_var1!=0) { cl_var1_con = paste0(cl_var1_con, 'CL_',var1,i,rel,round(temp_var1[i+1]-temp_var2[i],prec),
                                             ' == CL_',var1,o,rel,round(temp_var1[o+1]-temp_var2[o],prec),'\n') }
        if (cl_var2!=0) { cl_var2_con = paste0(cl_var2_con, 'CL_',var2,i,rel,round(temp_var2[i+1]-temp_var1[i],prec),
                                             ' == CL_',var2,o,rel,round(temp_var2[o+1]-temp_var1[o],prec),'\n') }
      }
      
    }
    constraints = paste0(constraints, ar_var1_con, ar_var2_con, cl_var1_con, cl_var2_con)
  }
  
  # Paste everything together 
  f = paste0(unit_effects, impulses, regressions, comevement, restrictions, constraints)
  
  if (verbose) { message('MODEL SYNTAX:'); cat(f) }
  
  return(f)
}

## Save model output as a graph in .png format (used in run_sem)
save_semgraph <- function(model, n_ocs, name) {
  
  args = list('eta_dep', '', '', '', '', '', 
              'dep1', 'dep2', 'dep3', 'dep4', 'dep5', 'dep6', 
              ' ', 'u_dep2', 'u_dep3', 'u_dep4', 'u_dep5', 'u_dep6', 
              ' ', 'u_cmr2', 'u_cmr3', 'u_cmr4', 'u_cmr5', 'u_cmr6', 
              'cmr1', 'cmr2', 'cmr3', 'cmr4', 'cmr5', 'cmr6', 
              'eta_cmr', '', '', '', '', '', rows=6)
  
  lay = rlang::exec(tidySEM::get_layout, !!!args)
  
  lay = as.matrix(lay[, 1:n_ocs])
  
  # Remove edges estimated or set to 0.00 for readability
  e = get_edges(model); e$show[e$label=='0.00'] <- FALSE
  # Create graph 
  p = tidySEM::graph_sem(model = model, layout = lay, edges = e)
  # Save image
  dir.create(paste0(out_folder,'_graphs'))
  ggplot2::ggsave(paste0(out_folder,'_graphs/',name,'.png'), p, device = "png", width= 20, height = 15)
}

## Fit a single SEM model (used in run_all_models)
run_sem <- function(mod, data, times, name1, name2, plot_semgraph=FALSE) {
  
  cat('- ', which(names(mat)==mod), mod)
  
  # Run model
  set.seed(310896)
  
  m <- lavaan::lavaan(gclpm_formula(meas_time = times, mod=mod,
                                    rel=time_constraints[1],
                                    stationarity = stationarity,
                                    fix_lambdas=fix_lambdas,
                                    # strata=group_levels,
                                    verbose=verbose),
                      data = data, 
                      se = 'robust',
                      missing = 'fiml', 
                      fixed.x = FALSE)
                      # group = group, 
                      # meanstructure = TRUE
                      # int.ov.free = TRUE
  
  if (lavInspect(m, "converged")){
    if (plot_semgraph) {
      # Print out the layout 
      save_semgraph(m, n_ocs=ncol(data)/2, name=paste(substr(name1,1,4),name2,mod, sep='_'))
    }
  } else {
    # Return warning message ( overwrites the model object )
    m <- paste(gsub('lavaan WARNING:', '', names(warnings())), collapse='\n') }
  
  return(m)
}

## Run the 64 SEM models in parallel :)
run_all_models <- function(param, which_models=names(mat), 
                           transform=TRUE,
                           normalize=TRUE,
                           verbose=verbose){
  
  cat('\n------------------------------------------------------\n', 
      param[[1]],' ~ ',param[[3]],
      '\n------------------------------------------------------\n')
  
  dc <- make_df(param[[1]], param[[2]], param[[3]], param[[4]], 
                # group = group,
                transform=transform,
                normalize=normalize, 
                return_time_scale=time_constraints[2],
                verbose=verbose)
  
  times <- dc[c('med_ages1','med_ages2')]
  print(times)
  
  # save variable names and summary for dashboard output 
  correrations <- dc$corrs
  data_summary <- dc$summary
  
  set.seed(310896)
  
  start <- Sys.time(); cat(' started at: ', as.character(start), '\n\n') 
  
  # Run this bitch (in parallel)
  fits <- mapply(run_sem, mod=which_models, # parallel::mc
                 MoreArgs=list(data=dc$data, 
                               times=times, 
                               name1=param[[1]], name2=param[[3]]), 
                 SIMPLIFY=FALSE) # mc.cores=12)
  end <- Sys.time()
  
  cat('\nDone! Runtime: ', difftime(end, start, units = 'mins'),
      'mins\n------------------------------------------------------\n')
  
  # Initialize output dataframes
  fit_measures = data.frame(matrix(nrow = 55, ncol = 0)); estimates = data.frame(); failed = data.frame()
  
  for (f in names(fits)){
    if( is.character(fits[[f]]) ) {
      message(f) # , '\n\t', fits[[f]]) 
      problem <- data.frame(fits[[f]], row.names=f); names(problem) <- 'problem'
      failed <- rbind(failed, problem)           
    } else { 
      # Inspect # print(lavaan::summary(fits[[f]]))
      
      # Fit measures
      fm <- as.data.frame(lavaan::fitmeasures(fits[[f]])); names(fm) <- f
      fit_measures <- cbind(fit_measures, fm)
      # Parameters
      estimates_un <- lavaan::parameterEstimates(fits[[f]]) # unstandardized estimates
      standardized <- lavaan::standardizedSolution(fits[[f]], type="std.all")[,c('est.std','ci.lower','ci.upper','pvalue')] # same but standardized
      
      estimates <- rbind(estimates, cbind(rep(f, nrow(estimates_un)), estimates_un, standardized)) 
    } 
  }
  save(fit_measures, estimates, failed, correrations, data_summary, times,
       file = paste0(out_folder,'/',substr(param[[1]],1,4),'_',param[[3]],'.RData'))
  
  # return(fits)
}

# ==============================================================================
# -------------------------------- Analysis ------------------------------------
# ==============================================================================

models <- list(
  # Self-reported depression ===================================================
  # - total fat mass/ FMI : 10/11 – 12/13 – 14 – 15½ / 16½ - 18 – 24
  list('sDEP_score', c(10.6, 12.8, 13.8, 16.6, 17.8, 23.8),
       'FMI',        c( 9.8, 11.8, 13.8, 15.4, 17.8, 24.5)),

  # - total lean mass/ LMI : 10/11 – 12/13 – 14 – 15½ / 16½ - 18 – 24
  list('sDEP_score', c(10.6, 12.8, 13.8, 16.6, 17.8, 23.8),
       'LMI',        c( 9.8, 11.8, 13.8, 15.4, 17.8, 24.5)),
  
  # – BMI: 10 – 13 – 14 – 16 - 18 - 24 # (with mean BMI 15.5 – 18)?
  list('sDEP_score', c(10.6, 12.8, 13.8, 16.6, 17.8, 23.8),
       'BMI',        c(10.7, 12.8, 13.8, 16,   17.8, 24.5)), # there is 17 too

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
  
  # – total lean mass/ LMI : 10 – 12 – 13.5 – 17.5 (or with mean lm 15.5 – 18)
  list('mDEP_score', c(9.6, 11.7, 13.1, 16.7),
       'LMI',        c(9.8, 11.8, 13.8, 17.8)),
  
  # – BMI: 10 – 12 – 13 – 16  # (with mean BMI 15.5 – 18) 
  list('mDEP_score', c(9.6, 11.7, 13.1, 16.7),
       'BMI',        c(9.8, 11.8, 12.8, 16)), # there is 17 too but lower correlation (0.2 vs. 0.4)

  list('mDEP_score', c(9.6, 11.7, 13.1, 16.7),
    'total_fatmass', c(9.8, 11.8, 13.8, 17.8)),

  list('mDEP_score', c(9.6, 11.7, 13.1, 16.7),
   'total_leanmass', c(9.8, 11.8, 13.8, 17.8)),
  
  # – waist circumference: 10 – 12 – 13 – 16
  list('mDEP_score', c(9.6, 11.7, 13.1, 16.7),
       'waist_circ', c(9.8, 11.8, 12.8, 15.4))
)

# parallel::detectCores()


# Run each depression / CMR marker combination in parallel =====================
foreach(i=1:length(models)) %dopar% {
  param = models[[i]]
  run_all_models(param[[1]], param[[2]], param[[3]], param[[4]])
}
