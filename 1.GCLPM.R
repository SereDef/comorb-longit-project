# ==============================================================================
# ============= 1. GENERALIZED CROSS-LAG PANEL MODEL (GCLPM) ===================
# ==============================================================================

# args = commandArgs(trailingOnly = TRUE) # provide access to a copy of the command line arguments supplied
# 
# if (length(args) == 0) {
#   stop("Supply output folder name!")
# } else {
#   out_folder <- args[1] # e.g. 'mod1_lfree_pstat' # lambda free, parameters stationary
# }

# Set parameters
out_folder = '../results/mod1_lfree_pstat'
dir.create(out_folder)

fix_lambdas = !(grepl('lfree', out_folder)) # heterogeneity allowed to vary over time 
stationarity = !(grepl('pfree', out_folder)) # model long term effects stationary over time 

# Load dependencies
invisible(lapply(c('lavaan','tidySEM','foreach'), require, character.only = TRUE));
# Note: I also tried parallel and pbapply for parallel processing but foreach worked best

# Read in data
data <- readRDS('../mats/raw_data.rds')

# ==============================================================================
# -------------------------- Set-up and functions ------------------------------
# ==============================================================================

# All possible combinations of the 8 main parameters
m = t(expand.grid(lapply(numeric(8), function(x) c(1, 0)))) 
m = m[, colSums(m)>1] # remove some combinations (e.g. less than 2 parameters in the model)
# order by number of parameters included and rename 
mat = data.frame( m[,order(colSums(-m))], 
                  row.names = c('maCL_dep','maCL_cmr','maAR_dep','maAR_cmr',
                                'ltCL_dep','ltCL_cmr','ltAR_dep','ltAR_cmr'))
# generate names for model structure (indicating the parameters which were excluded)
nmat = mat 
for (r in 1:nrow(mat)){ # TODO: can't make it work with apply (can't get row.name)
  nmat[r, ] <- replace(nmat[r, ], which(nmat[r, ]==0), rownames(nmat[r, ]))
}
nms = gsub('1-','', gsub('-1', '', lapply(nmat, paste0, collapse='-'))); nms[1]='full'
names(mat) <- nms
rm(m, nmat, nms, r)

# Remove some comparisons, to speed up processing: e.g., always estimate long term AR terms
mat <- mat[, -grep('ltAR', names(mat))] # --> 64 models

write.csv(mat, '../mats/model_structure.csv')

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
clpm_formula <- function(var1='dep', var2='cmr', n_ocs=NULL, meas_time=list(c(), c()), 
                         mod='full', fix_lambdas=FALSE, stationarity=FALSE) {
  
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
    
  } else { temp_var1 = temp_var2 = rep(1, n_ocs) }
  
  if (is.null(n_ocs)) { n_ocs = length(temp_var1) }
  
  # Define model structure according to mat matrix
  for (r in row.names(mat)) {
    if (mat[mod][r,]==0) { assign(gsub('lt','', tolower(r)),'0') } else { assign(gsub('lt','', tolower(r)), gsub('lt','', r)) }
  }
  
  unit_effects = paste0('# Unit effects\n',
                        'eta_',var1,' =~ ', paste0('l',1:(n_ocs),'*',var1,1:(n_ocs), collapse=' + '), '\n',
                        'eta_',var2,' =~ ', paste0('l',1:(n_ocs),'*',var2,1:(n_ocs), collapse=' + '), '\n') 
  
  impulses = '# Impulses\n'
  for (i in 2:n_ocs) { impulses = paste0(impulses, 'u_',var1,i,' =~ ',var1,i,'\n',var1,i,' ~~ 0*',var1,i,'\n') }
  for (i in 2:n_ocs) { impulses = paste0(impulses, 'u_',var2,i,' =~ ',var2,i,'\n',var2,i,' ~~ 0*',var2,i,'\n') }
  
  regressions = '# Regressions\n'
  
  for (i in n_ocs:2) { 
    # assign different name to each parameter unless it was set to 0
    if (ar_dep!=0) { ar_depi = paste0(ar_dep,i-1) } else { ar_depi = ar_dep }
    if (cl_dep!=0) { cl_depi = paste0(cl_dep,i-1) } else { cl_depi = cl_dep }
    if (ar_cmr!=0) { ar_cmri = paste0(ar_cmr,i-1) } else { ar_cmri = ar_cmr }
    if (cl_cmr!=0) { cl_cmri = paste0(cl_cmr,i-1) } else { cl_cmri = cl_cmr }
    
    if (i > 2) { # No MA terms because no impulses at first occasion
      
      if (maar_dep!=0) { maar_depi = paste0(maar_dep,i-1) } else { maar_depi = maar_dep }
      if (macl_dep!=0) { macl_depi = paste0(macl_dep,i-1) } else { macl_depi = macl_dep }
      if (maar_cmr!=0) { maar_cmri = paste0(maar_cmr,i-1) } else { maar_cmri = maar_cmr }
      if (macl_cmr!=0) { macl_cmri = paste0(macl_cmr,i-1) } else { macl_cmri = macl_cmr }
      
      ma_terms_var1 = paste0(' + ',maar_depi,'*u_',var1,i-1,' + ',macl_depi,'*u_',var2,i-1) 
      ma_terms_var2 = paste0(' + ',maar_cmri,'*u_',var2,i-1,' + ',macl_cmri,'*u_',var1,i-1) 
      
    } else { ma_terms_var1 = ma_terms_var2 = '' }
    
    # Add depression regressions
    regressions = paste0(regressions, 
                         var1,i,' ~ ',ar_depi,'*',var1,i-1,' + ',cl_depi,'*',var2,i-1, ma_terms_var1,'\n',
                         var2,i,' ~ ',ar_cmri,'*',var2,i-1,' + ',cl_cmri,'*',var1,i-1, ma_terms_var2,'\n') 
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
  } # else { constraints = paste0(constraints, 'l',n_ocs,' == 1\n') } # only last lambda is set to 1
  
  ar_dep_con = ar_cmr_con = cl_dep_con = cl_cmr_con = ''
  
  if (stationarity) {
    rel='*'
    for (i in 1:(n_ocs-1)) {
      for (o in c(1:(n_ocs-1))[-i]) {
        
        ar_dep_con = paste0(ar_dep_con, 'AR_',var1,i,rel,temp_var1[i+1]-temp_var1[i],' == AR_',var1,o,rel,temp_var1[o+1]-temp_var1[o],'\n')
        ar_cmr_con = paste0(ar_cmr_con, 'AR_',var2,i,rel,temp_var2[i+1]-temp_var2[i],' == AR_',var2,o,rel,temp_var2[o+1]-temp_var2[o],'\n')
        
        if (cl_dep!=0) { cl_dep_con = paste0(cl_dep_con, 'CL_',var1,i,rel,temp_var1[i+1]-temp_var2[i],
                                             ' == CL_',var1,o,rel,temp_var1[o+1]-temp_var2[o],'\n') }
        if (cl_cmr!=0) { cl_cmr_con = paste0(cl_cmr_con, 'CL_',var2,i,rel,temp_var2[i+1]-temp_var1[i],
                                             ' == CL_',var2,o,rel,temp_var2[o+1]-temp_var1[o],'\n') }
      }
      
    }
    constraints = paste0(constraints, ar_dep_con, ar_cmr_con, cl_dep_con, cl_cmr_con)
  }
  f = paste0(unit_effects, impulses, regressions, comevement, restrictions, constraints)
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
run_sem <- function(mod, data, dep_temp, cmr_temp, dep_name, cmr_name, 
                    normalize=TRUE, plot_semgraph=FALSE) {
  cat('- ', which(names(mat)==mod), mod)
  
  if (normalize) { 
    minmax_norm <- function(x, ...) { return( (x - min(x,...)) / (max(x,...) - min(x,...)) )}
    data <- sapply(data, minmax_norm, na.rm = TRUE) 
  }
  
  # Run model
  set.seed(310896)
  m <- lavaan::sem(clpm_formula(meas_time=list(dep_temp, cmr_temp), mod=mod, 
                                stationarity=stationarity, fix_lambdas=fix_lambdas), 
                   data, missing='fiml', se ='robust',
                   verbose=FALSE)
  
  if (lavInspect(m, "converged")){
    if (plot_semgraph) {
      # Print out the layout 
      save_semgraph(m, n_ocs=ncol(data)/2, name=paste(substr(dep_name,1,4),cmr_name,mod, sep='_'))
    }
  } else {
    # Return warning message ( overwrites the model object )
    m <- paste(gsub('lavaan WARNING:', '', names(warnings())), collapse='\n') }
  
  return(m)
}

## Run the 64 SEM models in parallel :)
run_all_models <- function(dep_name, dep_times, cmr_name, cmr_times, which_models=names(mat), 
                           verbose=FALSE){
  # Make dataframe
  d <- make_df(dep_name, dep_times, cmr_name, cmr_times, verbose=verbose)
  if (verbose) { print(summary(d$data)) }
  
  # save variable names and summary for dashboard output 
  dep_summ <- do.call(cbind, lapply(data[,sel(dep_name, dep_times)], summary))
  cmr_summ <- do.call(cbind, lapply(data[,sel(cmr_name, cmr_times)], summary))
  dat_summ <- cbind(dep_summ, cmr_summ) 
  dat_summ <- rbind(dat_summ, 'N_obs'= nrow(data) - dat_summ["NA's",])
  
  cat('\n------------------------------------------------------\n',cmr_name,' ~ ',dep_name)
  cat('\n------------------------------------------------------\nFitting the models...')
  
  # stat_seq <- c(T, rep(F, 2)) # ncol(mat)-1)) # Stationary is true in the full models but false in all others
  
  # Run this bitch (in parallel)
  start <- Sys.time(); cat(' started at: ', as.character(start), '\n') 
  fits <- mapply(run_sem, mod=which_models, # parallel::mc
                 MoreArgs=list(data=d$data, dep_temp=dep_times, cmr_temp=cmr_times, dep_name=dep_name, cmr_name=cmr_name), 
                 SIMPLIFY=FALSE) # mc.cores=12)
  end <- Sys.time()
  
  cat('\nDone! Runtime: ',round(end - start, 2),' hours\n------------------------------------------------------\n')
  
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
  save(dat_summ, fit_meas, estimates, failed, file = paste0(out_folder,'/',substr(dep_name,1,4),'_',cmr_name,'.RData'))
  
  # return(fits)
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

# parallel::detectCores()


# Run each depression / CMR marker combination in parallel =====================
foreach(i=1:length(models)) %dopar% {
  param = models[[i]]
  run_all_models(param[[1]], param[[2]], param[[3]], param[[4]])
}
