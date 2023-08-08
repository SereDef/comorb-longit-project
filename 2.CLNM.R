library(psychonetrics)
library(dplyr)
library(tidyr)
library(psych)
library(gdata)


data <- readRDS('../Mats/raw_data.rds')

which(sapply(data, is.character)) 

sel <- function(var, times=NULL){
  subs = names(data)[grep(paste(var, collapse='|'), names(data))]
  if (!is.null(times)) { subs = subs[grep(paste(paste0('_',times), collapse='|'), subs)] }
  return(subs)
}
# sel(c(paste0('sDEP',0:2)), 23.8)
# ==============================================================================
library(mgm) # Estimating Time-Varying Mixed Graphical Models in High-Dimensional Data

d <- data[,sel(c('9.6y','9.8y'))]
d = d[, -which(names(d)%in%sel(c('score','age')))]

dc <- d[complete.cases(d), ] # doesnt take missign values 

vtype <- c(rep('c',13), rep('g',18))
vlevl <- c(rep( 3, 13), rep( 1 ,18))

set.seed(1)
fit_mgm <- mgm::mgm(data = dc,
                   type = vtype,
                   levels = vlevl,
                   k = 2,
                   lambdaSel = "CV")

qgraph::qgraph(fit_mgm$pairwise$wadj,
       edge.color = fit_mgm$pairwise$edgecolor,
       layout = 'spring',
       nodeNames = names(d),
       legend = TRUE)
# fit_mgm$pairwise # = weighted adjacency matrix and the signs (if defined) of the parameters
# fit_mgm$interactions # list of all recovered interactions (cliques) and a list of parameters associated with all cliques; 
# fit_mgm$intercepts # all estimated thresholds/intercepts
# fit_mgm$nodemodels # list with the p glmnet objects from which all above results are computed. 
# We inspect the weigthed adjacency matrix stored in fit_mgm$pairwise$wadj

# ==============================================================================
library(qgraph)
library(lavaan)

d <- data[,sel(c('9.6y','9.8y'))]
d = d[, -which(names(d)%in%sel(c('score','age')))]

# cor_auto detects ordinal variables (with up to 7 unique integer values) and uses
# lavaan to estimate polychoric, polyserial and Pearson correlations.
mat <- qgraph::cor_auto(data[,sel('DEP')], 
                        forcePD=T # search for nearest positive-definite correlation matrix
                        # note: this can lead to unstable results
                        )

gra <- qgraph::qgraph(mat, 
                      graph = 'glasso', # with EPIC model selection
                      sampleSize = nrow(data),
                      layout = 'spring', # node placement
                      tuning = 0, # EBIC hyperparameter (gamma): [0-0.5] higher=fewer connections
                      lambda.min.ratio = 0.05 # minimal lambda ratio used in EBICglasso, defaults to 0.01.
                      )
# For more control on the estimation procedure, use EBICglasso() [= qgraph(..., graph="glasso")]
getWmat(gra) # estimated weights matrix (can be obtained directly using EBICglasso)

centralityPlot(gra) # centrality indices

d = data[,sel('17.8y')]
d = d[, -which(names(d)%in%sel(c('height','weight','tot_chol','score','age'), '17.8'))]

# =====================

library(bootnet)

res = estimateNetwork(d, default = 'EBICglasso', threshold=T)
plot(res, cut=0, theme='colorblind', layout='spring', labels=names)

names = c('felt\nmiserable\nor unhappy', 
          'had\nfun',
          'did not\nenjoy\nanything', 
          'so tired\njust sat\naround', 
          'was very\nrestless',
          'felt they\nwere no good\nanymore', 
          'cried\na lot',
          'felt\nhappy', 
          'hard to\nthink or\nconcentrate', 
          'hated\nthemselves',
          'enjoyed doing\nlots of things',
          'felt they\nwere a\nbad person', 
          'felt\nlonely', 
          'nobody\nreally loved\nthem', 
          'never\nas good as\nothers', 
          'felt did\neverything\nwrong', 
          'had a\ngood\ntime',
          
          'total\nfat', 'total\nlean','trunk\nfat','android\nfat', 
          'liver\nfat',
          'SBP', 'DBP',
          'heart\nrate', 
          'cf-PWV', 
          'ejection\nduration',
          'IMT',
          'intra-ventr.\nthickness\n(dyastolic)', 
          'intra-ventr.\nthickness\n(systolic)',
          'left\natrial\nsize',
          'left ventr.\nvolume\n(dyastolic)', 
          'left ventr.\nvolume\n(systolic)', 
          'poster. wall\nthickness\n(dyastolic)', 
          'poster. wall\nthickness\n(systolic)', 
          'HDL-c', 'LDL-c',
          'insulin',
          'triglycerides',
          'glucose',  
          'CRP',
          'BMI',  
          'FMI')
# ==============================================================================

# TODO: figure out liver measure 
# Liver Scan: Fat in liver (0-1) vs. Fatty liver result (CAP value; dB/m) 100-400 corr=.15

check <- function(var, times=c('10.6y','17.8y','23.8y')) {
  vars = paste0(var,'_',times)
  for (v in vars) { print(summary(as.factor(data[,v]))); cat('\n') }
  return(data[,vars])
}

felt_miserable_or_unhappy <- check('sDEP01')
did_not_enjoy_anything    <- check('sDEP03') # didn't enjoy anything at all
so_tired_just_sat_around  <- check('sDEP04') # felt so tired that they just sat around and did nothing
was_very_restless         <- check('sDEP05') 
they_were_no_good_anymore <- check('sDEP06') # felt they were no good any more
cried_a_lot               <- check('sDEP07') 
felt_happy                <- check('sDEP08')
hard_to_think_concentrate <- check('sDEP09') # found it hard to think properly or concentrate
hated_themselves          <- check('sDEP10')
felt_were_a_bad_person    <- check('sDEP12')  # felt they were a bad person
felt_lonely               <- check('sDEP13')
nobody_really_loved_them  <- check('sDEP14') # thought nobody really loved them 
never_as_good_as_others   <- check('sDEP15') # thought they could never be as good as other kids/people
felt_did_everything_wrong <- check('sDEP16') # felt they did everything wrong

# Other positive items are not matching between timepoints
# has_been_having_fun      <- check(c('sDEP02_17.8y'))
# enjoyed_lots_of_things   <- check(c('sDEP11_17.8y')) # enjoyed doing lots of things
# has_had_a_good_time      <- check(c('sDEP17_17.8y'))
# laughed_a_lot            <- check(c('sDEP19_23.8y')) 
# looked_forward_to_future <- check(c('sDEP20_23.8y')) # looked forward to the day ahead
# positive_about_future.   <- check(c('sDEP21_23.8y')) # felt really positive about the future 
# felt_valued              <- check(c('sDEP22_23.8y'))

# height <- data[,sel('height', c(17.8,24.5))]
# weight <- data[,sel('weight', c(17.8,24.5))]
check2 <- function(var, times=c(9.8, 17.8, 24.5)) {
  d = data[,sel(var, times)]
  print(summary(d))
  print(round(cor(d, use='pairwise.complete.obs'),2))
  return(d)
}
bmi        <- check2('BMI')
waist_circ <- check2('waist_circ', c(9.8, 15.4, 24.5))
totl_fatm  <- check2('total_fatmass')
totl_leanm <- check2('total_leanmass')
# andr_fatm <- check2('android_fatmass')
# trun_fatm <- check2('trunk_fatmass')
# liver_fat <- check2('liver_fat')
sbp       <- check2('SBP', c(10.6, 17.8, 24.5))
dbp       <- check2('DBP', c(10.6, 17.8, 24.5))
# imt       <- check2('IMT')
cf_pwv    <- check2('PWV', c(10.6, 17.8, 24.5))
# hrt_rate  <- check2('heart_rate')
# eject_dur <- check2('eject_dur')
# heart echo
# iv_thick_s  <- check2('sIVS') # interventricular septum thickness in systole  (cm)
# iv_thick_d  <- check2('dIVS') #               ' " "                  diastole (cm)
# pw_thick_s  <- check2('sPWT') # posterior wall thickness in systole  (cm)
# pw_thick_d  <- check2('dPWT') #               " " "         diastole  (cm)
# lv_volum_s  <- check2('sLVID') # internal dimension (volume?) end-systolic  (cm)
# lv_volum_d  <- check2('dLVID') #               " " "          end-diastolic (cm)
# metabolites
# tot_chol    <- check2('tot_chol')
hdl_chol    <- check2('HDL_chol')
ldl_chol    <- check2('LDL_chol')
insulin     <- check2('insulin')
triglyc     <- check2('triglyc')
# glucose     <- check2('glucose') # highly correlated but less missing 
crp         <- check2('CRP')

# only relevant variable for design matrix
d = cbind(cried_a_lot, did_not_enjoy_anything, felt_did_everything_wrong, felt_happy, 
          felt_lonely, felt_miserable_or_unhappy, felt_were_a_bad_person, hated_themselves,
          hard_to_think_concentrate, never_as_good_as_others, nobody_really_loved_them, 
          so_tired_just_sat_around, they_were_no_good_anymore, was_very_restless,
          bmi, waist_circ, 
          totl_fatm, totl_leanm, 
          sbp, dbp, cf_pwv, # hrt_rate, imt,
          # eject_dur,lv_volum_d, lv_volum_s, pw_thick_d, pw_thick_s, iv_thick_d, iv_thick_s, 
          hdl_chol, ldl_chol, insulin, triglyc, crp) # glucose, 

names = c('cried\na lot', 'did not\nenjoy\nanything', 'felt did\neverything\nwrong', 'felt\nhappy', 
          'felt\nlonely', 'felt\nmiserable\nor unhappy', 'felt they\nwere a\nbad person', 'hated\nthemselves',
          'hard to\nthink or\nconcentrate', 'never\nas good as\nothers', 'nobody\nreally loved\nthem', 
          'so tired\njust sat\naround', 'felt they\nwere no good\nanymore', 'was very\nrestless',
          'BMI', 'waist\ncircumference',
          'total\nfat\nmass', 'total\nlean\nmass', 
          'SBP', 'DBP', 'cf-PWV', # 'heart\nrate', 'IMT',
          # 'ejection\nduration','left ventr.\nvolume\n(dyastolic)', 'left ventr.\nvolume\n(systolic)', 
          # 'poster. wall\nthickness\n(dyastolic)', 'poster. wall\nthickness\n(systolic)', 
          # 'intra-ventr.\nthickness\n(dyastolic)', 'intra-ventr.\nthickness\n(systolic)', 
          'HDL-c', 'LDL-c','insulin', 'triglycerides', 'CRP') # 'glucose', 

# TODO: add sex variable 
sex <- rbinom(nrow(data), 1, 0.5)
# random binomial deviate generator function creating a vector of length nrow(data) 
# containing '0' and '1' with success probability of 0.5
df <- cbind(data$sex, d)

# ---------------------------------------------------------------------------
nT = 3 # number of timepoints 
nS = ncol(d) # number of symptoms/markers
# Design matrix with nT columns and nS rows
des <- matrix(as.vector(colnames(d)), nrow=nS/nT, byrow = T) 

# Unpruned model, get fit
# panelgvar() specifies a graphical vector-autoregression (GVAR) model on panel data
# When using only observed variables in the network, the panel data model takes the form of a random intercept cross-lagged panel model, 
# except the contemporaneous and between-subjects structures are modeled as networks. 
# This allows for the fixed-effects decomposition into temporal, contemporaneous, and between-subjects networks. 
# NOTE: when full-information maximum likelihood (FIML) estimation is used, the panel data model is a multi-level GVAR model 
# (as implemented in mlVAR) with only random intercepts (no random network parameters). 
start_time <- Sys.time()
umod <- psychonetrics::panelgvar(d, vars = des, 
                                  #lambda=
                                 estimator = 'FIML', 
                                 missing ='pairwise', 
                                 verbose = T) %>% 
  psychonetrics::runmodel() # run unpruned 
# Warning messages:
# 1: In psychonetrics::runmodel(.) :
#   Information matrix or implied variance-covariance matrix was not positive semi-definite. This can happen because the model is not identified, or because the optimizer encountered problems. Try running the model with a different optimizer using setoptimizer(...).
# 2: In psychonetrics::runmodel(.) :
#   One or more parameters were estimated to be near its bounds. This may be indicative of, for example, a Heywood case, but also of an optimization problem. Interpret results and fit with great care. For unconstrained estimation, set bounded = FALSE.
# 3: In psychonetrics::runmodel(.) :
#   Model might not have converged properly: mean(abs(gradient)) > 1.
Sys.time() - start_time
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

namedMatrix <- function(model, vrbls, type){
  # Type gets values = 'PDC' (namedTemporal), 
  # 'omega_zeta_within' (namedContemporaneous),
  # 'omega_zeta_between' (namedBetween)
  mat <- psychonetrics::getmatrix(model, type) # extract an estimated matrix from model
  rownames(mat) <- colnames(mat) <- vrbls
  return(mat)
}

prunedtemporal <- namedMatrix(pruned, names, 'PDC')
prunedcontemporaneous <- namedMatrix(pruned, names, 'omega_zeta_within')
prunedbetween <- namedMatrix(pruned, names, 'omega_zeta_between')

layout <- qgraph::averageLayout(prunedtemporal, prunedcontemporaneous, prunedbetween, layout = "spring")
groups <- c(rep('dep',14), rep('cmr',12))

library(qgraph)
qgraph(prunedtemporal, labels=names, threshold=0.01,
       title='Temporal network',
       theme = "colorblind", 
       repulsion = 0.7, 
       maximum = 1, 
       layout= layout, 
       groups = groups, 
       color = c('#FFEAF5','#E1EDFF'),
       vsize = 8, 
       label.scale.equal=T,
       vsize2=7,
       shape ='ellipse',
       # borders=F,
       border.width=0.5,
       label.cex = 1.1,
       legend = FALSE,
       filetype= "pdf",
       filename= "tempPlot")

qgraph(prunedcontemporaneous, labels=names, threshold=0.01,
       title='Contemporaneous network',
       theme = "colorblind", 
       repulsion = 0.7, 
       maximum = 1, 
       layout= layout, 
       groups = groups, 
       color = c('#FFEAF5','#E1EDFF'),
       vsize = 8, 
       label.scale.equal=T,
       vsize2=7,
       shape ='ellipse',
       # borders=F,
       border.width=0.5,
       label.cex = 1.1,
       legend = FALSE,
       filetype= "pdf",
       filename= "contempPlot")

qgraph(prunedbetween, labels=names, threshold=0.01,
       title='Between-person network',
       theme = "colorblind", 
       repulsion = 0.7, 
       maximum = 1, 
       layout= layout, 
       groups = groups, 
       color = c('#FFEAF5','#E1EDFF'),
       vsize = 8, 
       label.scale.equal=T,
       vsize2=7,
       shape ='ellipse',
       # borders=F,
       border.width=0.5,
       label.cex = 1.1,
       legend = FALSE,
       filetype= "pdf",
       filename= "betweenPlot")


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
genderComp <- df %>% tidyr::drop_na(sex)

Pconfigural <-  psychonetrics::panelgvar(genderComp, vars = des, estimator = 'FIML', missing = 'pairwise', 
                                         groups= 'sex',
                                         beta = prune.adjacencyT, 
                                         omega_zeta_within = prune.adjacencyC, 
                                         omega_zeta_between = prune.adjacencyB) %>% 
  psychonetrics::runmodel()
prunedconfig.fit <- Pconfigural@fitmeasures

Pconstrained <- Pconfigural %>% psychonetrics::groupequal('beta') %>% 
  psychonetrics::groupequal('omega_zeta_within') %>% 
  psychonetrics::groupequal("omega_zeta_between") %>% 
  psychonetrics::runmodel()
prunedconstrained.fit <- Pconstrained@fitmeasures

prunedcomparison <- psychonetrics::compare(Pconfigural = Pconfigural, 
                                           Pconstrained = Pconstrained)

aicDif <- Pconfigural@fitmeasures$aic.ll-Pconstrained@fitmeasures$aic.ll
bicDif <- Pconfigural@fitmeasures$bic-Pconstrained@fitmeasures$bic
prunedICs <- cbind(aicDif, bicDif)

#stepup model, get fit, parameters and matrices
stepup <- pruned %>% psychonetrics::stepup()
stepupfit <- stepup@fitmeasures
stepupparameters <- stepup %>% parameters()

stepuptemporal <- namedTemporal(stepup, names)
stepupcontemporaneous <- namedContemporaneous(stepup, names)
stepupbetween <- namedBetween(stepup, names)
stepuptemp.desc <- temporalDescript(stepuptemporal )
stepupcont.desc <- crossNetDescript(stepupcontemporaneous )
stepupbet.desc <- crossNetDescript(stepupbetween )

#stepup centrality
PDCcents<- centrality_auto(stepuptemporal )
PDCcentS <- PDCcents$node.centrality
contCents<- centrality_auto(getmatrix(stepup, "omega_zeta_within"))
contCentS <- contCents$node.centrality
betCents<- centrality_auto(getmatrix(stepup, "omega_zeta_between"))
betCentS <- betCents$node.centrality

#stepup adjacency matrices
stepupadjacencyT <- 1*(getmatrix(stepup, "beta")!=0)
stepupadjacencyB<- 1*(getmatrix(stepup, "omega_zeta_between")!=0)
stepupadjacencyC<- 1*(getmatrix(stepup, "omega_zeta_within")!=0)

# pruned invariance
Sconfigural<- panelgvar(genderComp, vars= des, estimator= "FIML", missing= "pairwise", groups= "Gender",
                        beta= stepupadjacencyT, omega_zeta_within= stepupadjacencyC, omega_zeta_between= stepupadjacencyB)%>% runmodel
stepupconfig.fit <- Sconfigural@fitmeasures

Sconstrained<- Sconfigural%>% groupequal("beta")%>% groupequal("omega_zeta_within")%>% groupequal("omega_zeta_between")%>% runmodel()
stepupconstrained.fit <- Sconstrained@fitmeasures

stepupcomparison<- psychonetrics::compare(Sconfigural= Sconfigural, Sconstrained= Sconstrained)

SaicDif<- Sconfigural@fitmeasures$aic.ll-Sconstrained@fitmeasures$aic.ll
SbicDif<- Sconfigural@fitmeasures$bic-Sconstrained@fitmeasures$bic
stepupICs <- cbind(SaicDif, SbicDif)

# Mother reported depression and fat mass 10-13 years
m_totfat3 <- data.frame('dep1'= data[,'mDEP_score_9.6y'], 'fat1'= data[,'total_fatmass_9.8y'],
                        'dep2'= data[,'mDEP_score_11.7y'],'fat2'= data[,'total_fatmass_11.8y'],
                        'dep3'= data[,'mDEP_score_13.1y'],'fat3'= data[,'total_fatmass_12.8y'])# sDEP_score_12.8y
# Mother reported depression and fat mass 10-17 years (including average fat between 15.4 and 17.8 measures)
m_totfat4 <- cbind(m_totfat3, data.frame('dep4'= data[,'mDEP_score_16.7y']))
m_totfat4$fat4 <- ifelse(rowSums(is.na(data[,sel('total_fatmass',c(15,17))]))<2, 
                                 rowMeans(data[,sel('total_fatmass',c(15,17))], na.rm=T), NA)

# Self reported depression and fat mass 12-24 years
m_totfat <- data[,c('sDEP_score_12.8y','total_fatmass_11.8y', # mDEP_score_11.7y
                    'sDEP_score_17.8y','total_fatmass_17.8y',
                    'sDEP_score_23.8y','total_fatmass_24.5y')]

names(m_totfat) <- c('dep1','fat1','dep2','fat2','dep3','fat3')


samp = m_totfat[rowSums(is.na(m_totfat)) != ncol(m_totfat), ]

# mf <- cov(samp, use='pairwise.complete.obs')

# ==============================================================================
# ==============================================================================
sel <- function(var, times=NULL, d=data){
  subs = names(d)[grep(paste(var, collapse='|'), names(d))]
  if (!is.null(times)) { subs = subs[grep(paste(paste0('_',times), collapse='|'), subs)] }
  return(subs)
}

library(bootnet)

names = data.frame('DEP01'='felt\nmiserable\nor unhappy', 
                   'DEP02'='had fun', # only few
                   'DEP03'='did not\nenjoy\nanything', 
                   'DEP04'='so tired\njust sat\naround', 
                   'DEP05'='was very\nrestless',
                   'DEP06'='felt they\nwere no good\nanymore', 
                   'DEP07'='cried\na lot', 
                   'DEP08'='felt\nhappy', 
                   'DEP09'='hard to\nthink or\nconcentrate', 
                   'DEP10'='hated\nthemselves',
                   'DEP11'='enjoyed\ndoing lots of\nthings', # only few
                   'DEP12'='felt they\nwere a\nbad person',
                   'DEP13'='felt\nlonely', 
                   'DEP14'='nobody\nreally loved\nthem', 
                   'DEP15'='never\nas good as\nothers', 
                   'DEP16'='felt did\neverything\nwrong', 
                   'DEP17'='had a\ngood time', 
                   #'DEP18'='had a\ngood time', 
                   #'DEP19'='had a\ngood time', 
                   #'DEP20'='had a\ngood time', 
                   #'DEP21'='had a\ngood time', 
                   #'DEP22'='had a\ngood time', 
                   'BMI'='BMI',
                   'waist_circ'='waist\ncircumference',
                   'hip_circ'='hip\ncircumference',
                   'total_fatmass'='total\nfat', 
                   'android_fatmass'='android\nfat', 
                   'trunk_fatmass'='trunk\nfat', 
                   'liver_fat'='liver\nfat',
                   'SBP'='SBP', 
                   'DBP'='DBP', 
                   'IMT'='IMT', 
                   'cfPWV'='cf-PWV', 
                   'heart_rate'='heart\nrate', 
                   'eject_dur'='ejection\nduration',
                   'dLVID'='left ventr.\ndiameter\n(dyastolic)', 
                   'dLVID'='left ventr.\ndiameter\n(systolic)', 
                   'dPWT'='poster. wall\nthickness\n(dyastolic)', 
                   'sPWT'='poster. wall\nthickness\n(systolic)', 
                   'dIVS'='intra-ventr.\nseptum\n(dyastolic)', 
                   'sIVS'='intra-ventr.\nseptum\n(systolic)', 
                   'insulin'='insulin', 
                   'glucose2'='glucose', 
                   'HDL_chol'='HDL-c', 
                   'LDL_chol'='LDL-c', 
                   'triglyc'='triglycerides', 
                   'CRP'='CRP',
                   'IL_6'='IL-6')

subdata = data[, sel('_9.')]
names(subdata)
# rm agez, total scores, height and weight
subdata = subdata[, sel(names(names), d=subdata)]
names(subdata)

labels = list()
for (v in names(subdata)) { 
  print(v)
  for (sv in names(names)) {
    print(sv)
    if (grep(sv, v)==T) { print( c(v, sv)) 
      break() }
  }
}
unlist(lapply(names(names), function(x) grep(x, names(subdata))))
labels
# subdata = data[, sel('_9.')]
# subdata = data[, sel('_10.')]
# subdata = data[, sel('_11.')]

subdata = subdata[complete.cases(subdata),]

res = bootnet::estimateNetwork(subdata, default = 'EBICglasso', threshold=T)

plot(res, cut=0, labels=labels, theme='colorblind', layout='spring')

library("qgraph")

centralityPlot(res, orderBy='Strength')


