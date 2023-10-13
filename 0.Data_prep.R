require(foreign)

# File name and location 
filename <- 'EarlyCause_AHupdated_CIDB2957_21OCT21.sav' # change if needed
datapath <- dirname(file.choose()) # volumes[...] (EarlyCause/WP3/SD/Data)

# Read in the file
full <- foreign::read.spss(file.path(datapath, filename), 
                           use.value.labels=F, to.data.frame=T)
names(full) <- tolower(names(full)) # all column names to lower case

# Initiate dataframe and set up unique child ID
data <- data.frame('IDC' = paste(full$cidb2957, # mother ID
                                 gsub('\\s+', '', full$qlet), # sibling ID (A or B)
                                 sep = "_"), # --> format "1_A"
                   'sex' = as.factor(full$kz021)) 

# ==============================================================================
# -------------------------------- Functions -----------------------------------

# Find variables in the dataset
f <- function(pref) { print(sort(names(full)[grep(pref, names(full))]))}

# Extract and clean depression items
extr_dep <- function(pref, nrs, remove=c(), age, plus=NULL, reporter='s'){
  # remove positively worded items (not in SMFQ)
  if (length(remove) > 0) { nrs = nrs[-remove] }
 
  n_adj <- sprintf('%02d', nrs) # add leading 0s to the number sequence
  if(!is.null(plus)) { n_adj <- paste0(n_adj, plus) } # if necessary, add one element at the end
  
  # identify variables 
  dep <- full[, paste0(pref, n_adj)]
  message(ncol(dep), ' items.\nScore:')
  # any other variables left ? 
  # print(setdiff(names(full)[grep(pref, names(full))], names(dep)))
  
  # rename items 
  i_adj <- sprintf('%02d', 1:13)
  names(dep) <- paste0(reporter,'DEP', i_adj)

  # Transform response to 0=not true; 1=sometimes true; 2=true
  # note: 22.9y and 23.8y depression scores are reversed for some reason
  if (age %in% c('ypb9992','ypc2650')) { resc=1 } else { resc=3 }
  dep <- as.data.frame(sapply(dep, function(x) abs(x-resc) )) 
  
  # Check the direction is as expected
  if (any(colMeans(dep, na.rm=T) > 1)) { cat('ATTENTION: check direction of answer scoring!'); print(summary(dep)) }
  # Take care of reverse phrased items 
  # rev = which(names(dep) %in% paste0(reporter,'DEP',c('02','08','11',17:22)))
  # dep[,rev] <- as.data.frame(sapply(dep[,rev], function(x) abs(x-2) ))
  
  # Compute total score 
  score_var <- paste0(reporter,'DEP_score')
  dep[,score_var] <- ifelse(rowSums(is.na(dep))!=ncol(dep), rowSums(dep, na.rm=T), NA)
  print(summary(dep[,score_var]))
  # Add age variable (in months) and display its descriptives
  age_var = paste0(reporter,'DEP_age')
  dep[,age_var]<- full[,age]/12
  message('\nAge:'); print(summary(dep[,age_var]))
  
  # Histogram 
  hist(dep[,score_var], main=paste(score_var,round(median(dep[,age_var], na.rm=T),1)), 
       breaks=max(dep[,score_var], na.rm=T), xlim=c(0, 27)); 
  
  # Add age mark to variable 
  names(dep) <- paste0(names(dep),'_',round(median(dep[,age_var], na.rm=T),1),'y')
  
  return(dep)
}

# Compute BMI given height in cm and weight in kg 
bmi_calc <- function(h, w, title='BMI', age=''){
  h = h/100 # transform to meters
  bmi <- w / h^2
  # print(summary(bmi))
  hist(bmi, main=paste(title, age), breaks=max(bmi, na.rm=T))
  return(bmi)
}

# Remove outliers
rm_outliers <- function(var, cutoff = 5, verbose=TRUE) {
  quants <- quantile(var, probs=c(.25, .75), na.rm=TRUE);
  iqr <- quants[2]-quants[1]
  too_high <- which(var > quants[2]+cutoff*iqr)
  too_low  <- which(var < quants[1]-cutoff*iqr)
  var[ c(too_high, too_low) ] <- NA
  if (verbose){
    # message(names(var))
    cat(length(too_low)+length(too_high), ' outliers removed (', 
        length(too_low), ' below and ', length(too_high), ' above the cutoff).\n', sep='')
  }
  return(var)
}

# Extract and clean cardio-metabolic items
extr_cmr <- function(df, age){
  # identify variables 
  cmr <- as.data.frame(full[, names(df)])
  # rename them
  names(cmr) <- df
  
  # Remove outliers
  except = which(names(cmr) %in% c('insulin','glucose','triglyc','CRP','IL_6','canabis'))
  if (identical(except, integer(0))) { 
    cmr <- as.data.frame(sapply(cmr, rm_outliers))
  } else { cmr[,-except] = as.data.frame(sapply(cmr[,-except], rm_outliers)) }
  
  # cmr[, except] = as.data.frame(sapply(cmr[, except], rm_outliers, cutoff=10))
  # browser()
  
  # Add age variable (in months) and display its descriptives
  cmr[,'CMR_age']<- full[,age]/12
  # print(summary(cmr$CMR_age))
  
  # Compute BMI if all info is available
  if (length(grep('height|weight', names(cmr))) > 1) {
    message('Calculating BMI')
    cmr[,'BMI'] <- bmi_calc(cmr$height, cmr$weight)
    cmr[,'BMI'] <- rm_outliers(cmr[,'BMI'])
  }
  # Compute FMI if all info is available
  if (length(grep('height|total_fatmass', names(cmr))) > 1) {
    message('calculating FMI, TFI and LMI')
    new = c('FMI','TFI','LMI'); old = c('total_fatmass','trunk_fatmass','total_leanmass')
    for (i in 1:3) { cmr[,new[i]] <- bmi_calc(cmr[,'height'], cmr[,old[i]]/1000, 
        title=new[i], age=round(median(cmr$CMR_age,na.rm=T),1)) }
    cmr[,new] <- rm_outliers(cmr[,new])
  }
  # Add age mark to variable 
  names(cmr) <- paste0(names(cmr),'_',round(median(cmr$CMR_age,na.rm=T),1),'y')
  
  print(summary(cmr))
  return(cmr)
}

# ---- SMFQ mother-reported (13 items) / self-reported (17-16) -----------------

# 1	    Felt miserable or unhappy in the past 2 weeks
# 2	    Didn't enjoy anything at all in the past 2 weeks
# 3   	Felt so tired they just sat around and did nothing in the past 2 weeks
# 4	    Was very restless in the past 2 weeks
# 5	    Felt they were no good any more in the past 2 weeks
# 6	    Cried a lot in the past 2 weeks
# 7	    Found it hard to think properly or concentrate in the past 2 weeks
# 8     Hated themselves in the past 2 weeks
# 9	    Felt they were a bad person in the past 2 weeks
# 10    Felt lonely in the past 2 weeks
# 11    Thought nobody really loved them in the past 2 weeks
# 12	  Thought they would never be as good as other people in the past 2 weeks
# 13	  Felt they did everything wrong in the past 2 weeks

# ---------------------- Removed from the scoring ------------------------------
# 2 Had fun; 8 Felt happy; 11 Enjoyed doing lots of things; 17 Had a good time; 18 Laughed a lot; 
# 19 Looked forward to the day ahead; 20 Felt really positive about the future; 21 Felt valued; 22 Felt unhappy

pdf('histogr_data.pdf') # Note: saved in working directory :) 
# ------------------------------------------------------------------------------
mdep_09.6y <- extr_dep('ku6', 60:72, reporter='m', age = 'ku991a') # tot scores ku673a / ku673b + ku647 (child argues with brothers and sisters)

cmr_09.6y <- extr_cmr(data.frame('pub203'='height'), # cm (parent-reported) # 'pub204'='weight'), # kg (parent-reported)
                      age = 'pub295')

# Combine two weight measure (from two DXA visits)
full[,'f9_weight'] <- rowMeans(full[,c('f9dx010','f9ms057')], na.rm=T)
# Calculate waist to hip ratio
full[,'f9_whr'] <- full[,'f9ms018'] / full[,'f9ms020'] # Hip circumference (cm)

cmr_09.8y <- extr_cmr(data.frame('f9_weight'='weight', # kg (DXA)
                                   'f9ms018'='waist_circ', # cm
                                    'f9_whr'='waist_hip_ratio', # cm
                                   'f9dx135'='total_fatmass', # grams (--> trasformed later)
                                   'f9dx136'='total_leanmass',# grams (--> trasformed later)
                                   'f9dx126'='trunk_fatmass', # grams (--> trasformed later)
                                   'chol_f9'='tot_chol', # mmol/l
                                    'hdl_f9'='HDL_chol', # mmol/l
                                    'ldl_f9'='LDL_chol', # mmol/l
                                'insulin_f9'='insulin',  # mU/L
                                   'trig_f9'='triglyc',  # mmol/l,
                                    'crp_f9'='CRP', # mg/L
                                    'il6_f9'='IL_6'), # Interleukin 6 (pg/ml)
                      age = 'f9003c') # f9006c REVISIT? 196 children came twice, mostly 1 month apart

# Add BMI, FMI, TMI and LMI
add = data.frame('BMI'='weight','FMI'='total_fatmass','TFI'='trunk_fatmass','LMI'='total_leanmass')
for (v in names(add)) {
  if (v=='BMI') {div=1} else {div=1000} # weight is already in kg but fat/lean mass is in grams
  at='_9.8y'
  # Index to height (m) squared
  cmr_09.8y[,paste0(v,at)] <- bmi_calc(cmr_09.6y[,'height_9.6y'], cmr_09.8y[,paste0(add[v],at)]/div, title=v, age='~9.7')
  # Remove outliers
  cmr_09.8y[,paste0(v,at)] <- rm_outliers(cmr_09.8y[,paste0(v,at)], cutoff = 5)
}
# Correct insulin levels 
cmr_09.8y[,'insulin_9.8y'] <- rm_outliers(cmr_09.8y[,'insulin_9.8y'], cutoff = 10)

summary(cmr_09.8y)

# EXTRA ------------------------------------------------------------------------
# * Puberty
#   pub209 = frequency of participation in vigorous physical activity during past month
#   pub211 = how old was child when she had her first period (but few cases!)
#   pub230 = development stage of breasts
#   pub235 = development stage of pubic hair
#   pub250 = development stage of testes scrotum and penis
#   pub255 = development stage of pubic hair
# * Blood samples 
#   apoai_f9 = Apolipoprotein Al; apob_f9 = Apolipoprotein B mg/dl
#   hb_f9 = haemoblobin; hba1c_f9 = hemoglobin A1c (HbA1c), amount of blood sugar (glucose) attached to hemoglobin.
#   igf1_f9 = Insulin-like growth factor 1 (IGF-1)
# ------------------------------------------------------------------------------
dep_10.6y <- extr_dep('fddp1', 10:25, remove=c(2,8,11), age = 'fd003c') # Depression score = fddp130 (prorated?)

cmr_10.7y <- extr_cmr(data.frame('pub303'='height',  # cm (parent-reported)
                                 'pub304'='weight'), # kg (parent-reported)
                      age = 'pub397a')

cmr_10.6y <- extr_cmr(data.frame('fdms018'='waist_circ', # cm
                                 'fdar117'='SBP',  # Systolic blood pressure (mmHg)
                                 'fdar118'='DBP',  # Diastolic blood pressure (mmHg)
                                 'fdar114'='PWV'), # Pulse Wave Velocity (m/s) # radial - carotid!!
                      age = 'fd003c')
# EXTRA ------------------------------------------------------------------------
# + alcohol (fdaa492-93) + diet (fddd2-4..)

# MISSING ----------------------------------------------------------------------
#   fdms010 = Height (cm)
#   fdms026 = Weight (kg)
#   fdar111 = Brachial diameter in millimetres
#   fdar115 = Brachial distensibility coefficient
#   fdar116 = Brachial artery compliance
#   FDAA483 = Freq Child smoked cigarettes - FDAA484 = No. cigarettes Child smoked per week

# ------------------------------------------------------------------------------
mdep_11.7y <- extr_dep('kw60', 0:12, reporter='m', age = 'kw9991a') # tot scores kw6100a / kw6100b  + SDQ scores kw6600 - kw6605

# Calculate waist to hip ratio
full$fems_whr <- full[,'fems018'] / full[,'fems020'] # Hip circumference (cm)
cmr_11.7y <- extr_cmr(data.frame( 'pub403'='height', # cm (parent-reported) # Note: age slightly different ('pub497a')
                                 'fedx016'='weight', # kg --> 'pub404'='weight', # kg (parent-reported)
                                 'fems018'='waist_circ', # cm 
                                'fems_whr'='waist_hip_ratio', # cm 
                                 'fedx135'='total_fatmass', # grams (--> trasformed later)
                                 'fedx136'='total_leanmass',# grams (--> trasformed later)
                                 'fedx126'='trunk_fatmass'),# grams (--> trasformed later)
                       age = 'fe003c')

# MISSING 
#   fems010 = Height (cm)
#   fems026 = Weight (kg)
#   fems028a = Fat percentage, impedance

# ------------------------------------------------------------------------------
dep_12.8y <- extr_dep('ff65', 0:15, remove=c(2,8,11), age = 'ff0011a') # (+ date attendance b)

cmr_12.8y <- extr_cmr(data.frame('ff2000'='height', # cm
                                 'ff2030'='weight', # kg
                                 'ff2020'='waist_circ'), # cm
                      age = 'ff0011a')
# MISSING 
# ff2036 = Fat percentage (%)
# ff2620-27 = BP and pulse measures
# ------------------------------------------------------------------------------
mdep_13.1y <- extr_dep('ta50', 20:32, reporter='m', age = 'ta9991a') 
# + SDQ scores ta7025a - f 
# + Teenager upset or distressed about weight/body shape (ta6160)

cmr_13.1y <- extr_cmr(data.frame('pub503'='height',  # cm (parent-reported)
                                 'pub504'='weight'), # kg (parent-reported)
                      age = 'pub597a')

# ------------------------------------------------------------------------------
dep_13.8y <- extr_dep('fg72', 10:25, remove=c(2,8,11), age = 'fg0011a')

cmr_13.8y <- extr_cmr(data.frame('fg3100'='height',    # cm
                                 'fg3207'='weight',# kg (DXA)
                                 'fg3254'='total_fatmass',  # grams (--> trasformed later)
                                 'fg3255'='total_leanmass', # grams (--> trasformed later)
                                 'fg3245'='trunk_fatmass',  # grams (--> trasformed later)
                                 'fg3257'='android_fatmass',# grams (--> trasformed later)
                                 'fg1020'='heart_rate'), # at rest (BPM)
                      age = 'fg0011a')
# (+ fg4822-24,29: smoked cigarettes)
# (+ fg4873,77-80, fg4915,17,19: drinking)
# (+ fg5422-23,26,28-29 cannabis )
# (+ fg7363 = Factor IV Emotional Stability: Big 5 factor marker: Personality Questionnaire)

# MISSING
#   fg3120 = Waist circumference (cms)
#   fg3130 = Weight (Kgs)
#   fg3136 = Fat percentage (%)
#   fg3260 = Gynoid - fat mass (g)
#   fg1021-22 = BP at rest (V1 only)
#   gf6120-27 = BP and pulse measures

# ------------------------------------------------------------------------------
# DAWBA 
# fh6630 = DAWBA YP has been afraid of meeting new people, last 4 weeks
# fh6868 = DAWBA Whether any disorder present (any informant 6-band computer prediction, ICD-10 & DSM-IV)
# fh6869 = DAWBA Whether any externalising disorder present (ODD CD ADHD parent 6-band computer prediction, ICD-10 & DSM-IV)
# 70 = ADHD, 71 = Hyperkinesis, 72 = Any behavioural disorder (ODD or CD), 73 = ODD, 74 = CD, 
# 75 = Any emotional disorder (self-report), 76 = depression, 
# 77 = any anxiety, 78 = Generalised anxiety, 79 = Panic disorder, 80 = Agoraphobia, 81 = PTSD, 82 = Social Phobia, 83 = Specific Phobia, 
# 84 = any disorder present (any informant), 85 = externalising disorder present, 86 = ADHD, 88 = ODD or CD, 89 = ODD, 90 = CD,
# 91 = any emotional disorder, 92 = depressive, 93 = anxiety, 97 = PTSD

cmr_15.4y <- extr_cmr(data.frame('fh3000'='height', # cm
                              'fh3010'='weight', # kg
                            # 'fh2209'='weight_dxa', # more missing values
                              'fh4020'='waist_circ', # cm [fh4021 = code]
                              'fh2254'='total_fatmass',  # grams (--> trasformed later)
                              'fh2255'='total_leanmass', # grams (--> trasformed later)
                              'fh2245'='trunk_fatmass',  # grams (--> trasformed later)
                              'fh2257'='android_fatmass',# grams (--> trasformed later)
                           #  'fh2006'='heart_rate', # maximum # only 4 values: 203,204,205,206
                            'chol_tf3'='tot_chol', # mmol/l
                             'hdl_tf3'='HDL_chol', # mmol/l
                             'ldl_tf3'='LDL_chol', # mmol/l
                         'insulin_tf3'='insulin',  # mU/L
                            'trig_tf3'='triglyc',  # mmol/l,
                         'glucose_tf3'='glucose',  # mmol/l
                           # 'glc_tf3'='glucose1', # more missing values
                             'crp_tf3'='CRP'), 
                   age = 'fh0011a') # fh5307 = Age in months of YP at clinic visit: TF3, but similar to fh0011a?
# (fh5100-10, actigraphy data )
# (fh5333 Number of times YP usually wakes up at night:)
# (fh8230-35 arguing with parents)
# (fh8430-32,40-41,50-51,55-56 smoking)
# (fh8512-13,42,45,65,71,79,81,87 drinking)
# (fh8610,12-13,16,21,31,61 cannabis

# apoa1_tf3 = Apolipoprotein AI; apob_tf3 = Apolipoprotein B mg/dl
# dha_tf3 = Docosahexaenoic acid (DHA) omega-3 fatty acid
# hb_tf3 = haemoblobin; testosterone_tf3 ?

# MISSING
# fh2260 = gynoid fat mass
# fh3016 = Fat percentage (%)
# fh2030-37 = BP and pulse measures

# ------------------------------------------------------------------------------
dep_16.6y <- extr_dep('ccs45', 0:15, remove=c(2,8,11), age = 'ccs9991a')
# (ccs2000-2221) life events 
# (ccs3500,10, 40-49 alcohol)
# (ccs4000,05,10,30,35,40 smoking)
# (ccs4060,65,70,75,4141 cannabis)
# (ccs4141,53,54,60-68 other drugs)
# ccs5010 asthma?

# MISSING 
# DEPRESSION SYMPTMONS ccs2660-76 maybe in psychosis quest
# eating disorders (ccs5500-5600)

mdep_16.7y <- extr_dep('tc40',30:42, reporter='m', age = 'tc9991a') #, + SDQ scores tc4025a - f 

cmr_16.0y <- extr_cmr(data.frame('pub803'='height',  # cm (parent-reported)
                                 'pub804'='weight'), # kg (parent-reported)
                      age = 'pub897a')
cmr_17.0y <- extr_cmr(data.frame('pub903'='height',  # cm (parent-reported)
                                 'pub904'='weight'), # kg (parent-reported)
                      age = 'pub997a')

# ------------------------------------------------------------------------------
dep_17.8y <- extr_dep('ccxd9', 0:15, remove=c(2,8,11), age = 'ccxd006') # (05 in years), CCXD917 total score

# Combine some three BP measurements
full$fjel_sbp <- rowMeans(full[,c('fjel036','fjel040')], na.rm=T) # excluding first measurement 'fjel032'
full$fjel_dbp <- rowMeans(full[,c('fjel037','fjel041')], na.rm=T) # excluding first measurement 'fjel033'

# Calculate cardiac measurements: 
# Left ventricular mass (LVM) (g) = 0.8⋅[1.04⋅((LVEDD + IVSd + PWd)^3 − LVEDD^3)] + 0.6, where:
# LVEDD – LV end-diastolic dimension, given in centimeters (cm);
# IVSd – Interventricular septal end-diastole, given in centimeters (cm);
# PWd – Posterior wall thickness at end-diastole, given in centimeters (cm); and
# 1.04 – Heart muscle density in g/cm³.
full$fjgr_LVM = 0.8 * ( 1.04 * ( rowSums(full[,c('fjgr043',  # Average left ventricular internal diameter during diastole (cm)
                                                'fjgr031',  # Average interventricular septum in diastole (cm)
                                                'fjgr053')] # Average posterior wall thickness in diastole (cm)
                                     )^3 - full[,'fjgr043']^3 ) ) + 0.6

# Relative wall thickness (RWT), defined as 2 times posterior wall thickness divided by the LV diastolic diameter
full$fjgr_RWT <- 2 * full[,'fjgr053'] / full[,'fjgr043']
# Ejection fraction (%) is estimated by fractional shortening (%), calculated using Teichholz's formula:
# FS = (LV end-diastolic dimension – LV end-systolic dimension) / LV end-diastolic dimension
# fjgr047 = Average left ventricular internal diameter during systole (cm)
full$fjgr_FS <- (full[,'fjgr043'] - full[,'fjgr047']) / full[,'fjgr043'] * 100

# cigarette smoking 
full$fjsm_tot <- full[,'fjsm500'] # weekly number of cigarettes 
full$fjsm_tot[is.na(full$fjsm_tot) & (full[,'fjsm350']==2 & full[,'fjsm450']==2)] <- 0 # smokes every day = no & smokes every week = no
full$fjsm_tot[is.na(full$fjsm_tot)] <- full[is.na(full$fjsm_tot), 'fjsm400']*7 # daily cigarettes present but weekly NA

cmr_17.8y <- extr_cmr(data.frame('fjmr020'='height', # cm
                              'fjmr022'='weight', # kg
                              'fjdx135'='total_fatmass',  # grams (--> transformed later)
                              'fjdx136'='total_leanmass', # grams (--> transformed later)
                              'fjdx126'='trunk_fatmass',  # grams (--> transformed later)
                              'fjdx138'='android_fatmass',# grams (--> transformed later)
                            # 'fjli100'='liver_fat', # only binary (0,1) with 43 cases
                                       # arteries
                             'fjel_sbp'='SBP', # Omron BP (seated after 5min rest)
                             'fjel_dbp'='DBP', # Omron BP (seated after 5min rest)
                             'fjar079d'='IMT', # Baseline (end diastole) intima-media thickness (mm)
                             'fjar083d'='PWV', # carotid to femoral (m/s) average
                           # 'fjar088d'='crPWV', # carotid to radial (m/s) average
                           # 'fjel115'='eject_dur', # Ejection Duration, ms (? --> from radial artery tonometry )
                                       # heart 
                              'fjel116'='heart_rate',
                             'fjgr_LVM'='LVM', # Left ventricular mass (LVM) (g)
                             'fjgr_RWT'='RWT', # Relative wall thickness (RWT) (cm)
                              'fjgr_FS'='FS',  # Fractional shortening (%)
                              #'fjgr039'='LAS',  # Average left atrial size, m-mode (cm)
                              #'fjgr043'='LVIDD',# Average left ventricular internal diameter during diastole (cm)
                             # assuming the age is the same here 
                             'chol_tf4'='tot_chol', # mmol/l
                              'hdl_tf4'='HDL_chol', # mmol/l
                              'ldl_tf4'='LDL_chol', # mmol/l
                          'insulin_tf4'='insulin',  # mU/L
                             'trig_tf4'='triglyc',  # mmol/l,
                          'glucose_tf4'='glucose',  # mmol/l
                            # 'glc_tf4'='glucose1', # more missing values
                              'crp_tf4'='CRP', 
                                     # other
                             'fjal4000'='alcohol', # Alcohol use disorder identification test (AUDIT) score
                             'fjsm_tot'='smoking', # Average number of cigarettes per week
                             'fjdr4500'='canabis'),# Cannabis Abuse Screening Test (CAST) score
                   age = 'fj003a') # (b in years)

# Also available, bu not included: ---------------------------------------------
# * {CIS-R} Depression 
#   fjci101-3 = Number of primary / secondary / somatic depressive symptoms 
#   fjci350 = Depression score; 603 = Mild, 608 = Moderate, 609 = severe depressive episode
#   fjci366 = Suicide risk score; fjci155 = Neurotic symptom score; fjci604 = Panic disorder syndrome
#   fjci1000 = Deptot:Sum of all the 5 depression symptom scores (depression, fatigue, concentration, sleep, depressive thoughts) + 1001-3 binary
# * Pulscor SBP/DBP (seated after 5min rest)
#   fjel022-23 --> 265 obs only (measured at the same time as Omron for comparison)
# Echocardiography: 
#   fjgr031 # Average interventricular septum in diastole (cm); fjgr035 # Average interventricular septum in systole (cm)
#   fjgr053 # Average posterior wall thickness in diastole (cm); fjgr057 # Average posterior wall thickness in systole (cm)
#   fjgr048 # Old ejection fraction [not collected systematically]: GRACE substudy
# * Medical history of vascular conditions
#   fjar042 = treatment for hypertension; 43-46 = diabetes...; 47 = high cholesterol
# * Alcohol use 
#   fjal4001 = AUDIT score categorized in risk levels defined by the Australian government
# * Smoking
#   fjsm050 = ever smoked a whole cigarette; 100 = age first smoked; 150 = n cigarettes lifetime; 250 = smoked in past 30 days; 300 = age last smoked;
#   fjsm1000-1001 = Fagerstrom Test for Nicotine Dependence (FTND) total score + FTND classification of dependence
# * Cannabis
#   fjdr100 = age first tried cannabis; 250 = frequency use cannabis; 500 = how many per day; 1050 = ever used cannabis alone
# * Psychotic experiences --> fjpl161-72 
# * Antisocial behavior --> fjaa2000-2207,6700 = knifes, guns and so on 
# * Relationships --> fjpc2100 = how easy to discuss problems in family
# * Life events --> fjle100-75
# ------------------------------------------------------------------------------
# apoa1_tf4 = Apolipoprotein AI; apob_tf4 = Apolipoprotein B mg/dl
# dha_tf4 = Docosahexaenoic acid (DHA) omega-3 fatty acid
# hb_tf4 = haemoblobin

# MISSING 
# - CIS-R symptoms! e.g. weight changes, somatic symptoms, fatigue, concentration, sleep problems, irritability...
# FJMR025a = Fat percentage (%) 
# FJDX141 = Gynoid: fat mass (g) 
# (Baseline and peak vessel diameter, arterial distensibility)

# ------------------------------------------------------------------------------
dep_18.7y <- extr_dep('cct27',  0:12, age = 'cct9991a') # (c in years), cct2715 = total score
# cct4105 exercise frequency 
# (cct5000,01,05,10-15 smoking)
# (cct5020,25,30-39 drinking)
# (cct5050-52, 55,56,71 canabis) and other drugs...
# MISSING eating disorders (cct41..)
# ------------------------------------------------------------------------------
dep_21.9y <- extr_dep('ypa2', 0:12, age = 'ypa9020', plus=0)
# ypa1018 if they had a child
# ------------------------------------------------------------------------------
dep_22.9y <- extr_dep('ypb5',  0:17, remove=c(3,8,12,15,17), age = 'ypb9992', plus=0) # YPB5180 total score 

# ypb1213_imputeno = Ever been diagnosed with hypertension (high blood pressure) - Silent no's included
# ypb1214 = Ever been diagnosed with heart attack/myocardial infarction
# ypb1229 = Ever been diagnosed with type 1 diabetes (juvenile onset diabetes)
# ypb1230 = Ever been diagnosed with type 2 diabetes (adult onset diabetes)
# ypb1231 = schizophrenia / 32 bipolar / 33 depression
# received help for ADHD 52-55
# ypb2050 = Level of physical activity, relative to others of a similar age]
# ypb6210 = Since age 21, whether had an abortion (including partner) and affect this had
# ypb8000-... childhood trauma 
# missing(# smoking and alcohol cannabis)
# ------------------------------------------------------------------------------
dep_23.8y <- extr_dep('ypc16', 50:67, remove=c(3,8,12,15,17), age = 'ypc2650')

# ypc0593 = Respondent is always optimistic about their future 
# ypc1828 = Respondent has ever witnessed a sudden, violent death (eg murder, suicide, or aftermath of an accident)
# ypc1829 = Respondent has ever experienced the sudden, unexpected death of someone close to them
# YPC1830 =  Respondent has ever experienced any other very traumatic or extremely stressful event
# ypc2330 / 80 = miscarriage
# MISSING other depression items  + eating , smoking 
# ------------------------------------------------------------------------------

# Combine left and right mean IMT measurements (also min and ax available)
full$fkcv_imt <- rowMeans(full[,c('fkcv1131','fkcv2131')], na.rm=T)
# 24.5 anthropometry needs to be recoded into cm
full$fkms_height    <- full[,'fkms1000']/10
full$fkms_waistcirc <- full[,'fkms1052']/10
# Calculate waist to hip ratio
full$fkms_whr <- full[,'fkms_waistcirc'] / (full[,'fkms1062']/10) # Hip circumference (trans in cm)

# Calculate cardiac measurements (see formulas above)
full$fkec_LVM = 0.8 * ( 1.04 * ( rowSums(full[,c('fkec5080',  # Average LV internal End Diastolic Dimension (cm)
                                                'fkec5050',  # Average thinkness LV interventricular septum in diastole (cm)
                                                'fkec5180')] # Average LV posterior wall Diastolic thickness (cm)
                                         )^3 - full[,'fjgr043']^3 ) ) + 0.6
full$fkec_RWT <- 2 * full[,'fkec5180'] / full[,'fkec5080']
full$fkec_FS <- abs( (full[,'fkec5080'] - full[,'fkec5090']) / full[,'fkec5080'] * 100 ) 

cmr_24.5y <- extr_cmr(data.frame('fkms_height'='height', # cm
                                 'fkms1030'='weight',    # kg
                           'fkms_waistcirc'='waist_circ',# cm
                              'fkms_whr'='waist_hip_ratio',  # cm
                              'fkdx1001'='total_fatmass',  # grams (--> transformed later)
                              'fkdx1002'='total_leanmass', # grams (--> transformed later)
                              'fkdx1031'='trunk_fatmass',  # grams (--> transformed later)
                              'fkdx1041'='android_fatmass',# grams (--> transformed later)
                              'fkli1010'='liver_fat',
                                        # arteries
                              'fkbp1030'='SBP', # Average, seated (mmHg)
                              'fkbp1031'='DBP', # Average, seated (mmHg)
                              'fkcv_imt'='IMT', # Baseline (end diastole) intima-media thickness (mm)
                              'fkcv4200'='PWV', # carotid to femoral (m/s) average
                            # 'fksp1804'='eject_dur', # Ejection Duration, ms
                                        # heart 
                              'fksp1837'='heart_rate',
                              'fkec_LVM'='LVM', # Left ventricular mass (LVM) (g)
                              'fkec_RWT'='RWT', # Relative wall thickness (RWT) (cm)
                               'fkec_FS'='FS',  # Fractional shortening (%)
                                   # no Average left atrial size, m-mode (cm)
                                   # no Average left ventricular internal diameter during diastole (cm), dimention ok?
                              # assuming the age is the same here 
                              'chol_f24'='tot_chol', # mmol/l
                               'hdl_f24'='HDL_chol', # mmol/l
                               'ldl_f24'='LDL_chol', # mmol/l
                           'insulin_f24'='insulin',  # mU/L
                              'trig_f24'='triglyc',  # mmol/l,
                           'glucose_f24'='glucose',  # mmol/l
                             # 'glc_f24'='glucose1', # for consistency
                               'crp_f24'='CRP', 
                        # other
                        'fjal4000'='alcohol', # Alcohol use disorder identification test (AUDIT) score
                        'fjsm_tot'='smoking', # Average number of cigarettes per week
                        'fjdr4500'='canabis'),# Cannabis Abuse Screening Test (CAST) score
                   
                   age = 'fkar0010') # (11 in years)

# Also available, bu not included: ---------------------------------------------

# * Echocardiography:
#   fkec5050 # Average LV interventricular septum diastolic thickness (cm); fkec5060 # Average LV interventricular septum systolic thickness (cm)
#   fkec5100 = Average LV outflow tract (LVOT) diameter (cm)
#   fkec5180 # Average LV posterior wall diastolic thickness (cm); fkec5190 # Average LV posterior wall systolic thickness (cm)
# * Tissue Doppler imaging (TDI):
#   fkec5240: Average LV lateral wall a' velocity (cm/s); fkec5250: Average LV lateral wall e' velocity (cm/s); fkec5260: Average LV lateral wall s' velocity (cm/s)
#   fkec5270: Average RV a' (cm/s); fkec5280: Average RV e' (cm/s); fkec5290: Average RV s' (cm/s) # TDI
#   fkec5300: Average LV a' septum (cm/s); fkec5310: Average LV e' septum (cm/s); fkec5320: Average LV s' septum (cm/s) # TDI
# * fkco3002 = uses insuline treatment

# (fkac1010-1110 accelerometer data)
# (fkca1025 = canabis use score)
# (fkde1040,045,230,260 sexual abuse)

# (fkpl[...] psychotic experiences)

# MISSING 
# fkal[...] = alcohol use (AUDIT score)
# fksm[...] = smoking frequency
# FKDX1051: Gynoid fat mass (g): F@24

# CIS-R DEPRESSION SYMPTOMS! (fkdq...) ------------
# fkdq1000 = Mild, 1010 = moderate, 1020 severe depressive episode
# fkdq1030 = GAD, 50 = social phobia, 60 = phobia, 70 = panic disorder, 80 = panick attack symptoms tot

# FKDQ5000: In past month, YP had a spell of feeling sad/miserable/depressed: F@24
# FKDQ5010: In past 7 days, YP felt sad/miserable/depressed: F@24
# FKDQ5020: In past month, YP has been able to enjoy things as much as usual: F@24
# FKDQ5030: In past 7 days, YP has been able to enjoy things as much as usual: F@24
# FKDQ5040: In past 7 days, no. of days YP felt sad/depressed/uninterested in things: F@24
# FKDQ5050: In past 7 days, YP felt sad/depressed/uninterested for >3 hours in a day: F@24
# FKDQ5060: In past week, main reason for feeling sad/depressed/uninterested: F@24
# FKDQ5070: In past 7 days, YP felt happier when something nice happened/in company: F@24
# FKDQ5080: Length of time YP has felt sad/miserable/depressed/uninterested in things: F@24
# FKDQ5100: In past 7 days, time of day felt more sad/depressed/uninterested: F@24
# FKDQ5110: In past month, YP's interest in sex increased, decreased or stayed same: F@24
# FKDQ5120: In past 7 days, when sad/uninterested YP has been restless/not sit still: F@24
# FKDQ5130: In past 7 days, when sad/uninterested YP has done things slower than usual: F@24
# FKDQ5140: In past 7 days, YP felt guilty/blamed themselves, even if not their fault: F@24
# FKDQ5150: In past 7 days, YP felt they are not as good as other people: F@24
# FKDQ5160: In past 7 days, YP has felt helpless (e.g., about future): F@24
# FKDQ5170: In past 7 days, YP has felt that life is not worth living: F@24
# FKDQ5180: In past 7 days, YP has had thoughts of harming themselves: F@24
# FKDQ5190: In past 7 days, YP has thought about a way to kill themselves: F@24
# FKDQ6010: YP has had any worries In past month: F@24
# FKDQ6200: In past month, YP felt anxious about situations where no real danger: F@24
# FKDQ6410: Length of time YP had felt anxious/nervous/tense: F@24

# BUT MISSING
# other symptoms 

# ==============================================================================
#cat(ls(), sep=", ")

data <- cbind(data, 
     mdep_09.6y, cmr_09.6y, cmr_09.8y, 
      dep_10.6y, cmr_10.6y, cmr_10.7y, 
     mdep_11.7y, cmr_11.7y, #_1, cmr_11.7y_2, 
      dep_12.8y, cmr_12.8y, 
     mdep_13.1y, cmr_13.1y, 
      dep_13.8y, cmr_13.8y, 
      cmr_15.4y, 
      dep_16.6y, mdep_16.7y, cmr_16.0y, cmr_17.0y,
      dep_17.8y, cmr_17.8y, 
      dep_18.7y, 
      dep_21.9y, 
      dep_22.9y, 
      dep_23.8y, 
      cmr_24.5y)


# Select groups of variables  
sel <- function(var, times=NULL){
  subs = names(data)[grep(paste(var, collapse='|'), names(data))]
  if (!is.null(times)) { subs = subs[grep(paste(paste0('_',times), collapse='|'), subs)] }
  # print(summary(data[,subs]))
  return(subs)
}

# Convert fat/lean mass variables into kg (variances are too large otherwise)
data[, sel('mass')] <- data[, sel('mass')] / 1000

# ------------------------------------------------------------------------------
# remove all unnecessary attributes that may corrupt the file
for (var in colnames(data)) {
  attr(data[,deparse(as.name(var))], "names") <- NULL
  attr(data[,deparse(as.name(var))], "value.labels") <- NULL
}

# save and run -----------------------------------------------------------------

crm <- cor(data[,-c(1,2)], method='spearman', use='pairwise.complete.obs')
# vcm <- cov(data[,-c(1,2)], method='spearman', use='pairwise.complete.obs')

saveRDS(data, 'raw_data.rds'); write.csv(data,'raw_data.csv')

# write.csv(vcm, 'varcov_matrix.csv')
write.csv(crm, 'corr_matrix.csv')

dev.off()
