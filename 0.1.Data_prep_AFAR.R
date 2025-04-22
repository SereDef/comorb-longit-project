# devtools::install_github("SereDef/ALSPAC.helpR")
library(ALSPAC.helpR)
library(dplyr)

# File name and location 
filename <- 'Cecil_B4195_07Dec23_2.sav'

# Read in the file
data <- load_alspac( # open in interactive session
                    lower.case = TRUE,
                    keep.value.labels = TRUE,
                    load.metadata = FALSE)

# Initiate dataframe and set up unique child ID
dset <- data.frame('IDC' = make_idc(mom.id='cidb4195'),
                   'sex' = data$kz021, # Participant sex (Male, Female)
              'ethnicity'= data$c804)  # Child ethnic background (White, Non-white)

# ==============================================================================
# -------------------------------- Functions -----------------------------------
to_numeric <- function(vars) {
  data[,vars] <<- sapply(data[,vars], function(x) { 
    if (!is.numeric(x)) as.numeric(levels(x))[x] else x 
    })
}

winsorize <- function(x, xmin=min(x, na.rm=TRUE), xmax=max(x, na.rm=TRUE)) {
  message('Original min and max values: ', min(x, na.rm=TRUE), ' and ', max(x, na.rm=TRUE))
  
  xmin <- if (is.character(xmin)) quantile(x, as.numeric(xmin), na.rm=TRUE) else xmin
  xmax <- if (is.character(xmax)) quantile(x, as.numeric(xmax), na.rm=TRUE) else xmax
  
  cat('Winsorising', sum(x < xmin, na.rm=TRUE), 'values below', xmin, 'and',
      sum(x > xmax, na.rm=TRUE), 'values above', xmax,'\n\n')
  
  x[x < xmin] <- xmin
  x[x > xmax] <- xmax
  
  hist(x)
  print(summary(x))
}

# Extract and clean depression items
extr_dep <- function(pref, nrs, remove=c(), age, plus=NULL, reporter='s'){
  
  # Remove positively worded items (not in SMFQ)
  if (length(remove) > 0)  nrs <- nrs[-remove]
  
  # Add leading 0s to the number sequence & if necessary, add plus at the end
  n_adj <- paste0(sprintf('%02d', nrs), plus)
  
  dep <- data %>%
    # Identify variables 
    select(paste0(pref, n_adj)) %>%
    # any other variables left ? 
    # print(setdiff(names(data)[grep(pref, names(data))], names(dep)))
    
    mutate(across(everything(), function(x) { 
      # Drop empty levels
      droplevels(x)
      # Transform response to 0=not true; 1=sometimes true; 2=true
      x <- case_when((grepl('Not', x, ignore.case=TRUE)) ~ 0, # e.g., 'Not true', 'Not at all' ...
                     (grepl('Sometimes', x, ignore.case=TRUE)) ~ 1, # e.g. 'Sometimes', 'Sometimes true' ...
                     (grepl('True', x, ignore.case=TRUE)) ~ 2, # e.g. 'True', '3. True' ...
                     TRUE ~ NA)
    }))
  
  # Rename depression items 
  names(dep) <- paste0(reporter,'DEP', sprintf('%02d', 1:13))
  
  # Check recoding ( most answers == 0 --> average should be lower than 1 )
  if (any(colMeans(dep, na.rm=TRUE) > 1)) { 
    cat('ATTENTION: check the direction of answer scoring!')
    print(summary(dep)) 
  }
  # Check the number of items is correct
  if (ncol(dep) != 13) cat('ATTENTION: There should be 13 items in this score, but ',
                           ncol(dep), ' where included!')
  
  # Compute total score 
  score_var <- paste0(reporter,'DEP_score')
  dep[,score_var] <- ifelse(rowSums(is.na(dep))!=ncol(dep), rowSums(dep, na.rm=TRUE), NA)
  
  # Add age variable (in months) and display its descriptives
  age_var = paste0(reporter,'DEP_age')
  dep[,age_var] <- if (!is.numeric(data[,age])) as.numeric(levels(data[,age]))[data[,age]] / 12 else data[,age] / 12
  
  # Print summary
  print(summary(dep[,c(score_var, age_var)]))
  
  # Histogram 
  hist(dep[,score_var], 
       main=paste(score_var, round(median(dep[,age_var], na.rm=TRUE),1)), 
       breaks=max(dep[,score_var], na.rm=TRUE), 
       xlim=c(0, 27))
  
  # Add age mark to variable names
  names(dep) <- paste0(names(dep),'_',
                       round(median(dep[,age_var], na.rm=TRUE),1),'y')
  
  return(dep)
}


# Compute BMI given height in cm and weight in kg 
bmi_calc <- function(h, w, title='BMI', age=''){
  
  # Clean up a little
  #h <- winsorize(h, xmin='0.001', xmax='0.999') 
  #w <- winsorize(w, xmin='0.001', xmax='0.999') 
  
  h = h/100 # transform to meters
  bmi <- w / h^2
  # print(summary(bmi))
  hist(bmi, main=paste(title, age), breaks=max(bmi, na.rm=T))
  return(bmi)
}

# Extract and clean cardio-metabolic items
extr_cmr <- function(df, age, outlier_cutoff=5){
  # identify variables 
  cmr <- as.data.frame(data[, names(df)])
  
  cmr <- data %>%
    # Identify variables 
    select(names(df)) %>%
    
    mutate(across(everything(), function(x) { 
      if (!is.numeric(x)) {
        # Drop empty levels
        droplevels(x)
        # Transform to numeric
        x <- suppressWarnings(as.numeric(levels(x))[x])
      }
      return(x)
    }))
  
  # rename them
  names(cmr) <- df
  
  # Remove outliers (set to NA)
  # except <- which(names(cmr) %in% c('insulin','glucose','triglyc','CRP','IL_6','canabis'))
  # if (identical(except, integer(0))) { 
  #   cmr <- as.data.frame(sapply(cmr, rm_outliers, cutoff = outlier_cutoff))
  # } else { 
  #   cmr[,-except] <- as.data.frame(sapply(cmr[,-except], rm_outliers, cutoff = outlier_cutoff)) }
  
  # Add age variable (in months) and display its descriptives
  cmr[,'CMR_age'] <- if (!is.numeric(data[,age])) as.numeric(levels(data[,age]))[data[,age]] / 12 else data[,age] / 12
  # print(summary(cmr$CMR_age))
  
  # Transform mass from grams to kg 
  cmr[, grep('mass', names(cmr))] <- cmr[, grep('mass', names(cmr))] / 1000
  
  # Compute BMI if all info is available
  if ('height' %in% names(cmr)) {
    for (to_index in c('weight','total_fatmass','total_leanmass')) {
      indexed <- if (to_index=='weight') 'BMI' else if (to_index=='total_fatmass') 'FMI' else 'LMI'
     
      if (to_index %in% names(cmr)) {
        message('Calculating ', indexed, '...')
        
        cmr[,indexed] <- bmi_calc(cmr[,'height'], cmr[,to_index], 
                                  title=indexed, 
                                  age=round(median(cmr[,'CMR_age'], na.rm=TRUE),1))
      }
    }
  }
  

  # Add age mark to variable 
  names(cmr) <- paste0(names(cmr),'_',round(median(cmr$CMR_age,na.rm=T),1),'y')
  
  print(summary(cmr))
  return(cmr)
}

# Extract and combine puberty variables
extr_pub <- function(age, scale=12){
  
  pubid <- substr(age,1,4)
  
  pubnames <- names(data)[(grepl(paste0('^',pubid), names(data), ignore.case=TRUE))]
  
  renames <- c(pub_age = age, phys_activity = paste0(pubid,'09'))
  
  pubset <- data %>%
    # Identify variables + add sex for puberty stage
    select(any_of(c(pubnames, 'kz021'))) %>% 
    
    # Recode all columns to numeric 
    mutate_at(vars(matches('09$|30$|35$|50$|55$')), ~ as.numeric(.)) %>%
    mutate_at(vars(age), ~ as.numeric(as.character(.))/scale) %>% 
    
    # Compute puberty stage for male and female participants and age
    mutate(puberty_stage = na_if(case_when(
      kz021 == 'Female' ~ rowSums(select(., matches('30$|35$')), na.rm = TRUE), # Breasts + pubic hair
      kz021 == 'Male' ~ rowSums(select(., matches('50$|55$')), na.rm = TRUE), # Testes,scrotum penis + pubic hair
      TRUE ~ NA), 0)) %>% 
    
    rename(all_of(renames)) %>%
    select(pub_age, phys_activity, puberty_stage)
  
  # Add age mark to variable
  names(pubset) <- paste0(names(pubset),'_',
                          round(median(pubset$pub_age, na.rm=TRUE), 1),'y')

  print(summary(pubset))
  return(pubset)
}

extr_sleep <- function(bedvar, wakevar, age, scale=12) {
  
  dt <- data.frame(matrix(nrow=nrow(data), ncol=0))
  
  for (t in c(bedvar, wakevar)) {
    # Transform to numeric
    multiv = strsplit(t, '/')[[1]]
    if (length(multiv)>1) {
      t.h <- as.numeric(as.character(data[,multiv[1]]))
      t.m <- as.numeric(as.character(data[,multiv[2]]))
    } else {
      t.h <- as.numeric(as.character(data[,paste0(t,'a')]))
      t.m <- as.numeric(as.character(data[,paste0(t,'b')]))
    }
    
    # Reformat hours for bed and wake times
    if (t == bedvar) {
      t.h <- ifelse(is.na(t.h), NA, 
                    ifelse(t.h >= 24, t.h-24, ifelse(t.h==12, 0, # midnight
                                                     # 4 assumed 4pm, 3 assumed 3am
                                                     ifelse(t.h >= 4 & t.h < 12, t.h + 12, t.h))))
      
      day <- ifelse(t.h < 10, '2020-01-02', '2020-01-01') # random day
      name <- 'bed_time'
    } else {
      t.h <- ifelse(is.na(t.h), NA, # 2 assumed 2pm, 3 assumed 3am
                    ifelse(t.h < 3, t.h + 12, sprintf('%02d', t.h)))
      day = '2020-01-02' # random day + 1 day
      name <- 'wake_time'
    }
    # Reformat minutes
    t.m <- ifelse(is.na(t.m) | (t.m > 59), NA, sprintf('%02d', t.m))
    
    time <- ifelse(is.na(t.h), NA, 
                   ifelse(is.na(t.m), paste(paste(day,t.h), '00', sep = ':'), 
                          paste(paste(day,t.h), t.m, sep = ':')))
    
    dt[, name] <- time
  }
  
  dt[,'sleep_dur'] <- -as.numeric(difftime(dt[,'bed_time'], dt[,'wake_time'], units='hours'))
  
  # Strip random date from bed and wake times: remove everything before the first space
  dt[,'bed_time'] <- sub(".*\\s", "", dt[,'bed_time'])
  dt[,'wake_time'] <- sub(".*\\s", "", dt[,'wake_time'])
  
  # Add age mark to variable
  names(dt) <- paste0(names(dt),'_',
                      round(median(as.numeric(as.character(data[,age]))/scale, na.rm=TRUE), 1),'y')
  
  print(summary(dt))
  
  return(dt)
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

# ----- Removed from the scoring -----
# 2 Had fun; 8 Felt happy; 11 Enjoyed doing lots of things; 17 Had a good time; 
# 18 Laughed a lot; 19 Looked forward to the day ahead; 20 Felt really positive 
# about the future; 21 Felt valued; 22 Felt unhappy

# ------------------------------------------------------------------------------

# pdf('histogr_data_AFAR.pdf') # NOTE: saved in working directory :) 

# f7ms010 Height (cm): F7  7.42 years [6.83 - 9.42] f7003c
# f7ms018 Waist circumference (cm): F7  
# f7ms020 Hip circumference (cm): F7  
# f7ms026 Weight (kg): F7 

sink('./Data_prep_AFAR_report.txt')

cat('\n================================ 8 YEARS =====================================\n')

pub_08.2y <- extr_pub(age='pub195')

# ALSO : kt6.. diet (FFQ)
# f8lf020 Child height (cm): lung function: F8  8.58 years [7.42 - 10.58] f8003c/f8006c
# f8lf021 Child weight (kg): lung function: F8 

cat('\n================================ 9 YEARS =====================================\n')

mdep_09.6y <- extr_dep('ku6', 60:72, reporter='m', age = 'ku991a') 

drugs_09.6y <- data %>%
            # Child has tried cigarettes
  transmute(smoking_9.6y = case_when(ku431 =='Do not think so' ~ 0, 
                                     ku431 =='Possibly' ~ 0.25, 
                                     ku431 =='Probably' ~ 0.75,
                                     ku431 =='Yes, parent knows' ~ 1, TRUE ~ NA), 
            # Child has tried alcohol
            alcohol_9.6y = case_when(ku430 =='Do not think so' ~ 0, 
                                     ku430 =='Possibly' ~ 0.25, 
                                     ku430 =='Probably' ~ 0.75,
                                     ku430 =='Yes, parent knows' ~ 1, TRUE ~ NA))

lapply(drugs_09.6y, table, useNA='ifany')

sleep_09.6y <- extr_sleep(bedvar='ku341', wakevar='ku340', age = 'ku991a') # school days

mdep_09.6y <- cbind(mdep_09.6y, drugs_09.6y, sleep_09.6y) 

# ALSO : SDQ - symptoms + scores (ku6-7...)
#        Sleep - bed and wake time on weekends (ku342-3) + sleep quality questions
#        Child has tried drugs (ku432) [ONLY 5 CASES]
#        Age child tried alcohol / cigarettes / drugs (ku435-7) [years]
#        Discipline (ku2-3...)

pub_09.6y <- extr_pub(age = 'pub295')

# ------------------------------------------------------------------------------
cmr_09.8y <- extr_cmr(data.frame('pub203'='height',  # cm; NOTE: age = 9.6 years (pub295)
                                'f9dx010'='weight', # kg (DXA)
                                'f9dx135'='total_fatmass', # grams (--> trasformed later)
                                 'crp_f9'='CRP'), # mg/L
                      age = 'f9003c')

# ALSO : Parathyroid hormone (pth_f9); Vitamin D (vitd*_f9)

# MISSING : Height, weight, waist and hip circumference (f9ms...)
#           DXA lean mass (f9dx136)
#           BP and pulse (f9sa...)
#           Blood markers: chol, hdl, ldl, insulin, trig; il6 (*_f9)


# cmr_09.8y[,'BMIw_9.8y'] <- winsorize(cmr_09.8y[,'BMI_9.8y'], xmax='0.995') #, xmin='0.001', xmax='0.999')
# cmr_09.8y[,'FMIw_9.8y'] <- winsorize(cmr_09.8y[,'FMI_9.8y'], xmax='0.995') #, xmin='0.001', xmax='0.999')

cat('\n================================ 10 YEARS ====================================\n')

dep_10.6y <- extr_dep('fddp1', 10:25, remove=c(2,8,11), age = 'fd003c') # Depression score = fddp130 (prorated?)

# ------------------------------------------------------------------------------
data <- data %>%
         # Child smoked cigarettes
  mutate(fd_smoke = case_when(fdaa482 =='No' ~ 0, fdaa482 =='Yes' ~ 1, TRUE ~ NA), 
         # Child drunk alcohol
         fd_alcol = case_when(fdaa492 =='No' ~ 0, fdaa492 =='Yes' ~ 1, TRUE ~ NA))

cmr_10.6y <- extr_cmr(data.frame('fdms010'='height', # Height (cm)
                                 'fdms026'='weight', # Weight (kg)
                                 'fdms018'='waist_circ', # Waist circumference (cm)
                                 'fdar117'='SBP',  # Systolic blood pressure (mmHg)
                                 'fdar118'='DBP',  # Diastolic blood pressure (mmHg)
                                 'fdar114'='PWV',  # Pulse Wave Velocity (m/s) # radial - carotid!!
                                'fd_smoke'='smoking',
                                'fd_alcol'='alcohol'), 
                      age = 'fd003c')

pub_10.7y <- extr_pub(age = 'pub397a')

# ALSO : DAWBA - symptoms + diagnoses (kv...) [10.7 years] 
#        Freq child drunk alcohol / No. Times per week (fdaa493-4)
#        Freq child smoked cigarettes / No. cigarettes per week / parents permission (fdaa483-5)
#        Child smoked cannabis / Freq / No. times per week (fdaa522-4) [10.6 years] # ONLY 1 CASE
#        Friends drunk alcohol / freq (fdaa490-1); Alcohol intake g (fddd340) / g/day (fddd275: diet diary)

cat('\n================================ 11 YEARS ====================================\n')

mdep_11.7y <- extr_dep('kw60', 0:12, reporter='m', age = 'kw9991a') # tot scores kw6100a / kw6100b  + SDQ scores kw6600 - kw6605

mdep_11.7y <- cbind(mdep_11.7y, 
                    extr_sleep(bedvar='kw4061', wakevar='kw4060', age='kw9991a')) # school days

# ALSO : SDQ - symptoms + scores (kw65-66...)

# ------------------------------------------------------------------------------

to_numeric(c('fems026','fedx016','fems018','fems020'))

data <- data %>%
         # Combine two weight measures
  mutate(fe_weight = rowMeans(.[, c('fems026','fedx016')], na.rm=TRUE), 
         # Calculate waist to hip ratio
         fe_WHR = fems018 / fems020, # Waist / hip circumference (cm)
         # Child drunk alcohol
         fe_alcol = case_when(febp191 =='No' ~ 0, febp191 =='Yes' ~ 1, TRUE ~ NA))

cmr_11.7y <- extr_cmr(data.frame('fems010'='height', # Height (cm)
                               'fe_weight'='weight', # Weight (kg) also fedx016 DXA weight (Kg)
                                 'fems018'='waist_circ', # Waist circumference (cm)
                                  'fe_WHR'='waist_hip_ratio', # cm 
                                 'fedx135'='total_fatmass',  # grams (--> trasformed later)
                                 'fedx136'='total_leanmass', # grams (--> trasformed later)
                                'fems028a'='fat_percent', # Fat percentage, impedance
                                 'fesa021'='SBP',  # Mean BP systolic
                                 'fesa022'='DBP',  # Mean BP diastolic
                                 'fesa023'='heart_rate', #  Mean Pulse == heart rate? beats per minute
                                'fe_alcol'='alcohol'),
                       age = 'fe003c')

# ALSO : DXA trunk fat mass (fedx126)
#        Got really drunk on alcohol - How many times: F11 (febp192)

pub_11.7y <- extr_pub(age = 'pub497a')

cat('\n================================ 12 YEARS ====================================\n')

dep_12.8y <- extr_dep('ff65', 0:15, remove=c(2,8,11), age = 'ff0011a')

# ALSO : Depression Computer Task - symptoms (ff42...)

# ------------------------------------------------------------------------------

to_numeric(paste0('ff26', c(20:22,25:27)))

data <- data %>%
         # Combine two BP measures
  mutate(ff_SBP = rowMeans(.[, c('ff2620','ff2625')], na.rm=TRUE), 
         ff_DBP = rowMeans(.[, c('ff2621','ff2626')], na.rm=TRUE), 
         # Combine two pulse measures
         ff_pulse = rowMeans(.[, c('ff2622','ff2627')], na.rm=TRUE), 
         # Teenager has smoked cigarettes
         ff_smoke = case_when(ff6610 =='No' ~ 0, ff6610 =='Yes' ~ 1, TRUE ~ NA))

cmr_12.8y <- extr_cmr(data.frame('ff2000'='height', # Height (cms) TF1
                                 'ff2030'='weight', # Weight (Kgs) TF1
                                 'ff2020'='waist_circ', # Waist circumference (cms) TF1 
                                 'ff2036'='fat_percent', # Fat percentage (%)  TF1
                                 'ff_SBP'='SBP',  # Mean BP systolic
                                 'ff_DBP'='DBP',  # Mean BP diastolic
                               'ff_pulse'='heart_rate', #  Mean pulse (bpm)
                               'ff_smoke'='smoking'),
                      age = 'ff0011a')

# ALSO : Teenager has smoked cigarettes, past 6 months / Freq / N. times / at least once a week /
#        Number cigarettes teenager smokes per week / parents permission / age first and last smoked (ff6611-18)

cat('\n================================ 13 YEARS ====================================\n')

mdep_13.1y <- extr_dep('ta50', 20:32, reporter='m', age = 'ta9991a') 

# ------------------------------------------------------------------------------

to_numeric(c('ta6000','pub503','ta6001','pub504','ta9991a','pub597a'))

data <- data %>%
  # Combine (self-reported) height and weight measures
  mutate(ta_height = rowMeans(.[, c('ta6000','pub503')], na.rm=TRUE), 
         ta_weight = rowMeans(.[, c('ta6001','pub504')], na.rm=TRUE), 
         # Combine two measurement ages 
         ta_mean_age = rowMeans(.[, c('ta9991a','pub597a')], na.rm=TRUE), 
         # Frequency teenager drinks wine
         ta_alcol1 = case_when(ta8280 =='Not at all' ~ 0, 
                               ta8280 =='< Once a week' ~ 0.5, 
                               ta8280 =='Once a week' ~ 1, 
                               ta8280 =='> Once a week' ~ 2, TRUE ~ NA), 
         # Frequency teenager drinks beer/lager
         ta_alcol2 = case_when(ta8281 =='Not at all' ~ 0, 
                               ta8281 =='< Once a week' ~ 0.5, 
                               ta8281 =='Once a week' ~ 1, 
                               ta8281 =='> Once a week' ~ 2, TRUE ~ NA),
         # Frequency teenager drinks spirits
         ta_alcol3 = case_when(ta8282 =='Not at all' ~ 0, 
                               ta8282 =='< Once a week' ~ 0.5, 
                               ta8282 =='Once a week' ~ 1, 
                               ta8282 =='> Once a week' ~ 2, TRUE ~ NA), 
         # Frequency teenager drinks other alcohol
         ta_alcol4 = case_when(ta8283 =='Not at all' ~ 0, 
                               ta8283 =='< Once a week' ~ 0.5, 
                               ta8283 =='Once a week' ~ 1, 
                               ta8283 =='> Once a week' ~ 2, TRUE ~ NA)) %>%
  # Combine alcohol measurements 
  mutate(ta_alcol = rowSums(.[, c('ta_alcol1','ta_alcol2','ta_alcol3','ta_alcol4')], na.rm=TRUE) * 
           ifelse(rowSums(is.na(.[, c('ta_alcol1','ta_alcol2','ta_alcol3','ta_alcol4')])) == 4, NA, 1))

cmr_13.1y <- extr_cmr(data.frame('ta_height'='height', # cm (parent-reported)
                                 'ta_weight'='weight', # kg (parent-reported)
                                  'ta_alcol'='alcohol'),
                      age = 'ta_mean_age')

pub_13.1y <- extr_pub(age = 'pub597a')

# ALSO : SDQ - symptoms + scores (ta7000-7025)
#        Health - e.g. sleep quality and tired/lacking energy (ta11-12..)
#        Best description of teenager's alcohol drinking (ta8285)

# ------------------------------------------------------------------------------

dep_13.8y <- extr_dep('fg72', 10:25, remove=c(2,8,11), age = 'fg0011a')

# ------------------------------------------------------------------------------

to_numeric(c('fg3130','fg3207', # weight
             paste0('fg10',20:22), paste0('fg61', c(20:22,25:27)), # BP and pulse
             'fg1550')) # alcohol

data <- data %>%
         # Combine multiple measurements
  mutate(fg_weight = rowMeans(.[, c('fg3130','fg3207')], na.rm=TRUE), # visit and DXA
         fg_SBP = rowMeans(.[, c('fg6120','fg6125','fg1021')], na.rm=TRUE), 
         fg_DBP = rowMeans(.[, c('fg6121','fg6126','fg1022')], na.rm=TRUE), 
         fg_pulse = rowMeans(.[, c('fg6122','fg6127','fg1020')], na.rm=TRUE), # Pulse and heart rate at rest
         # Alcohol intake (g) diet diary mean
         fg_alcol = case_when(fg1550 == 0 ~ 0, fg1550 > 0 ~ 1, TRUE ~ NA), 
         # Teenager has smoked cigarettes
         fg_smoke = case_when(fg4822 =='No' ~ 0, fg4822 =='Yes' ~ 1, TRUE ~ NA))
        
cmr_13.8y <- extr_cmr(data.frame('fg3100'='height', # Height (cms): TF2
                              'fg_weight'='weight', # Weight (Kgs): TF2 and DXA
                                 'fg3120'='waist_circ', # Waist circumference (cms): TF2
                                 'fg3254'='total_fatmass',  # grams (--> trasformed later)
                                 'fg3255'='total_leanmass', # grams (--> trasformed later)
                                 'fg3257'='android_fatmass',# grams (--> trasformed later)
                                 'fg3260'='gynoid_fatmass', # grams (--> trasformed later)
                                 'fg3136'='fat_percent', # Fat percentage (%): TF2
                                 'fg_SBP'='SBP', # TF2
                                 'fg_DBP'='DBP', # TF2
                               'fg_pulse'='heart_rate', # at rest (BPM)
                               'fg_alcol'='alcohol',
                               'fg_smoke'='smoking'),
                      age = 'fg0011a')

# ALSO : DAWBA - symptoms + diagnoses (tb4-8...) [13.8 years]
#        Personality (tb30..)
# ALSO : DXA trunk fat mass (fg3245), total bone mass (fg3253)
#        Average HR during exercise (fg1026); Post exercise BP and HR (fg1024,25,27)
#        Smoked in past 6 months / Freq / N. cigarettes per week / parents permission /
#        / still smokes / ... (fg4823-30)

cat('\n================================ 14 YEARS ====================================\n')

cmr_14.2y <- data %>%
  # Combine smoking measurements
  transmute(smoking_14.2y = case_when(ccr705 == 'Usually 1+ every day' ~ 6, 
                                      ccr705 == 'Usually >6 a week, but not every day' ~ 5, 
                                      ccr705 == 'Usually 1-6 a week' ~ 4, 
                                      ccr705 == 'Sometimes, but < 1 a week' ~ 3, 
                                      ccr705 == 'Used to sometimes, but never now' ~ 2, 
                                      ccr705 == 'Only ever tried once or twice' ~ 1, 
                                      ccr700 == 'No' ~ 0,
                                      TRUE ~ NA), 
            canabis_14.2y = case_when(ccr765 == 'Usually every day' ~ 6, 
                                      ccr765 == 'Usually >6 times a week, but not every day' ~ 5, 
                                      ccr765 == 'Usually 1-6 times a week' ~ 4, 
                                      ccr765 == 'Sometimes, but < once a week' ~ 3, 
                                      ccr765 == 'Used to sometimes but never now' ~ 2, 
                                      ccr765 == 'Once ever tried once or twice' ~ 1, 
                                      ccr760 == 'No' ~ 0,
                                      TRUE ~ NA),
            CMR_age_14.2y = as.numeric(as.character(ccr991a)) /12 )
summary(cmr_14.2y)

pub_14.7y <- extr_pub(age = 'pub697a')

# ALSO : Psychosis (ccr3-4...)
#        Age when first smoked a cigarette (ccr710) / Total number of cigarettes (ccr730) /
#        / smoked any cigarettes since 14th birthday (ccr735) / N cigarettes per day (ccr740)
#        / nicotine patches or gums...(ccr751-3)
#        Age when first tried cannabis (ccr770) / Total number of times used cannabis (ccr775) /
#        / type of cannabis and method (ccr780-795) / other cannabis questions... (ccr8...)
#        Other drugs: amphetamines, ecstasy, LSD ...(ccr85...)

cat('\n================================ 15 YEARS ====================================\n')

pub_15.3y <- extr_pub(age = 'pub797a') 

to_numeric(c('fh3000','fh2207', 'fh3010','fh2209', # height and weight
             paste0('fh20', c(30:32,35:37)), # BP and pulse
             'fh5440','fh5441')) # sleep

data <- data %>%
         # Combine height and weight measures 
  mutate(fh_height = rowMeans(.[, c('fh3000','fh2207')], na.rm=TRUE), 
         fh_weight = rowMeans(.[, c('fh3010','fh2209')], na.rm=TRUE), 
         # Combine two BP measures
         fh_SBP = rowMeans(.[, c('fh2030','fh2035')], na.rm=TRUE), 
         fh_DBP = rowMeans(.[, c('fh2031','fh2036')], na.rm=TRUE), 
         # Combine two pulse measures
         fh_pulse = rowMeans(.[, c('fh2032','fh2037')], na.rm=TRUE), 
         # Length of time YP sleeps on normal school night: hour & minutes: TF3
         fh_sleep = ifelse(is.na(fh5440), NA, 
                           rowSums(cbind(fh5440, fh5441/60), na.rm=TRUE))
         )

cmr_15.4y <- extr_cmr(data.frame('fh_height'='height', # cm
                                 'fh_weight'='weight', # kg
                                 'fh4020'='waist_circ', # Waist circumference (cms): TF3
                                 'fh2254'='total_fatmass',  # grams (--> trasformed later)
                                 'fh2255'='total_leanmass', # grams (--> trasformed later)
                                 'fh2257'='android_fatmass',# grams (--> trasformed later)
                                 'fh2260'='gynoid_fatmass', # grams (--> trasformed later)
                                 'fh3016'='fat_percent', # Fat percentage (%): TF3
                                 'fh_SBP'='SBP', # BP - systolic TF3
                                 'fh_DBP'='DBP', # BP - diastolic TF3
                               'fh_pulse'='heart_rate', # pulse TF3
                                'crp_tf3'='CRP', # C-reactive Protein mg/l, TF3
                           'cortisol_tf3'='cortisol', # Plasma cortisol (nmol/l) - From morning samples, TF3
                            'serumtg_tf3'='triglycerides', # Serum total triglycerides (mmol/l): TF3
                             'serumc_tf3'='cholesterol', # Serum total cholesterol (mmol/l): TF3 # Free cholesterol (mmol/l): TF3
                              # 'hdlc_tf3'='HDL_chol', # Total cholesterol in HDL (mmol/l): TF3
                              # 'ldlc_tf3'='LDL_chol', # Total cholesterol in LDL (mmol/l): TF3
                              #  'glc_tf3'='glucose', # Glucose (mmol/l): TF3
                               'fh_sleep'='sleep_dur'),
                      age = 'fh0011a') # fh5307 = Age in months of YP at clinic visit: TF3, but similar to fh0011a?

# ALSO : DAWBA - symptoms + diagnoses (fh63-8...); Adult DAWBA (fh95...) [15.4 years]
#        Sleep questionnaire: sleep duration in the weekend (fh5480-1) /
#        / Time go to sleep / Time takes to fall asleep } week and weekend (fh53-4...) /
#        / N times wakes up at night
#        WASI IQ test (fh62...)
# ALSO : DXA trunk fat mass (fh2245)
#        Maximum HR (fh2006)
#        N cigarettes every day TF3 (fh8451) # Only 350 observations

# MISSING : Blood markers: chol, hdl, ldl, insulin, trig (*_tf3)
#           Smoking, alcohol and cannabis 

cat('\n================================ 16 YEARS ====================================\n')

pub_16.0y <- extr_pub(age = 'pub897a') 

cmr_16.0y <- extr_cmr(data.frame('pub803'='height',  # cm (parent-reported)
                                 'pub804'='weight'), # kg (parent-reported)
                      age = 'pub897a')

# ------------------------------------------------------------------------------

cmr_16.6y <- data %>%
  # Combine smoking measurements
  transmute(smoking_16.6y = case_when(ccs4005 == 'YP smokes 1 or more cigarettes every day' ~ 6, 
                                      ccs4005 == 'YP smokes more than 6 cigarettes a week, but not every day' ~ 5, 
                                      ccs4005 == 'YP smokes between 1 and 6 cigarettes a week' ~ 4, 
                                      ccs4005 == 'YP sometimes smokes cigarettes but less than once a week' ~ 3, 
                                      ccs4005 == 'YP used to smoke cigarettes but doesn\'t now' ~ 2, 
                                      ccs4005 == 'YP has only ever smoked cigarettes once or twice' ~ 1, 
                                      ccs4000 == 'No' ~ 0,
                                      TRUE ~ NA), 
            canabis_16.6y = case_when(ccs4065 == 'YP usually uses cannabis everyday' ~ 6, 
                                      ccs4065 == 'YP uses cannabis more than 6 times a week but not every day' ~ 5, 
                                      ccs4065 == 'YP uses cannabis between 1 and 6 times a week' ~ 4, 
                                      ccs4065 == 'YP Sometimes uses cannabis but less than once a week' ~ 3, 
                                      ccs4065 == 'YP used to take cannabis but doesn\'t now' ~ 2, 
                                      ccs4065 == 'YP has tried cannabis once or twice' ~ 1, 
                                      ccs4060 == 'No' ~ 0,
                                      TRUE ~ NA),
            CMR_age_16.6y = as.numeric(as.character(ccs9991a)) /12 )
summary(cmr_16.6y)

dep_16.6y <- extr_dep('ccs45', 0:15, remove=c(2,8,11), age = 'ccs9991a')

# ALSO : PLIKS depression symptoms (ccs2660-76)
#        Parents divorce (ccs2050)
#        N cigarettes every day TF3 (fh8451) # Only 350 observations
#        Age when first smoked a cigarette (ccs4010) / Total number of cigarettes (ccs4030) /
#        / smoked any cigarettes since 15th birthday (ccs4035) / N cigarettes per day (ccs4040)
#        / nicotine patches or gums (ccs4050)
#        Age when first tried cannabis (ccs4070) / Total number of times used cannabis (ccs4075) /
#        / type of cannabis and method (ccs4080-95) / other cannabis questions... (ccs41...)
#        Other drugs: amphetamines, ecstasy, LSD ...(ccs4150-70)
# MISSING : alcohol (ccs35...)

mdep_16.7y <- extr_dep('tc40',30:42, reporter='m', age = 'tc9991a') 

# ALSO : SDQ - symptoms + scores (tc4000-7025)
#        Study teenager has tried alcohol /  cigarettes / cannabis / ecstasy / other (tc2010-4)
#        Age at which tried alcohol /  cigarettes / cannabis / ecstasy / other (tc2020-4)
# 

cat('\n================================ 17 YEARS ====================================\n')

pub_17.0y <- extr_pub(age = 'pub997a') 

cmr_17.0y <- extr_cmr(data.frame('pub903'='height',  # cm (parent-reported)
                                 'pub904'='weight'), # kg (parent-reported)
                      age = 'pub997a')

# ------------------------------------------------------------------------------

dep_17.8y <- extr_dep('ccxd9', 0:15, remove=c(2,8,11), age = 'ccxd006') # (05 in years), CCXD917 total score

# ALSO : Bachman self esteem symptoms + score (ccxd8)

to_numeric(c('fjel024','fjel116', # Pulse (seated after 5mins rest) and Heart Rate: ELBA
             'fjal4000', 'fjsm1000')) # Alcohol and smoking scores

data <- data %>%
  # Combine two pulse measures
  mutate(fj_pulse = rowMeans(.[, c('fjel024','fjel116')], na.rm=TRUE))

cmr_17.8y <- extr_cmr(data.frame('fjmr020'='height', # cm
                                 'fjmr022'='weight', # kg
                                 'fjdx135'='total_fatmass',  # grams (--> transformed later)
                                 'fjdx136'='total_leanmass', # grams (--> transformed later)
                                 'fjdx138'='android_fatmass',# grams (--> transformed later)
                                 'fjdx141'='gynoid_fatmass', # grams (--> transformed later)
                                'fjmr025a'='fat_percent', # Fat percentage (%)
                                 'fjel022'='SBP', # Systolic Pulscor blood pressure (seated after 5mins rest): ELBA
                                 'fjel023'='DBP', # Diastolic Pulscor blood pressure (seated after 5mins rest): ELBA: TF4
                                'fjar083d'='PWV', # carotid to femoral (m/s) average
                                'fj_pulse'='heart_rate', # Pulse / Heart Rate: ELBA
                             # assuming the age is the same here 
                              'serumc_tf4'='tot_chol', # Serum total cholesterol (mmol/l): TF4
                                'hdlc_tf4'='HDL_chol', # Total cholesterol in HDL (mmol/l): TF4
                                'ldlc_tf4'='LDL_chol', # Total cholesterol in LDL (mmol/l): TF4
                             'serumtg_tf4'='triglyc',  # Serum total triglycerides (mmol/l): TF4
                                 'glc_tf4'='glucose',  # Glucose (mmol/l): TF4
                                'fjal4000'='alcohol',  # Total score for the alcohol use disorders identification test (AUDIT): TF4
                                'fjsm1000'='smoking'), # FTND total score: TF4 smoking
                   age = 'fj003a') # Age in months ar clinic visit [F17] (b in years)

# ALSO : CIS-R subdomains + total score (FJCI050)
#        DXA trunk fat mass (FJDX126) and bone mass (FJDX134)
#        Carotid-radial PWV (FJAR088d)
#        Fat in liver (FJLI100) > only binary (0,1) with 43 cases
#        Parents divorce (fjle112-3)
#        Alcohol: ever had a whole drink (FJAL050) / age first had a full drink (FJAL100)
#        / Freq (FJAL1000) / N of drinks a day FJAL1050 / AUDIT categorical (FJAL4001)
#        Smoking: Number of cigarettes per day (FJSM400)

# MISSING : CRP (crp_tf4); Insulin
#           IMT (fjar079d), Heart measures (GRACE: FJGR...)
#           Cannabis and other drugs (FJDR...)

cat('\n================================ 18 YEARS ====================================\n')

dep_18.7y <- extr_dep('cct27',  0:12, age = 'cct9991a') # (c in years), cct2715 = total score

cmr_18.7y <- data %>%
            # Combine smoking, alcohol and cannabis measurements
  transmute(smoking_18.7y = case_when(cct5012 == 'Yes' ~ 3, # Every day
                                      cct5014 == 'Yes' ~ 2, # Every week
                                      cct5000 == 'Yes' ~ 1, # Not every week
                                      cct5000 == 'No' ~ 0, # Never smoked 
                                      TRUE ~ NA),
            alcohol_18.7y = case_when(cct5030 == '4 or more times a week' ~ 4,
                                      cct5030 == '2-3 times a week' ~ 3,
                                      cct5030 == '2-4 times a month' ~ 2,
                                      cct5030 == 'Monthly or less' ~ 1,
                                      cct5030 == 'Never' ~ 0,
                                      cct5020 == 'No' ~ 0, # Never a whole alcoholic drink
                                      TRUE ~ NA),
            canabis_18.7y = case_when(cct5055 == 'Daily or almost daily' ~ 5, # in past 12 months
                                      cct5055 == 'Weekly' ~ 4, 
                                      cct5055 == 'Monthly (but less than weekly)' ~ 3, 
                                      cct5055 == 'Less than monthly' ~ 2, 
                                      cct5055 == 'Once or twice' ~ 1, 
                                      cct5050 == 'No' ~ 0,
                                      TRUE ~ NA))

lapply(cmr_18.7y, table, useNA='ifany')

# ALSO : Age when smoked first (cct5001) / N cigarettes in lifetime (cct5005) /
#        / smoked in past month (cct5010) / Age when smoked last (cct5011) /
#        / N cigarettes daily (cct5013) / N cigarettes weekly (cct5015) /
#        Age when drink first (cct5025) / N units per day when drinking (cct5031) /
#        / >6 units in one occasion (cct5032) / Freq unable to stop (cct5033) /
#        / other drinking questions (cct5034-9)
#       Age when first tried cannabis (cct5051) / Last time (cct5052) / Age when last tried (cct5053)
#        / Number of joints/spliffs per day (cct5056) / other cannabis questions (cct5070-5)
#       Other drugs [cocaine, amphetamine, sedatives or sleeping pills, hallucinogens...] (cct5100-82)

# MISSING : Exercise frequency (cct4105)
#           Eating disorders (cct41..)

cat('\n================================ 19 YEARS ====================================')
cat('\n================================ 20 YEARS ====================================')
cat('\n================================ 21 YEARS ====================================\n')

dep_21.9y <- extr_dep('ypa2', 0:12, age = 'ypa9020', plus=0)

cat('\n================================ 22 YEARS ====================================\n')

dep_22.9y <- extr_dep('ypb5',  0:17, remove=c(3,8,12,15,17), age = 'ypb9992', plus=0) # YPB5180 total score 

# ALSO : Diagnosis of medical conditions (YPB1210-35, *_imputeno)

# MISSING : Level of physical activity, relative to others of a similar age (ypb2050)
#           Smoking, alcohol and cannabis

cat('\n================================ 23 YEARS ====================================\n')

dep_23.8y <- extr_dep('ypc16', 50:67, remove=c(3,8,12,15,17), age = 'ypc2650')

# ALSO : Childhood and current trauma (YPC1810-69) and Life events (YPC2150-2410)

# MISSING : smoking, eating behavior, education and employment 

cat('\n================================ 23 YEARS ====================================\n')

to_numeric(c('fkcv1131','fkcv2131', # cIMT: left and right
             'fksp1837','fkbp1032', # PWA: Heart Rate average (bpm) & Average seated pulse rate (bpm) 
             'fkms1000','fkms1052','fkms1062')) # Height, waist and hip circumference
             

data <- data %>%
         # Combine left and right mean IMT measurement
  mutate(fk_imt = rowMeans(.[, c('fkcv1131','fkcv2131')], na.rm=TRUE), 
         # Combine two pulse measures
         fk_pulse = rowMeans(.[, c('fksp1837','fkbp1032')], na.rm=TRUE), 
         # Recode anthropometry measures 
         fk_height = fkms1000 / 10, # Standing height (mm) --> cm
         fk_wc = fkms1052 / 10, # Average waist circumference (mm) --> cm
         # Calculate waist to hip ratio
         fk_whr = fk_wc / (fkms1062 / 10) # Average hip circumference (mm) --> cm
  )

cmr_24.5y <- extr_cmr(data.frame('fk_height'='height', # Height (cm): F@24
                                  'fkms1030'='weight', # Weight (kg): F@24
                                     'fk_wc'='waist_circ',# Waist circumference (cm): F@24
                                    'fk_whr'='waist_hip_ratio',  # cm
                                  'fkdx1001'='total_fatmass',   # grams (--> transformed later)
                                  'fkdx1002'='total_leanmass',  # grams (--> transformed later)
                                  'fkdx1041'='android_fatmass', # grams (--> transformed later)
                                  'fkdx1051'='gynoid_fatmass',  # grams (--> transformed later)
                                  'fkli1010'='liver_fat', # Fatty liver result (CAP value; dB/m): F@24
                                  'fkbp1030'='SBP', # Average seated systolic blood pressure (mmHg): F@24
                                  'fkbp1031'='DBP', # Average seated diastolic blood pressure (mmHg): F@24
                                    'fk_imt'='IMT', # cIMT average mean value (mm)
                                  'fkcv4200'='PWV', # Average carotid-femoral pulse wave velocity (m/s): F@24
                                  'fk_pulse'='heart_rate', # bpm
                                  'fkal1500'='alcohol', # Alcohol use disorder identification test (AUDIT) score : F@24
                                  'fksm1150'='smoking', # Fagerstrom score: F@24
                                 
                                  # assuming the age is the same here 
                                'serumc_f24'='tot_chol', # Serum total cholesterol (mmol/l): F24
                               'serumtg_f24'='triglyc',  # Serum total triglycerides (mmol/l): F24
                                   'glc_f24'='glucose',  # Glucose (mmol/l): F24
                                   'crp_f24'='CRP'), # C-Reactive Protein mg/L, Focus@24
                   
                   age = 'fkar0010') # (11 in years)

# ALSO : CIS-R diagnoses (FKDQ...)
#        Stressful experience (FKDE...)
#        DXA trunk fat mass (FKDX1031)
#        PLINKS interview (FKPL...)

# MISSING : Echocardiography and tissue Doppler imaging (FKEC...)
#           Cannabis (FKCA1025) and other drugs (FJDR...)
#           Insulin, HLD and LDL cholesterol (f24)


# ==============================================================================
# ==============================================================================
# Parental education (measured at 5.1 years)

m_edu = ifelse( !is.na(data$k6280), 0, # Mother has no educational qualification
        ifelse( !is.na(data$k6292)     # Mother has a university degree
              | !is.na(data$k6286)     # Mother is a state enrolled nurse
              | !is.na(data$k6287), 2, # Mother is a state registered nurse
        ifelse( !is.na(data$k6281)     # Mother has CSE/GCSE (D,E,F,G)
              | !is.na(data$k6282)     # Mother has O-level/GCSE (A,B,C)
              | !is.na(data$k6283)     # Mother has A-levels
              | !is.na(data$k6284), 1, # Mother has vocational qualification
        NA))) 

p_edu = ifelse( !is.na(data$k6300), 0, # Partner has no educational qualification
        ifelse( !is.na(data$k6312)     # Partner has a university degree
              | !is.na(data$k6306)     # Partner is a state enrolled nurse
              | !is.na(data$k6307), 2, # Partner is a state registered nurse
        ifelse( !is.na(data$k6301)     # Partner has CSE/GCSE (D,E,F,G)
              | !is.na(data$k6302)     # Partner has O-level/GCSE (A,B,C)
              | !is.na(data$k6303)     # Partner has A-levels
              | !is.na(data$k6304), 1, # Partner has vocational qualification
              NA)))

# NOTE: other education levels present but don't know what to do with them ...
#       k6285 / k6305: Mother / Partner has done apprenticeship
#       k6288 / k6308: Mother / Partner has City & Guilds intermediate technical qualification
#       k6289 / k6309: Mother / Partner has City & Guilds final technical qualification
#       k6290 / k6310: Mother / Partner has City & Guilds full technical qualification
#       k6291 / k6311: Mother / Partner has a teaching qualification

parent_edu <- data.frame('m_edu' = m_edu, 
                         'p_edu' = p_edu)

# Display 
cat('\n=========================== Parental education ===============================\n')

addmargins(table(parent_edu, useNA = 'ifany'))

sink()

# ==============================================================================

#cat(ls(), sep=", ")
d <- cbind(dset, parent_edu, pub_08.2y,
     mdep_09.6y, cmr_09.8y, pub_09.6y,
      dep_10.6y, cmr_10.6y, pub_10.7y,
     mdep_11.7y, cmr_11.7y, pub_11.7y,
      dep_12.8y, cmr_12.8y, 
     mdep_13.1y, cmr_13.1y, pub_13.1y,
      dep_13.8y, cmr_13.8y, 
                 cmr_14.2y, pub_14.7y,
                 cmr_15.4y, pub_15.3y,
      dep_16.6y, cmr_16.0y, pub_16.0y,
     mdep_16.7y, cmr_17.0y, pub_17.0y,
      dep_17.8y, cmr_17.8y, 
      dep_18.7y, cmr_18.7y,
      dep_21.9y, 
      dep_22.9y, 
      dep_23.8y, cmr_24.5y)

# ------------------------------------------------------------------------------
# Remove all unnecessary attributes that may corrupt the file
for (var in colnames(d)) {
  attr(d[,deparse(as.name(var))], "names") <- NULL
  attr(d[,deparse(as.name(var))], "value.labels") <- NULL
}

# Save -------------------------------------------------------------------------

crm <- cor(d[,-c(1:3, grep('time', names(d)))], method='spearman', use='pairwise.complete.obs')
# vcm <- cov(d[,-c(1:3)], method='spearman', use='pairwise.complete.obs')

saveRDS(d, 'raw_data_AFAR.rds'); write.csv(d,'raw_data_AFAR.csv')

write.csv(crm, 'corr_matrix_AFAR.csv')
# write.csv(vcm, 'varcov_matrix.csv')

dev.off()
