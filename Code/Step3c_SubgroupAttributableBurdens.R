################################################################################
# Program Name: Step3c_Subgroup_Attributable_Burdens.R
# Program Purpose: Use MRP results to estimate CBG-subgroup level burdens
# NOTE: modification of code from Step3_AttributableBurdenAnalysis

# Author: Katherine Burley Farr
# Contact: kburley@ad.unc.edu
# Affiliation: UNC Department of Public Policy, Data-Driven EnviroLab
################################################################################

rm(list = ls())

{
  library(tidyverse)
  library(dplyr)
  library(readxl)
  library(writexl)
  library(openxlsx)
  library(tidytext)
  library(igraph)
  library(ggraph)
  library(tokenizers)
  library(ggrepel)
  library(ggpubr)
  library(scales)
  library(RColorBrewer)
  library(gridExtra)
  library(lubridate)
  library(sf)
  library(lme4)
  library(dlnm)
  library(data.table)
  library(splines)
  library(wesanderson)
}
library(colorspace)

## 1. Bring in data ----

final_predictions <- read_csv("../Data/Analysis/2_Predictions_Multilevel.csv")
p12_all <- read_csv("../Data/Original/Census/Census_P12_Tables.csv")

census_data <- p12_all %>%
  dplyr::select(`Geography`, `GEOID`, `Geographic Area Name`, `race_factor`, `!!Total:!!Female:!!65 and 66 years`,
         `!!Total:!!Female:!!67 to 69 years`, `!!Total:!!Female:!!70 to 74 years`, `!!Total:!!Female:!!75 to 79 years`,
         `!!Total:!!Female:!!80 to 84 years`, `!!Total:!!Female:!!85 years and over`, `!!Total:!!Male:!!65 and 66 years`,
         `!!Total:!!Male:!!67 to 69 years`, `!!Total:!!Male:!!70 to 74 years`, `!!Total:!!Male:!!75 to 79 years`,
         `!!Total:!!Male:!!80 to 84 years`, `!!Total:!!Male:!!85 years and over`) %>%
  group_by(`Geography`, `GEOID`, `Geographic Area Name`, `race_factor`) %>%
  summarise(Female.65_66 = sum(`!!Total:!!Female:!!65 and 66 years`),
            Female.67_69 = sum(`!!Total:!!Female:!!67 to 69 years`),
            Female.70_74 = sum(`!!Total:!!Female:!!70 to 74 years`),
            Female.75_79 = sum(`!!Total:!!Female:!!75 to 79 years`),
            Female.80_84 = sum(`!!Total:!!Female:!!80 to 84 years`),
            Female.85_114 = sum(`!!Total:!!Female:!!85 years and over`),
            Male.65_66 = sum(`!!Total:!!Male:!!65 and 66 years`),
            Male.67_69 = sum(`!!Total:!!Male:!!67 to 69 years`),
            Male.70_74 = sum(`!!Total:!!Male:!!70 to 74 years`),
            Male.75_79 = sum(`!!Total:!!Male:!!75 to 79 years`),
            Male.80_84 = sum(`!!Total:!!Male:!!80 to 84 years`),
            Male.85_114 = sum(`!!Total:!!Male:!!85 years and over`)) %>%
  ungroup() %>%
  pivot_longer(Female.65_66:Male.85_114, names_to="sex_age", values_to = "pop_count") %>%
  separate(sex_age, into=c("Sex","Age_census"), sep="\\.") %>%
  mutate(Age = case_when(Age_census %in% c("65_66","67_69","70_74") ~ "(64,74]",
                         Age_census %in% c("75_79", "80_84") ~ "(74,84]",
                         Age_census == "85_114" ~ "(84,110]")) %>%
  mutate(age_factor = as.factor(Age),
         sex_factor = as.factor(Sex)) %>%
  dplyr::select(-Age_census, Age, Sex) %>%
  group_by(`Geography`, `GEOID`, `Geographic Area Name`, `race_factor`, `sex_factor`, `age_factor`) %>%
  summarise(pop_count = sum(pop_count)) %>%
  ungroup()

rm(p12_all)

### 2. Estimate Daily Attributable Number by CBG and Subgroup ----

subgroup_total_CBG <- census_data %>%
  left_join(final_predictions, by=c("age_factor", "sex_factor", "race_factor","GEOID")) %>% 
  mutate(total_CVD = pred_weighted*pop_count) %>% # daily estimate of TOTAL CVD BURDEN is PREDICTED PROBABILITY * POPULATION
  filter(pop_count>0) %>%
  filter(!is.na(total_CVD))

# Sum by Subgroup:
gender <- subgroup_total_CBG %>%
  group_by(`GEOID`, `sex_factor`) %>%
  summarise(pop_count = sum(pop_count),
            total_CVD = sum(total_CVD)) %>%
  ungroup() %>%
  mutate(subgroup_type = "sex") %>%
  rename(subgroup_value = sex_factor)

race <- subgroup_total_CBG %>%
  group_by(`GEOID`, `race_factor`) %>%
  summarise(pop_count = sum(pop_count),
            total_CVD = sum(total_CVD)) %>%
  ungroup() %>%
  mutate(subgroup_type = "race") %>%
  rename(subgroup_value = race_factor)

age <- subgroup_total_CBG %>%
  group_by(`GEOID`, `age_factor`) %>%
  summarise(pop_count = sum(pop_count),
            total_CVD = sum(total_CVD)) %>%
  ungroup() %>%
  mutate(subgroup_type = "age") %>%
  rename(subgroup_value = age_factor)

# Append - this is what will be used in attributable burdens calculation
subgroup_totals <- gender %>%
  bind_rows(race) %>%
  bind_rows(age) 

rm(subgroup_total_CBG, age, race, gender, final_predictions)
  
### 3. Estimate Heat-Attributable Burdens and Rates ----

# Daily Temperatures
rtp_temp_df_cbg <- read_csv("../Data/Original/Temperature/2_RTP_CBG_Temp_2018.csv") %>% # CBG daily temperatures during summer
  dplyr::select(-"...1") %>%
  # Convert fahrenheit to celsius 
  mutate(tmean_c_1km = (tmean_f_1km-32)*(5/9),
         tmean_c_ws_nn = (tmean_f_ws_nn-32)*(5/9),
         tmean_c_ws_tps = (tmean_f_ws_tps-32)*(5/9)) %>%
  group_by(GEOID) %>%
  mutate(avg.temp = mean(tmean_c_1km, na.rm = T), # just use the hi-res temp for now
         range.temp = diff(range(tmean_c_1km, na.rm = T))) %>%
  ungroup()

# Reduced Var-Cov Matrix
load("../Data/Analysis/reduced_coef_vcov_CVD_RTP_2575.nsarglag_ZIPZCTA_6DayLag_YD.RData") 

# Compare temperatures (exposures) and hospitalization
rtp_df_cbg_subgroup <- expand.grid(GEOID = unique(subgroup_totals$GEOID),
                                   date = unique(rtp_temp_df_cbg$date),
                                   subgroup_value = unique(subgroup_totals$subgroup_value)) %>%
  left_join(rtp_temp_df_cbg[c("GEOID", "date", "tmean_c_1km")], by=c("GEOID", "date")) %>%
  left_join(subgroup_totals, by=c("GEOID", "subgroup_value")) %>%
  filter(!is.na(pop_count)) # drop out any CBG-subgroups with no 

# Gasparrini attrdl function:
source("attrdl.R") 
  
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ### Prepare for Loop  ----
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # Values to Loop Over
  cbg.idx=0
  num_cbg <- length(unique(rtp_df_cbg_subgroup$GEOID))
  cbg.id <- unique(rtp_df_cbg_subgroup$GEOID)
  subgroup.id <- unique(rtp_df_cbg_subgroup$subgroup_value)
  
  # Model coefficients needed for calculations
  coefs<-overall.coef
  vcovs<-overall.vcov

  # Set Reference Values for the calculations!
  ref.val = 25
  range.min = 25 
  
  # Empty Data frame to save CBG level results
  cbg_an_info_sg <- data.frame(GEOID=double(),
                               date=as.Date(x = integer(0), origin = "2000-01-01"),
                               subgroup_type=as.character(),
                               subgroup_value=as.character(),
                               tmean_c_1km=double(),
                               pop_count=double(),
                               total_CVD=double(),
                               af_c_1km=double(),
                               an_c_1km=double())
  
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ### LOOP  ----
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # loop through each CBG to estimate forward AF and AN for each temperature measure
  
  # Loop through CBGS
  for (z in 1:length(cbg.id)) {
    
    cbg.idx=cbg.idx+1
    q = cbg.id[z]
    print(q)
    
    # Loop through subgroups
    for (s in subgroup.id) {
      print(s)
    
    # Data only for this CBG and subgroup
    # input.data<-subset(rtp_df_cbg, GEOID==q)
      input.data <- rtp_df_cbg_subgroup %>% filter(GEOID == q & subgroup_value == s)
      
      if(length(which(!is.na(input.data$tmean_c_1km)))==0) {
        next
      }
    
    # Create crossbasis for each temperature
    cb_TMean.ns.lag6.1km <- crossbasis(input.data$tmean_c_1km, lag = c(0, 6), #constrain the lag effects to account for correlation across lag days
                                       argvar =  list(fun = "ns", knots=quantile(input.data$tmean_c_1km,c(0.25,0.75),na.rm=T)),
                                       arglag = list(fun="ns",knots=logknots(c(0,6),df=5))) # 3 internal evenly spaced knots across lag days 0-6

    
    # Calculate DAILY Forward Attributable Fraction and Attributable Number using attrdl
    
    # FORMAT:
    #   - x: AN EXPOSURE VECTOR OR (ONLY FOR dir="back") A MATRIX OF LAGGED EXPOSURES
    #     input.data$tmean_c_[temp source] if only for backwards, not sure why this is required
    
    #   - basis: THE CROSS-BASIS COMPUTED FROM x
    #     cb_TMean.ns.lag6.[temp source]
    
    #   - cases: THE CASES VECTOR OR (ONLY FOR dir="forw") THE MATRIX OF FUTURE CASES
    #     input.data$total_CVD
    
    #   - model: THE FITTED MODEL
    #     NOT USING
    
    #   - coef, vcov: COEF AND VCOV FOR basis IF model IS NOT PROVIDED
    #     reduced coefs and vcovs provided by EPA/Melissa from RTP_CVD_WARM_ZCTA_dlnmModel.R
    
    
    #   - model.link: LINK FUNCTION IF model IS NOT PROVIDED
    #     "log" - based on previous attr_burdens code
    
    
    #   - type: EITHER "an" OR "af" FOR ATTRIBUTABLE NUMBER OR FRACTION
    #     will estimate both
    
    #   - dir: EITHER "back" OR "forw" FOR BACKWARD OR FORWARD PERSPECTIVES
    #     using forward here only
    
    
    #   - tot: IF TRUE, THE TOTAL ATTRIBUTABLE RISK IS COMPUTED
    #     for DAILY estimates, must specific "tot=F"!!!!
    
    
    #   - cen: THE REFERENCE VALUE USED AS COUNTERFACTUAL SCENARIO
    #     cen.val.50.[temp source]
    
    
    #   - range: THE RANGE OF EXPOSURE. IF NULL, THE WHOLE RANGE IS USED
    #     from 50th percentile to max temp value...
    #     I think this one just means it will only estimate for temperatures in this range?
    #     IF we want to also calculate negative burdens, will need to adjust to be min to max
    
    
    #   - sim: IF SIMULATION SAMPLES SHOULD BE RETURNED. ONLY FOR tot=TRUE
    #   - nsim: NUMBER OF SIMULATION SAMPLES
    #     since tot=F, will not be able to use
    
    
    estimates <- input.data %>%
      # ATTRIBUTABLE FRACTIONS
      mutate(af_c_1km = attrdl(x=input.data$tmean_c_1km, #I thought x was only required for backward..
                               basis=cb_TMean.ns.lag6.1km,
                               cases=input.data$total_CVD,
                               model.link="log",
                               coef=coefs, 
                               vcov=vcovs, 
                               tot=F, # need this to get daily! 
                               type="af",
                               dir="forw",
                               cen=ref.val, # cen.val.50.1km
                               range=c(range.min,  max(input.data$tmean_c_1km,na.rm=T)))) %>%
      # ATTRIBUTABLE NUMBERS
      mutate(an_c_1km = attrdl(x=input.data$tmean_c_1km, #I thought x was only required for backward..
                               basis=cb_TMean.ns.lag6.1km,
                               cases=input.data$total_CVD,
                               model.link="log",
                               coef=coefs, 
                               vcov=vcovs, 
                               tot=F, # need this to get daily! 
                               type="an",
                               dir="forw",
                               cen=ref.val,
                               range=c(range.min, max(input.data$tmean_c_1km,na.rm=T)))) 
  
    
    cbg_an_info_sg <- cbg_an_info_sg %>%
      rbind(estimates)
    
    rm(cb_TMean.ns.lag6.1km, input.data, estimates)
    } # end of subgroup loop
  } # end of LOOP
  
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### 4. Plot Summer Subgroup Burdens and Rates  ----
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  subgroup_burdens <- cbg_an_info_sg %>% 
    mutate(month = month(date)) %>%
    filter(month %in% c(5,6,7,8,9)) %>%
    
    # CBG level summer totals
    group_by(GEOID, subgroup_type, subgroup_value, pop_count) %>% 
    summarise(tmean_c_1km = mean(tmean_c_1km),
              total_CVD = sum(total_CVD),
              an_c_1km = sum(an_c_1km)) %>%
    ungroup() %>%
    
    # Summer Totals by Subgroup!
    group_by(subgroup_type, subgroup_value) %>%
    summarise(pop_count = sum(pop_count),
              tot_cvd = sum(total_CVD),
              tot_an = sum(an_c_1km)) %>%
    ungroup() %>%
    
    # Calculate Rates
    mutate(rate_an = (tot_an/pop_count)*100000) %>% # per 100,000
    mutate(rate_cvd = (tot_cvd/pop_count)*100000) %>% # per 100,000
    
    # What about ratios?
    mutate(pct = (tot_an/tot_cvd)*100) # the same
  
  # EXPORT:
  write.csv(subgroup_burdens, "../Data/Results/3c_Subgroup_Attributable_Burdens.csv", row.names=FALSE)
