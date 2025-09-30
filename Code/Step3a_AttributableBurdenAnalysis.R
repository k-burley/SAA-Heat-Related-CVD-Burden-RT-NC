###################################################################
## Program Name: Step3_AttributableBurdenAnalysis.R
## Program Purpose: Calculate Excess Heat Attributable CVD Hospitalization Burden
## Relies on attrdl function from Gasparrini and Leone 2014: https://github.com/gasparrini/2014_gasparrini_BMCmrm_Rcodedata

## Created: March 9, 2025
## Author: Melissa McInroe, Katherine Burley Farr based on code by Cassandra O'Lenick
###################################################################

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#               FORWARD ATTRIBUTABLE BURDEN
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#       ATTRIBUTABLE NUMBER, ATTRIBUTABLE FRACTION

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 1. SET UP  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###~~~~~~~~~~~~~~~~~~~~~~~~
###  Packages  ----
###~~~~~~~~~~~~~~~~~~~~~~~~

library(dlnm)
library(data.table)
library(splines)
library(weathermetrics)
library(humidity)
library(gnm)
library(reshape2)
library(tidyr)
library(MASS)
library(lubridate)
library(stringr)
library(mixmeta)
library(mvmeta)
library(tsModel)
library(DescTools)
library(tidyverse)
library(tis)
library(dplyr)
library(arrow)
library(sf)
library(RColorBrewer)

## CBG FILES: ## 

###~~~~~~~~~~~~~~~~~~~~~~~~
### Bring in  Data  ----
###~~~~~~~~~~~~~~~~~~~~~~~~

  # Daily Temperatures
  rtp_temp_df_cbg <- read_csv("../Data/Original/Temperature/RTP_CBG_Temp_2018.csv") %>% # CBG daily avg air temperatures during summer 2018
    select(-"...1") %>%
    # Convert fahrenheit to celsius 
    mutate(tmean_c_1km = (tmean_f_1km-32)*(5/9), # Zhang et al. 2022 1km gridded 
           tmean_c_ws_nn = (tmean_f_ws_nn-32)*(5/9), # Interpolated Weather Station - Nearest neighbor interpolation
           tmean_c_ws_tps = (tmean_f_ws_tps-32)*(5/9)) %>% # Interpolated Weather Station - Thin plate spline regression interpolation
    group_by(GEOID) %>%
    mutate(avg.temp = mean(tmean_c_1km, na.rm = T), # just use the hi-res temp for now
           range.temp = diff(range(tmean_c_1km, na.rm = T))) %>%
    ungroup()

  # Reduced Var-Cov Matrix
  load("../Data/Analysis/reduced_coef_vcov_CVD_RTP_2575.nsarglag_ZIPZCTA_6DayLag_YD.RData") 

  # Estimated Total CBG Hospitalizations
  cbg_mrp_results <- read_csv("../Data/Analysis/2_Predicted_Incidence_Rates_ML.csv") %>%
    select(GEOID, pop_count, burden) %>%
    rename(total_CVD = burden)
  
  # Compare temperatures (exposures) and hospitalization
  rtp_df_cbg <- rtp_temp_df_cbg %>%
    left_join(cbg_mrp_results, by="GEOID") %>% # m:1 merge
    # order for the crossbasis and lags
    arrange(GEOID, date)
  

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2. Forward Attributable Analysis - USING attrdl ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Modification of the loop used before for ZCTAs


  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ### Bring in Gasparrini's attrdl function  ----
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  source("attrdl.R") # From: https://github.com/gasparrini/2014_gasparrini_BMCmrm_Rcodedata

  # Documentation from attrdl.R file:
  #   - x: AN EXPOSURE VECTOR OR (ONLY FOR dir="back") A MATRIX OF LAGGED EXPOSURES
  #   - basis: THE CROSS-BASIS COMPUTED FROM x
  #   - cases: THE CASES VECTOR OR (ONLY FOR dir="forw") THE MATRIX OF FUTURE CASES
  #   - model: THE FITTED MODEL
  #   - coef, vcov: COEF AND VCOV FOR basis IF model IS NOT PROVIDED
  #   - model.link: LINK FUNCTION IF model IS NOT PROVIDED
  #   - type: EITHER "an" OR "af" FOR ATTRIBUTABLE NUMBER OR FRACTION
  #   - dir: EITHER "back" OR "forw" FOR BACKWARD OR FORWARD PERSPECTIVES
  #   - tot: IF TRUE, THE TOTAL ATTRIBUTABLE RISK IS COMPUTED
  #   - cen: THE REFERENCE VALUE USED AS COUNTERFACTUAL SCENARIO
  #   - range: THE RANGE OF EXPOSURE. IF NULL, THE WHOLE RANGE IS USED
  #   - sim: IF SIMULATION SAMPLES SHOULD BE RETURNED. ONLY FOR tot=TRUE
  #   - nsim: NUMBER OF SIMULATION SAMPLES
  
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ### Prepare for Loop  ----
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  cbg.idx=0
  
  num_cbg <- length(unique(rtp_df_cbg$GEOID))
  
  coefs<-overall.coef
  vcovs<-overall.vcov
  
  cbg.id<-unique(rtp_df_cbg$GEOID)
  # test <- cbg.id[1:5]
  
  # Set Reference Values for the calculations! 
  ref.val = 25
  range.min = 25 # quantile(rtp_df_cbg$tmean_c_1km,0.95,na.rm=T)
  
  # Empty Data frame to save CBG level results
  cbg_an_info <- data.frame(matrix(ncol =8))
  cbg_an_info <- data.frame(df <- data.frame(GEOID=double(),
                                      date=as.Date(x = integer(0), origin = "2000-01-01"),
                                      af_c_1km=double(),
                                      af_c_ws_nn=double(),
                                      af_c_ws_tps=double(),
                                      an_c_1km=double(),
                                      an_c_ws_nn=double(),
                                      an_c_ws_tps=double()))
                     
  names(cbg_an_info) <- c("GEOID","date", "af_c_1km", "af_c_ws_nn", "af_c_ws_tps", "an_c_1km", "an_c_ws_nn", "an_c_ws_tps")
  
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ### LOOP  ----
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # loop through each CBG to estimate forward AF and AN for each temperature measure
  
  for (z in 1:length(cbg.id)) {
    
    cbg.idx=cbg.idx+1
    
    print(z)
    
    q = cbg.id[z]
    print(q)
    
    # Data only for this CBG
    input.data<-subset(rtp_df_cbg, GEOID==q)
    
    # if(length(which(!is.na(input.data$tmean_c_1km)))==0) {
    #   next
    # }
    
    # Create crossbasis for each temperature
    cb_TMean.ns.lag6.1km <- crossbasis(input.data$tmean_c_1km, lag = c(0, 6), #constrain the lag effects to account for correlation across lag days
                                   argvar =  list(fun = "ns", knots=quantile(input.data$tmean_c_1km,c(0.25,0.75),na.rm=T)),
                                   arglag = list(fun="ns",knots=logknots(c(0,6),df=5))) # 3 internal evenly spaced knots across lag days 0-6
    cb_TMean.ns.lag6.wsnn <- crossbasis(input.data$tmean_c_ws_nn, lag = c(0, 6),
                                       argvar =  list(fun = "ns", knots=quantile(input.data$tmean_c_1km,c(0.25,0.75),na.rm=T)),
                                       arglag = list(fun="ns",knots=logknots(c(0,6),df=5))) 
    cb_TMean.ns.lag6.wstps <- crossbasis(input.data$tmean_c_ws_tps, lag = c(0, 6),
                                       argvar =  list(fun = "ns", knots=quantile(input.data$tmean_c_1km,c(0.25,0.75),na.rm=T)),
                                       arglag = list(fun="ns",knots=logknots(c(0,6),df=5))) 
    
    
    # Estimate 50th percentile for each temperature source, each CBG
    cen.val.50.1km <- quantile(input.data$tmean_c_1km,0.50,na.rm=T)
    cen.val.50.wsnn <- quantile(input.data$tmean_c_ws_nn,0.50,na.rm=T)
    cen.val.50.wstps <- quantile(input.data$tmean_c_ws_tps,0.50,na.rm=T)
    
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
      mutate(af_c_ws_nn = attrdl(x=input.data$tmean_c_ws_nn, #I thought x was only required for backward..
                               basis=cb_TMean.ns.lag6.wsnn,
                               cases=input.data$total_CVD,
                               model.link="log",
                               coef=coefs, 
                               vcov=vcovs, 
                               tot=F, # need this to get daily! 
                               type="af",
                               dir="forw",
                               cen=ref.val,
                               range=c(range.min,  max(input.data$tmean_c_ws_nn,na.rm=T)))) %>% 
      mutate(af_c_ws_tps = attrdl(x=input.data$tmean_c_ws_tps, #I thought x was only required for backward..
                               basis=cb_TMean.ns.lag6.wstps,
                               cases=input.data$total_CVD,
                               model.link="log",
                               coef=coefs, 
                               vcov=vcovs, 
                               tot=F, # need this to get daily! 
                               type="af",
                               dir="forw",
                               cen=ref.val,
                               range=c(range.min,  max(input.data$tmean_c_ws_tps,na.rm=T)))) %>% 
      
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
                               range=c(range.min,  max(input.data$tmean_c_1km,na.rm=T)))) %>% 
      mutate(an_c_ws_nn = attrdl(x=input.data$tmean_c_ws_nn, #I thought x was only required for backward..
                                 basis=cb_TMean.ns.lag6.wsnn,
                                 cases=input.data$total_CVD,
                                 model.link="log",
                                 coef=coefs, 
                                 vcov=vcovs, 
                                 tot=F, # need this to get daily! 
                                 type="an",
                                 dir="forw",
                                 cen=ref.val,
                                 range=c(range.min,  max(input.data$tmean_c_ws_nn,na.rm=T)))) %>% 
      mutate(an_c_ws_tps = attrdl(x=input.data$tmean_c_ws_tps, #I thought x was only required for backward..
                                  basis=cb_TMean.ns.lag6.wstps,
                                  cases=input.data$total_CVD,
                                  model.link="log",
                                  coef=coefs, 
                                  vcov=vcovs, 
                                  tot=F, # need this to get daily! 
                                  type="an",
                                  dir="forw",
                                  cen=ref.val,
                                  range=c(range.min,  max(input.data$tmean_c_ws_tps,na.rm=T))))
    
    estimates_final <- estimates %>%
      select("GEOID","date", "af_c_1km", "af_c_ws_nn", "af_c_ws_tps", "an_c_1km", "an_c_ws_nn", "an_c_ws_tps")
      
    cbg_an_info <- cbg_an_info %>%
      rbind(estimates_final)
    
    rm(cb_TMean.ns.lag6.1km, cb_TMean.ns.lag6.wsnn, cb_TMean.ns.lag6.wstps,
       cen.val.50.1km, cen.val.50.wsnn, cen.val.50.wstps,
       input.data, estimates, estimates_final)
  
  } # end of LOOP
  
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 3. Clean Up for Export ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  
  cbg_an_final <- cbg_an_info %>%
    mutate(date = as.Date(date)) %>%
    filter(!is.na(GEOID)) %>%
    left_join(rtp_df_cbg, by=c("GEOID", "date")) %>%
    mutate(annual_an_rate10k_1km = (an_c_1km)/pop_count*10000,
           annual_an_rate10k_ws_nn = (an_c_ws_nn)/pop_count*10000,
           annual_an_rate10k_ws_tps = (an_c_ws_tps)/pop_count*10000) %>%
    mutate(reference_temp_c = ref.val)
  

  # Summer totals with 3 datasets - using 1km summer median as center
  sum(cbg_an_final$an_c_1km, na.rm=T)
  sum(cbg_an_final$an_c_ws_nn, na.rm=T)
  sum(cbg_an_final$an_c_ws_tps, na.rm=T)
  
  # Temps are slightly higher on avg in 
  summary(cbg_an_final$tmean_c_1km)
  summary(cbg_an_final$tmean_c_ws_nn)
  summary(cbg_an_final$tmean_c_ws_tps)
  
  ggplot(cbg_an_final) +
    geom_point(aes(x=tmean_c_1km, y=tmean_c_ws_nn)) +
    geom_smooth(aes(x=tmean_c_1km, y=tmean_c_ws_nn), method="lm")
  
  ggplot(cbg_an_final) +
    geom_point(aes(x=tmean_c_1km, y=tmean_c_ws_tps)) +
    geom_smooth(aes(x=tmean_c_1km, y=tmean_c_ws_tps), method="lm")
  
  write.csv(cbg_an_final, "../Data/Analysis/3_RTP_XCESS_ATTR_Burden_CBG_Summer2018_25C.csv")

  