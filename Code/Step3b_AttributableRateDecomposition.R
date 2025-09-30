################################################################################
# Program Name: Step3b_AttributableRateDecomposition.R
# Program Purpose: Decompose the attributable burdens
# NOTE: modification of code from Step3_AttributableBurdenAnalysis.R

# Author: Katherine Burley Farr
# Contact: kburley@ad.unc.edu
# Affiliation: UNC Department of Public Policy, Data-Driven EnviroLab
################################################################################

library(dlnm)
library(data.table)
library(splines)
library(weathermetrics)
library(humidity)
library(gnm)
library(reshape2)
library(tidyr)
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
library(waterfalls)

###~~~~~~~~~~~~~~~~~~~~~~~~
### Bring in  Data  ----
###~~~~~~~~~~~~~~~~~~~~~~~~

# Daily Temperatures
rtp_temp_df_cbg <- read_csv("../Data/Original/Temperature/RTP_CBG_Temp_2018.csv") %>% # CBG daily temperatures during summer
  select(-"...1") %>%
  # Convert fahrenheit to celsius - only using 1km for decomposition
  mutate(tmean_c_1km = (tmean_f_1km-32)*(5/9)) %>%
  group_by(GEOID) %>%
  ungroup()

# Reduced Var-Cov Matrix
load("../Data/Analysis/reduced_coef_vcov_CVD_RTP_2575.nsarglag_ZIPZCTA_6DayLag_YD.RData") 

# Estimated Total CBG Hospitalizations
cbg_mrp_results <- read_csv("../Data/Analysis/2_Predicted_Incidence_Rates_ML.csv") %>% 
  select(GEOID, pop_count, burden) %>%
  rename(total_CVD = burden) %>%
  # SPECIFICALLY FOR DECOMPOSITION
  mutate(total_cvd_pc = total_CVD/pop_count) %>% # back out CVD incidence per capita per CBG
  mutate(total_cvd_pc_avg = mean(total_cvd_pc)) %>%  # get the average CVD incidence per capita across all CBGs 
  mutate(total_cvd_avg = total_cvd_pc_avg*pop_count) %>% # multiply by CBG-level population 
  select(-c(total_cvd_pc, total_cvd_pc_avg))
  

# Compare temperatures (exposures) and hospitalization
rtp_df_cbg <- rtp_temp_df_cbg %>%
  left_join(cbg_mrp_results, by="GEOID") %>% # m:1 merge
  # order for the crossbasis and lags
  arrange(GEOID, date) %>%
  # SPECIFICALLY FOR DECOMPOSITION: 
  select(-c(tmean_f_ws_nn, tmean_f_ws_tps)) %>%
  group_by(date) %>%
  mutate(tmean_c_avg = mean(tmean_c_1km)) %>% # average daily mean temperature across CBGs
  filter(!is.na(total_CVD)) %>% # using these for the avg. 
  ungroup() 
  
source("attrdl.R") # From: https://github.com/gasparrini/2014_gasparrini_BMCmrm_Rcodedata

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2. Forward Attributable Analysis - USING attrdl - KB ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Modification of the loop used before for ZCTAs

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Prepare for Loop  ----
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cbg.idx=0

num_cbg <- length(unique(rtp_df_cbg$GEOID))

coefs<-overall.coef
vcovs<-overall.vcov

cbg.id<-unique(rtp_df_cbg$GEOID)

# Set Reference Values for the calculations! 
ref.val = 25
range.min = 25 # quantile(rtp_df_cbg$tmean_c_1km,0.95,na.rm=T)

# Empty Data frame to save CBG level results
cbg_an_info <- data.frame(matrix(ncol =6))
cbg_an_info <- data.frame(df <- data.frame(GEOID=double(),
                                           date=as.Date(x = integer(0), origin = "2000-01-01"),
                                           an_atac=double(),
                                           an_atgc=double(),
                                           an_gtac=double(),
                                           an_gtgc=double()))

names(cbg_an_info) <- c("GEOID","date", "an_atac", "an_atgc", "an_gtac", "an_gtgc")

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
  
  # Create crossbasis for each temperature
  cb_TMean.ns.lag6.dailyavg <- crossbasis(input.data$tmean_c_avg, lag = c(0, 6), #constrain the lag effects to account for correlation across lag days
                                     argvar =  list(fun = "ns", knots=quantile(input.data$tmean_c_avg,c(0.25,0.75),na.rm=T)),
                                     arglag = list(fun="ns",knots=logknots(c(0,6),df=5))) # 3 internal evenly spaced knots across lag days 0-6
  
  cb_TMean.ns.lag6.1km <- crossbasis(input.data$tmean_c_1km, lag = c(0, 6), #constrain the lag effects to account for correlation across lag days
                                     argvar =  list(fun = "ns", knots=quantile(input.data$tmean_c_1km,c(0.25,0.75),na.rm=T)),
                                     arglag = list(fun="ns",knots=logknots(c(0,6),df=5))) # 3 internal evenly spaced knots across lag days 0-6
  
  # Calculate DAILY Forward Attributable Fraction and Attributable Number using attrdl

  # COMBINATION CODES: 
  # atac = avg. temperature, avg. total CVD
  # gtgc = cbg temperature, cbg total CVD
  # atgc = avg. temperature, cbg total CVD
  # gtac = cbg temperature, avg. total CVD
  
  estimates <- input.data %>%

    # ATTRIBUTABLE NUMBERS
    # average temperature, average total CVD
    mutate(an_atac = attrdl(x=input.data$tmean_c_avg,
                                basis=cb_TMean.ns.lag6.dailyavg, # use CB with daily avg temperatures for every CBG
                                cases=input.data$total_cvd_avg,
                                model.link="log",
                                coef=coefs, 
                                vcov=vcovs, 
                                tot=F, # need this to get daily! 
                                type="an",
                                dir="forw",
                                cen=ref.val,
                                range=c(range.min,  max(input.data$tmean_c_avg,na.rm=T))))  %>%
    # average temperature, CBG-specific total CVD
    mutate(an_atgc = attrdl(x=input.data$tmean_c_avg,
                                basis=cb_TMean.ns.lag6.dailyavg, # use CB with daily avg temperatures for every CBG
                                cases=input.data$total_CVD,
                                model.link="log",
                                coef=coefs, 
                                vcov=vcovs, 
                                tot=F, # need this to get daily! 
                                type="an",
                                dir="forw",
                                cen=ref.val,
                                range=c(range.min,  max(input.data$tmean_c_avg,na.rm=T))))  %>%
    # CBG-specific temperature, avg total CVD
    mutate(an_gtac = attrdl(x=input.data$tmean_c_1km,
                                basis=cb_TMean.ns.lag6.1km, # use CB with daily avg temperatures for every CBG
                                cases=input.data$total_cvd_avg,
                                model.link="log",
                                coef=coefs, 
                                vcov=vcovs, 
                                tot=F, # need this to get daily! 
                                type="an",
                                dir="forw",
                                cen=ref.val,
                                range=c(range.min,  max(input.data$tmean_c_1km,na.rm=T))))  %>%
    # CBG-specific temperature, CBG-specific total CVD (normal calculation)
    mutate(an_gtgc = attrdl(x=input.data$tmean_c_1km, 
                             basis=cb_TMean.ns.lag6.1km,
                             cases=input.data$total_CVD,
                             model.link="log",
                             coef=coefs, 
                             vcov=vcovs, 
                             tot=F, # need this to get daily! 
                             type="an",
                             dir="forw",
                             cen=ref.val,
                             range=c(range.min,  max(input.data$tmean_c_1km,na.rm=T)))) 

  
  estimates_final <- estimates %>%
    select("GEOID","date", "an_atac", "an_atgc", "an_gtac", "an_gtgc")
  
  cbg_an_info <- cbg_an_info %>%
    rbind(estimates_final)
  
} # end of LOOP

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 3. Clean Up for Export ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

cbg_an_final <- cbg_an_info %>%
  mutate(date = as.Date(date)) %>%
  filter(!is.na(GEOID)) %>%
  left_join(rtp_df_cbg, by=c("GEOID", "date")) %>%
  mutate(reference_temp_c = ref.val) %>%
  mutate(an_cvddiff1 = (an_atgc - an_atac),
         an_tempdiff2 = (an_gtgc - an_atgc))

# Totals by CBG
sum_an_final_bycbg <- cbg_an_final %>%
  mutate(month = month(date)) %>%
  filter(month %in% c(5,6,7,8,9)) %>%
  group_by(GEOID, pop_count) %>%
  summarise(an_atac = sum(an_atac), 
            an_atgc = sum(an_atgc),
            an_gtac = sum(an_gtac),
            an_gtgc = sum(an_gtgc),,
            an_cvddiff1 = sum(an_cvddiff1),
            an_tempdiff2 = sum(an_tempdiff2)) %>%
  ungroup() %>%
  # Calculate Rates
  mutate(an_atac_rate = (an_atac/pop_count)*10000,
         an_gtac_rate = (an_gtac/pop_count)*10000,
         an_atgc_rate = (an_atgc/pop_count)*10000,
         an_gtgc_rate = (an_gtgc/pop_count)*10000,
         an_cvddiff1_rate = (an_cvddiff1/pop_count)*10000,
         an_tempdiff2_rate = (an_tempdiff2/pop_count)*10000)

write.csv(sum_an_final_bycbg, "../Data/Results/3b_Attributable_Rate_Decomposition.csv", row.names=F)
