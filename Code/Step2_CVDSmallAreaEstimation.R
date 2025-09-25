################################################################################
# Program Name: Step2_CVDSmallAreaEstimation.R
# Program Purpose: Multi-level regression and poststratification and prediction of CVD hospitalizations
# Estimates daily propensity of CVD hospitalization for an individual by age, sex, and race
# WITH community level variables

# Created by: Katherine Burley Farr
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
  library(reshape2)
  library(stargazer)
}

## 1. Bring in data ----

 cvd_counts2<-read_delim("O:\\PRIV\\IRBData\\Medicare\\Project_Folders\\SAEproject\\RTP_CVD_ZCTA_SubgroupCombination_20092019.csv",delim = "\t")
 
 CVD_data<-cvd_counts2 %>%
    mutate(ZCTA = str_pad(ZCTA_CC,5,"0",side = "left"),
           DATE = as.Date(paste0(substr(ADMSN_DT,1,4),"-",substr(ADMSN_DT,5,6),"-",substr(ADMSN_DT,7,8)),format="%Y-%m-%d")) %>%
   filter(YEAR == 2018) %>% 
   dplyr::select(ZCTA,DATE,MALE_WHT_65_74,MALE_WHT_75_84,MALE_WHT_85Plus,MALE_BLK_65_74,MALE_BLK_75_84,MALE_BLK_85Plus,MALE_OTH_65_74,
                 MALE_OTH_75_84,MALE_OTH_85Plus,MALE_API_65_74,MALE_API_75_84,MALE_API_85Plus,MALE_HPC_65_74,MALE_HPC_75_84,MALE_HPC_85Plus,
                 MALE_NTV_65_74,MALE_NTV_75_84,MALE_NTV_85Plus,FEMALE_WHT_65_74,FEMALE_WHT_75_84,FEMALE_WHT_85Plus,FEMALE_BLK_65_74,
                 FEMALE_BLK_75_84,FEMALE_BLK_85Plus,FEMALE_OTH_65_74,FEMALE_OTH_75_84,FEMALE_OTH_85Plus,FEMALE_API_65_74,FEMALE_API_75_84,
                 FEMALE_API_85Plus,FEMALE_HPC_65_74,FEMALE_HPC_75_84,FEMALE_HPC_85Plus,FEMALE_NTV_65_74,FEMALE_NTV_75_84,FEMALE_NTV_85Plus) %>%
   pivot_longer(cols = -c(ZCTA,DATE), names_to = "Variable", values_to = "Total_CVD") %>%
   ungroup()

 rm(cvd_counts2)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   

rtp_zcta<-read_csv("O:\\PRIV\\IRBData\\Medicare\\Project_Folders\\Multivariate\\Univariate_Temperature_CVD_Models\\Input_Data\\rtp_zip_AVGTemp.csv")
names(rtp_zcta)
rtp_zcta<-as.character(unique(rtp_zcta$ZCTA.x))

total_bene2<-read_delim("O:\\PRIV\\IRBData\\Medicare\\Project_Folders\\SAEproject\\RTP_TotalBENE_ZCTA_SubgroupCombination_20092019.csv",delim = "\t")
total_bene_df<-total_bene2 %>%
  group_by(ZCTA) %>%
  filter(BASF_YR_NUM == 2018) %>%
  summarise(across(FEMALE_API:MALE_UNK_85Plus, sum)) %>%
 dplyr::select(ZCTA,FEMALE_API_65_74,FEMALE_API_75_84,FEMALE_API_85Plus,FEMALE_BLK_65_74,FEMALE_BLK_75_84,FEMALE_BLK_85Plus,FEMALE_HPC_65_74,FEMALE_HPC_75_84,
FEMALE_HPC_85Plus,FEMALE_NAT_65_74,FEMALE_NAT_75_84,FEMALE_NAT_85Plus,FEMALE_OTH_65_74,FEMALE_OTH_75_84,FEMALE_OTH_85Plus,
FEMALE_WHT_65_74,FEMALE_WHT_75_84,FEMALE_WHT_85Plus,MALE_API_65_74,MALE_API_75_84,MALE_API_85Plus,MALE_BLK_65_74,MALE_BLK_75_84,MALE_BLK_85Plus,MALE_HPC_65_74,   
MALE_HPC_75_84,MALE_HPC_85Plus,MALE_NAT_65_74,MALE_NAT_75_84,MALE_NAT_85Plus,MALE_OTH_65_74,MALE_OTH_75_84,MALE_OTH_85Plus,   
MALE_WHT_65_74,MALE_WHT_75_84,MALE_WHT_85Plus) %>%
  pivot_longer(cols = -ZCTA, names_to = "Variable", values_to = "Value") %>%
  ungroup()

table(total_bene_df$ZCTA)

rm(total_bene2)
unique(total_bene_df$Variable)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   

pop_overlap_orig <- read.table("./Data/Original/Geography/RTP_zcta10CBG20Pop.csv", sep="\t", header=T)
focus_zctas_orig <- st_read("./Data/Original/Geography/RTP_CBSAZCTA.shp")
xwalk_orig <- read_csv("./Data/Original/Geography/zipcodeZCTA2009to2020_melt.csv")

p12_all <- read_csv("./Data/Original/Census/Census_P12_Tables.csv")
acs_zcta_orig <- read_csv("./Data/Original/Census/nhgis0024_ds244_20195_zcta_E.csv") # ACS 2015-2019
acs_cbg_orig <- read_csv("./Data/Original/Census/nhgis0025_ds249_20205_blck_grp_E.csv") # ACS 2016-2020


## 2. Format Data  ----

  totals <- total_bene_df %>%
    mutate(zcta = as.numeric(ZCTA))  %>%
    mutate(sex_factor = as.factor(case_when(grepl("FEMALE", Variable) ~ "Female",
                                            grepl("MALE", Variable) ~ "Male"))) %>%
    mutate(age_factor = as.factor(case_when(grepl("65_74", Variable) ~ "(64,74]",
                                            grepl("75_84", Variable) ~ "(74,84]",
                                            grepl("85Plus", Variable) ~ "(84,110]"))) %>%
    mutate(race_factor = as.factor(case_when(grepl("WHT", Variable) ~ "White",
                                             grepl("BLK", Variable) ~ "Black",
                                             grepl("API", Variable) ~ "Asian",
                                             grepl("NTV", Variable) ~ "NANative",
                                             grepl("HPC", Variable) ~ "Hispanic",
                                             grepl("OTH",Variable) ~ "Other"))) %>%
    mutate(race_factor = factor(race_factor, levels=c("White","Black","NANative","Hispanic","Other"))) %>% # EPA confirm this works to set White as base in regression
    # Keep only combos with non-NA
    filter(!is.na(sex_factor) & !is.na(age_factor) & !is.na(race_factor)) %>% # no observations for NANative Male 85+
    select(-c(Variable,ZCTA))

#---------------------------------------------------------------------------------------------------------------------------------

CVD_daily_counts <- CVD_data %>%
  mutate(zcta = as.numeric(ZCTA))  %>%

  mutate(sex_factor = as.factor(case_when(grepl("FEMALE", Variable) ~ "Female",
                                          grepl("MALE", Variable) ~ "Male"))) %>%
  mutate(age_factor = as.factor(case_when(grepl("65_74", Variable) ~ "(64,74]",
                                          grepl("75_84", Variable) ~ "(74,84]",
                                          grepl("85Plus", Variable) ~ "(84,110]"))) %>%
  mutate(race_factor = as.factor(case_when(grepl("WHT", Variable) ~ "White",
                                           grepl("BLK", Variable) ~ "Black",
                                           grepl("API", Variable) ~ "Asian",
                                           grepl("NTV", Variable) ~ "NANative",
                                           grepl("HPC", Variable) ~ "Hispanic",
                                           grepl("OTH",Variable) ~ "Other"))) %>%
  # Keep only combos with non-NA
  filter(!is.na(sex_factor) & !is.na(age_factor) & !is.na(race_factor)) %>% # no observations for NANative Male 85+
  select(-c(Variable,ZCTA))

# Create Df with all possible combinations
age <- unique(CVD_daily_counts$age_factor)
sex <- unique(CVD_daily_counts$sex_factor)
race <- unique(CVD_daily_counts$race_factor)
zcta <- unique(CVD_daily_counts$zcta)
dates <- unique(CVD_daily_counts$DATE)
#---------------------------------------------------------------------------------------------------------------------------------

### 2.3 Format ZCTA Level Variables -----

acs_zcta <- acs_zcta_orig %>%
  select(ZCTA5A, ALWGE001, ALWGE022, ALWGE023, ALWGE024, ALWGE025, ALWGE025, ALWYE001, ALWYE002, ALWYE008, ALWYE014, ALWYE019, ALWYE025, ALWYE030, ALWYE037, ALWYE043, ALWYE048, ALWYE054, ALWYE059, ALW1E001) %>%
  # Educational Attainment for the Population 25 Years and Over, ACS 2018-2022
  mutate(pct_w_bach = (ALWGE022+ALWGE023+ALWGE024+ALWGE025)/ALWGE001) %>%
  
  # Median Household Income in the Past 12 Months (in 2022 Inflation-Adjusted Dollars), ACS 2018-2022
  rename(median_hh_inc = ALW1E001) %>%
  
  # Households with Income in Past 12 Months Below Poverty Level, ACS 2018-22 (All and Those w/ householder 65+ )
  mutate(pct_hh_below_pl = ALWYE002/ALWYE001) %>%
  mutate(pct_hh_below_pl_65 = (ALWYE008+ALWYE014+ALWYE019+ALWYE025+ALWYE030)/(ALWYE008+ALWYE014+ALWYE019+ALWYE025+ALWYE030+ALWYE037+ALWYE043+ALWYE048+ALWYE054+ALWYE059)) %>%
  
  # Format
  mutate(zcta = as.numeric(ZCTA5A)) %>%
  select(zcta, pct_w_bach, median_hh_inc, pct_hh_below_pl, pct_hh_below_pl_65)


### 2.4 Create Regression Dataset -----

# Summary of the total # of beneficiaries and total # of CVD hospitalizations per group and date
combo_grid_daily <- expand.grid(age, sex, race, zcta, dates) %>%
  rename(age_factor = Var1,
         sex_factor = Var2,
         race_factor = Var3,
         zcta = Var4,
         date = Var5) %>%
  # Remove unknown
  filter(race_factor != "Unknown") %>%
  mutate(combo_id = row_number()) %>%
  
  # Bring in total beneficiary counts by group for each day
  left_join(totals, by=c("age_factor", "sex_factor", "race_factor", "zcta")) %>%
  
  # Bring in hospitalization totals by group and date
  left_join(CVD_daily_counts, by=c("date"="DATE","age_factor", "sex_factor", "race_factor", "zcta")) %>%

  # Check that all combos with CVD have at least 1 beneficiary
  mutate(total_CVD = replace_na(Total_CVD,0)) %>%
  mutate(total_ben = ifelse(Value==0 & total_CVD>0,total_CVD,Value)) %>%
    
  # NANative, Male, 85+ missing entirely - drop out
  filter(total_ben>0) %>%
  
  # Bring in zcta-level data (education, income, poverty)
  left_join(acs_zcta, by=c("zcta")) %>%

  # FOR AGGREGATED MODEL:
  mutate(cvd_pct = total_CVD/total_ben) # this is the dependent variable

rm(acs_zcta, acs_zcta_orig, xwalk, totals, pseudo_counts_daily)

## 3. Estimate the Multi-Level Regression (MR of MRP) ----

# For efficiency reasons, model is run with aggregated day-group data, weighted by total beneficiaries
# Predictions will estimate daily probability of CVD hospitalization, with cvd_pct as the DV
# If the data was in person-day format, we would run with model with an indicator for experiencing a CVD hospitalization for a given beneficiary-day (same result)
# Could incorporate seasonal/cyclical variables here such as holiday, weekend/weekday, season, etc. 

# FINAL MODEL SPECS:
# Level 1 (Individual) Vars: age, race, sex
# Level 2 (ZCTA) Vars: pct_hh_below_pl, pct_w_bach 
# ZCTA random effect

start.time <- Sys.time()

model_agg <- glmer(cvd_pct ~ age_factor + race_factor + sex_factor + pct_hh_below_pl + pct_w_bach + (1 | zcta), 
                   data = combo_grid_daily, weights = total_ben, family = "binomial")

end.time <- Sys.time()

start.time-end.time 

# Export Model Results
saveRDS(model_agg, file="./Data/Analysis/2_model_agg.RDS")

# Explore Model Results
model_agg <- readRDS("../Data/Analysis/2_model_agg.RDS.RDS")
summary(model_agg)
stargazer(model_agg, type="html", dep.var.labels = c("CVD"),
          covariate.labels = c("Age (74-83)", "Age (84+)", 
                               "Race (Black)", "Race (Hispanic)", "Race (Other)", "Race (White)",
                               "Sex (Male)", "Households Below PL (%)", "Pop with Bachelor's Degree+ (%)"), 
          out = "../Data/Results/2_MRP_Model_Results.htm")

#OR <- exp(coef(model_agg))
OR <- coef(model_agg)
OR 

fixef(model_agg)
exp(fixef(model_agg))

# Export OR and confidence intervals to CSV
ci_wald <- confint(model_agg,parm="beta_", method="Wald")  ## slow (~ 11 seconds)
# ci_profile <- confint(model_agg,parm="beta_", method="profile")  ## slow (~ 11 seconds)
ci_wald_or <- exp(ci_wald)
# ci_profile <- exp(ci_profile)

# Export odds ratios for table
or <- data.frame(exp(fixef(model_agg))) %>%
  cbind(ci_wald_or)

write.csv(or, "../Data/Results/2_MRP_Model_Results.csv")

sjPlot::plot_model(model_agg, show.values=TRUE, show.p=TRUE,
                   title="Odds Ratios for Experiencing Daily CVD Hospitalization") +
  theme_bw()


## 4. Predict  ----
# Using the model results, we will predict the daily probability of experiencing a CVD hospitalization 
# for an individual in a given age group, sex group, race group, and CBG.

### 4.1 Identify CBGs within Focus ZCTAs ----

# 86 ZCTAs within the RTP CBSAs
focus_zctas <- focus_zctas_orig %>%
  st_drop_geometry() %>%
  mutate(zcta = as.numeric(ZCTA5)) %>%
  select(zcta) 

# Format the population overlaps file 
pop_overlap <- pop_overlap_orig %>%
  rename(zcta = ZCTA5,
         pop_weight = PCT_CBGPop) %>%
  select(GEOID, zcta, pop_weight)

focus_cbgs_pop <- pop_overlap %>% 
  filter(zcta %in% focus_zctas$zcta) %>% 
  group_by(GEOID) %>% 
  summarise(total_pct = sum(pop_weight)) %>% # 1436
  ungroup() %>%
  mutate(total_pct = round(total_pct, 10)) %>% # something weird here
  filter(total_pct==100)

### 4.2 Format CBG Level Variables ----

# Education - AMRZ
# Median Inc - AMR8
# Poverty - AMR5

acs_cbg <- acs_cbg_orig %>%
  select(GEOID, AMRZE001, AMRZE022, AMRZE023, AMRZE024, AMRZE025, AMRZE025, AMR5E001, AMR5E002, AMR5E008, AMR5E014, AMR5E019, AMR5E025, AMR5E030, AMR5E037, AMR5E043, AMR5E048, AMR5E054, AMR5E059, AMR8E001) %>%
  separate(GEOID, into=c("Geography", "GEOID"), sep="US") %>%
  
  # Educational Attainment for the Population 25 Years and Over, ACS 2016-2020
  mutate(pct_w_bach = (AMRZE022+AMRZE023+AMRZE024+AMRZE025)/AMRZE001) %>%
  
  # Median Household Income in the Past 12 Months (in 2020 Inflation-Adjusted Dollars), ACS 2016-2020
  rename(median_hh_inc = AMR8E001) %>%
  
  # Households with Income in Past 12 Months Below Poverty Level, ACS 2016-2020 (All and Those w/ householder 65+ )
  mutate(pct_hh_below_pl = AMR5E002/AMR5E001) %>%
  mutate(pct_hh_below_pl_65 = (AMR5E008+AMR5E014+AMR5E019+AMR5E025+AMR5E030)/(AMR5E008+AMR5E014+AMR5E019+AMR5E025+AMR5E030+AMR5E037+AMR5E043+AMR5E048+AMR5E054+AMR5E059)) %>%
  
  # Format
  mutate(GEOID = as.numeric(GEOID)) %>%
  select(GEOID, pct_w_bach, median_hh_inc, pct_hh_below_pl, pct_hh_below_pl_65)

### 4.2 Create Grid of All Possible Combinations for Prediction ----

# ZCTA is required for the ZCTA fixed effect
# Because CBGs don't nest perfectly w/i ZCTAs, this is handled with population weighting in this section

# CREATE PREDICTION GRID with all relevant variables
age <- unique(combo_grid_daily$age_factor)
sex <- unique(combo_grid_daily$sex_factor)
race <- unique(combo_grid_daily$race_factor)
zcta <- unique(combo_grid_daily$zcta)
cbgs <- unique(focus_cbgs_pop$GEOID) 

combo_grid_cbg <- expand.grid(age, sex, race, zcta, cbgs) %>%
  rename(age_factor = Var1,
         sex_factor = Var2,
         race_factor = Var3,
         zcta = Var4,
         GEOID = Var5) %>%
  filter(zcta != 27709) %>% # not used in model due to missing data/low pop counts

  # Bring in population weights
  left_join(pop_overlap, by=c("GEOID", "zcta")) %>%
  filter(!is.na(pop_weight)) %>%
  
  # Bring in CBG poverty and educ
  left_join(acs_cbg, by="GEOID")

# PREDICT
# Predictions using age, sex, race, ZCTA random effect, and CBG-level education and poverty
combo_grid_cbg$prediction = predict(model_agg, newdata=combo_grid_cbg, type="response")
  
# FOR CBG-LEVEL PREDICTIONS: Calculate weighted avg preds using overlap pcts
final_predictions <- combo_grid_cbg %>%
  mutate(weight = pop_weight/100) %>%
  mutate(apply_weights = weight*prediction) %>%
  group_by(GEOID, age_factor, sex_factor, race_factor, median_hh_inc, pct_hh_below_pl, pct_hh_below_pl_65, pct_w_bach) %>%
  summarise(pred_weighted = sum(apply_weights)) # 94 still missing predictions

write.csv(final_predictions, "./Data/Analysis/2_Predictions_Multilevel.csv", row.names=F)

rm(combo_grid, educ_cbg, inc_cbg, pov_cbg, model_agg)
rm(age, dates, race, sex, zcta, cbgs)

## 5. Bring in Census Data + Geographic Data for Post-Stratification (P in MRP) ----

# final_predictions = read_csv("Data/Analysis/01_Predictions_Multilevel.csv")

# Bring in data - Shapefiles for CBGs
# cbg_sf <- st_read("Data/TIGER_Line_Shapefiles/tl_2020_37_bg/tl_2020_37_bg.shp") %>%
#   filter(GEOID %in% unique(final_predictions$GEOID))

# Bring in Data - Census - P12 - Sex by Age - CBG
# See "data-table-guide-dhc-dp" for table details
# Commented out code below to show where the P12 data frame is coming from

# # Hispanic or Latino
# p12h <- read_csv("Data/Census/CBG/DECENNIALDHC2020.P12H-Data.csv", skip=1) %>%
#   mutate(census_race = "Hispanic or Latino") %>%
#   mutate(race_factor = "Hispanic")
# 
# # White Alone, Not Hispanic or Latino
# p12i <- read_csv("Data/Census/CBG/DECENNIALDHC2020.P12I-Data.csv", skip=1) %>%
#   mutate(census_race = "White Alone") %>%
#   mutate(race_factor = "White")
# 
# # Black or African American Alone, Not Hispanic or Latino
# p12j <- read_csv("Data/Census/CBG/DECENNIALDHC2020.P12J-Data.csv", skip=1) %>%
#   mutate(census_race = "Black or African American Alone") %>%
#   mutate(race_factor = "Black")
# 
# # American Indian and Alaska Native Alone, Not Hispanic or Latino
# p12k <- read_csv("Data/Census/CBG/DECENNIALDHC2020.P12K-Data.csv", skip=1) %>%
#   mutate(census_race = "American Indian and Alaska Native Alone") %>%
#   mutate(race_factor = "NANative")
# 
# # Asian Alone, Not Hispanic or Latino
# p12l <- read_csv("Data/Census/CBG/DECENNIALDHC2020.P12L-Data.csv", skip=1) %>%
#   mutate(census_race = "Asian Alone") %>%
#   mutate(race_factor = "Asian")
# 
# # Native Hawaiian and Other Pacific Islander Alone, Not Hispanic or Latino
# p12m <- read_csv("Data/Census/CBG/DECENNIALDHC2020.P12M-Data.csv", skip=1) %>%
#   mutate(census_race = "Native Hawaiian and Other Pacific Islander Alone") %>%
#   mutate(race_factor = "Other")
# 
# # Some Other Race Alone, Not Hispanic or Latino
# p12n <- read_csv("Data/Census/CBG/DECENNIALDHC2020.P12N-Data.csv", skip=1) %>%
#   mutate(census_race = "Some Other Race Alone") %>%
#   mutate(race_factor = "Other")
# 
# # Two or More Races, Not Hispanic or Latino
# p12o <- read_csv("Data/Census/CBG/DECENNIALDHC2020.P12O-Data.csv", skip=1) %>%
#   mutate(census_race = "Two or More Races") %>%
#   mutate(race_factor = "Other")
# 
# p12_all <- bind_rows(p12h, p12i, p12j, p12k, p12l, p12m, p12n, p12o) %>%
#   separate(Geography, into=c("Geography", "GEOID"), sep="US") %>%
#   filter(GEOID %in% unique(final_predictions$GEOID)) %>%
#   mutate(GEOID = as.numeric(GEOID)) %>%
#   select(`Geography`, `GEOID`, `Geographic Area Name`, `race_factor`, everything())
# 
# rm(p12h, p12i, p12j, p12k, p12l, p12m, p12n, p12o)


### 5.1 Format Census Data -----
# Using 2020 CBG population counts.

census_data <- p12_all %>%
  select(`Geography`, `GEOID`, `Geographic Area Name`, `race_factor`, `!!Total:!!Female:!!65 and 66 years`,
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
  select(-Age_census, Age, Sex) %>%
  group_by(`Geography`, `GEOID`, `Geographic Area Name`, `race_factor`, `sex_factor`, `age_factor`) %>%
  summarise(pop_count = sum(pop_count)) %>%
  ungroup()

### 6. Estimate Burdens and Rates ----

burdens <- census_data %>%
  left_join(final_predictions, by=c("age_factor", "sex_factor", "race_factor","GEOID")) %>% 
  mutate(burden = pred_weighted*pop_count) %>% # daily estimate of TOTAL CVD BURDEN is PREDICTED PROBABILITY * POPULATION
  group_by(`GEOID`, `Geographic Area Name`) %>%
  summarise(pop_count = sum(pop_count, na.rm = T),
            burden = sum(burden, na.rm = T)) %>% # this is daily
  ungroup() %>%
  mutate(daily_incidence_10000 = (burden/pop_count)*10000) %>% # pop count is for population 65+
  filter(pop_count>10) # 8 CBGs, pop_count is for people 65+

write.csv(burdens, "./Data/Analysis/2_Predicted_Incidence_Rates_ML.csv", row.names=F)








