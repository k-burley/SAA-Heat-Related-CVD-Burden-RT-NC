################################################################################
# Program Name: 01b_MRP_Model_Validation.R
# Program Purpose: Internal validation for the results of the Step 2
# Compare aggregated, model-based ZCTA level estimates vs. direct ZCTA-level data from Medicare
# Following Zhang et al. 2014:
  # Measures: Desc stats, Pearson corr, MSE, MAD

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
  library(reshape2)
}

# 1. Read in Data -----

# Direct, Annual ZCTA Level CDV Hospitalization Totals from EPA
# Using data below as a placeholder!
zcta_annual_direct <- read.delim("O:\\PRIV\\IRBData\\Medicare\\Project_Folders\\SAEproject\\RTP_CVDIND_Year2018.csv", header = TRUE, sep = "\t") %>%
  rename("total_CVD_direct" = "TOTAL") %>%
  select(ZCTA, total_CVD_direct)

# Step 2 Results
cbg_mrp_results <- read.csv("../Data/Analysis/2_Predicted_Incidence_Rates_ML.csv") %>%
  mutate(total_CVD_model_cbg = burden*365) %>% # convert daily total to annual
  select(GEOID, pop_count, total_CVD_model_cbg)

# CBG to ZCTA Weight Options
focus_zctas <- st_read("../Data/Original/Geography/RTP_CBSAZCTA.shp")

zcta_cbg_pop <- read.table("../Data/Original/Geography/RTP_zcta10CBG20Pop.csv", sep="\t", header=T)

check_pop <- zcta_cbg_pop %>%
  filter(GEOID %in% unique(cbg_mrp_results$GEOID)) %>% # only ~40 ZCTAs fully covered
  group_by(ZCTA5) %>%
  summarise(total_pct = sum(PCT_ZCTAPop)) %>%
  filter(total_pct > 99.9) # only 41 ZCTAs fully covered

  # Note: we only predict burdens for CBGs fully covered by the focus ZCTAs, so not all ZCTAs are fully covered by these CBGs
  # thus there is a limited set of ZCTAs we can validate on using aggregation

# ASSIGN CBGs to ZCTA rather than weighting
assigned <- zcta_cbg_pop %>%
  filter(GEOID %in% unique(cbg_mrp_results$GEOID)) %>% 
  group_by(GEOID) %>%
  mutate(max_overlap = max(PCT_CBGPop)) %>%
  ungroup() %>%
  filter(PCT_CBGPop == max_overlap)

# 2. Combine direct CVD hospitalization counts and modeled estimates  -----

# USING CBG-ZCTA ASSIGNMENT VERSION:

zcta_annual_model_assn <- assigned %>%
  filter(ZCTA5 %in% unique(focus_zctas$ZCTA5)) %>%
  
  # Bring in CBG level estimated total CVD hospitalizations
  left_join(cbg_mrp_results, by=c("GEOID")) %>%
  
  # Aggregate to ZCTA
  group_by(ZCTA5) %>%
  summarise(total_CVD_model = sum(total_CVD_model_cbg)) %>%
  ungroup() %>%
  rename(zcta = ZCTA5)
  
  
# USING POPULATION WEIGHTED VERSION:

zcta_annual_model_pop <- zcta_cbg_pop %>%
  filter(ZCTA5 %in% unique(focus_zctas$ZCTA5)) %>% # to account for diff b/w pseudo data and actual data - shouldn't make a diff with the real data

  # Bring in the CBG level burdens
  left_join(cbg_mrp_results, by=c("GEOID")) %>% 
  filter(!is.na(total_CVD_model_cbg)) %>%

  # Aggregate to ZCTA 
  mutate(total_CVD_model_cbg_wgt = total_CVD_model_cbg*(PCT_CBGPop/100)) %>% # 1. allocate cbg estimates to the zcta-cbg overlaps with PCT_CBGPop
  group_by(ZCTA5) %>%
  summarise(total_CVD_model_wgt = sum(total_CVD_model_cbg_wgt), # 2. sum by ZCTA (do not need to apply PCT_ZCTAPop as weight, it is accounted for)
            total_PCT_ZCTAPop = sum(PCT_ZCTAPop)) %>%
  ungroup() %>%
  rename(zcta = ZCTA5) %>%
  filter(total_PCT_ZCTAPop>=99.9)

# Get ZCTA population
zcta_pop_2010 <- read_csv("./Data/Original/Census/ZCTA_SexbyAge_2010.csv") %>%
  select(ZCTA5A, H76020:H76025, H76044:H76049) %>%
  filter(ZCTA5A %in% unique(focus_zctas$ZCTA5)) %>%
  rowwise() %>%
  mutate(pop_count_zcta10 = sum(c_across(H76020:H76049))) %>%
  mutate(zcta = as.integer(ZCTA5A)) %>%
  select(zcta, pop_count_zcta10) 

# Create DF with ZCTA-level comparisons of aggregate modeled and direct estimates
validation <- zcta_annual_direct %>%
  left_join(zcta_annual_model_assn, by=c("ZCTA"="zcta")) %>% # zcta_annual_model_pop
  left_join(zcta_annual_model_pop, by=c("ZCTA"="zcta")) %>%
  
  # filter(!is.na(total_CVD_model)) %>%
  filter(ZCTA %in% unique(focus_zctas$ZCTA5)) %>%
  left_join(zcta_pop_2010, by=c("ZCTA"="zcta")) %>%
  
  # Rates of CVD Hospitalization Per Capita 
  mutate(total_CVD_direct_rate = total_CVD_direct/pop_count_zcta10) %>%

  # Mean Squared Error, Mean Absolute difference 
  mutate(mse_sumterm = (total_CVD_model-total_CVD_direct)^2,
         mad_sumterm  = abs(total_CVD_model-total_CVD_direct),
         mape_sumterm = mad_sumterm/total_CVD_direct) %>%
  mutate(mse_sumterm_wgt = (total_CVD_model_wgt-total_CVD_direct)^2,
         mad_sumterm_wgt  = abs(total_CVD_model_wgt-total_CVD_direct),
         mape_sumterm_wgt = mad_sumterm_wgt/total_CVD_direct)
  
# 3. Compare distributions of CVD hospitalization rates per capita @ CBG vs. ZCTA ----

cbg_cvd_rates <- cbg_mrp_results %>%
  mutate(total_CVD_model_rate = total_CVD_model_cbg/pop_count) # pop_count is census population 65+

# Descriptive Statistics
summary(validation$total_CVD_direct_rate)
summary(cbg_cvd_rates$total_CVD_model_rate)

rate_cbg_ds <- cbg_cvd_rates %>%
  mutate(group = "CBG Modeled Annual CVD Hospitalizations Per Capita") %>%
  group_by(group) %>%
  summarise(min = min(total_CVD_model_rate),
            p25 = quantile(total_CVD_model_rate, .25),
            mean = mean(total_CVD_model_rate),
            median = median(total_CVD_model_rate),
            p75 = quantile(total_CVD_model_rate, .75),
            max = max(total_CVD_model_rate))

rate_zcta_ds <- validation %>%
  mutate(group = "ZCTA Actual Annual CVD Hospitalizations Per Capita") %>%
  group_by(group) %>%
  summarise(min = min(total_CVD_direct_rate),
            p25 = quantile(total_CVD_direct_rate, .25),
            mean = mean(total_CVD_direct_rate),
            median = median(total_CVD_direct_rate),
            p75 = quantile(total_CVD_direct_rate, .75),
            max = max(total_CVD_direct_rate))

rates <- rate_cbg_ds %>% bind_rows(rate_zcta_ds)

write.csv(rates, "./Data/Results/Validation/Compare_Rates_CBG_vs_ZCTA_DescStats.csv")

# 4. Calculate Validation Metrics at ZCTA -----

# Correlation Coefficients
cor.test(validation$total_CVD_direct, validation$total_CVD_model, method="pearson") 
cor.test(validation$total_CVD_direct, validation$total_CVD_model, method="spearman")

cor.test(validation$total_CVD_direct, validation$total_CVD_model_wgt, method="pearson")
cor.test(validation$total_CVD_direct, validation$total_CVD_model_wgt, method="spearman")

# Mean Squared Error - how is Zhang expressing MSE and MAD as a percentage?
(1/nrow(validation[!is.na(validation$total_CVD_model),]))*sum(validation$mse_sumterm, na.rm=T)

# Mean Absolute Difference
(1/nrow(validation[!is.na(validation$total_CVD_model),]))*sum(validation$mad_sumterm, na.rm=T)

# Mean Absolute Percentage Error
(1/nrow(validation[!is.na(validation$total_CVD_model),]))*sum(validation$mape_sumterm, na.rm=T)


# Export: 
zcta_aggregated_comparison_stats <- data.frame(stat = c("Pearson Correlation Coefficient", "Spearman Correlation Coefficient",
                                                        "Pearson Correlation Coefficient", "Spearman Correlation Coefficient",
                                                        "Mean Squared Error", "Mean Absolute Difference", "Mean Absolute Percentage Error",
                                                        "Mean Squared Error", "Mean Absolute Difference", "Mean Absolute Percentage Error"),
                                               aggregation_type = c("Assigned", "Assigned", "Weighted", "Weighted", "Assigned", "Assigned", "Assigned",
                                                                    "Weighted", "Weighted", "Weighted"),
                                               value = c(cor.test(validation$total_CVD_direct, validation$total_CVD_model, method="pearson")$estimate,
                                                         cor.test(validation$total_CVD_direct, validation$total_CVD_model, method="spearman")$estimate,
                                                         cor.test(validation$total_CVD_direct, validation$total_CVD_model_wgt, method="pearson")$estimate,
                                                         cor.test(validation$total_CVD_direct, validation$total_CVD_model_wgt, method="spearman")$estimate,
                                                         (1/nrow(validation[!is.na(validation$total_CVD_model),]))*sum(validation$mse_sumterm, na.rm=T),
                                                         (1/nrow(validation[!is.na(validation$total_CVD_model),]))*sum(validation$mad_sumterm, na.rm=T),
                                                         (1/nrow(validation[!is.na(validation$total_CVD_model),]))*sum(validation$mape_sumterm, na.rm=T),
                                                         (1/nrow(validation[!is.na(validation$total_CVD_model_wgt),]))*sum(validation$mse_sumterm_wgt, na.rm=T),
                                                         (1/nrow(validation[!is.na(validation$total_CVD_model_wgt),]))*sum(validation$mad_sumterm_wgt, na.rm=T),
                                                         (1/nrow(validation[!is.na(validation$total_CVD_model_wgt),]))*sum(validation$mape_sumterm_wgt, na.rm=T)))

write.csv(zcta_aggregated_comparison_stats, "./Data/Results/Validation/Compare_ZCTA_Agg_vs_Actual_Statistics.csv")

# Descriptive Statistics - Overall Distributions 
summary(validation$total_CVD_direct) # 84 ZCTAs
summary(validation$total_CVD_model) # 79 ZCTAs
summary(validation$total_CVD_model_wgt) # 39 ZCTAs

# Descriptive Statistics - subset of 39 ZCTAs w/ total_CVD_model_wgt
subset <- validation %>% filter(!is.na(total_CVD_model_wgt))
summary(subset$total_CVD_direct)
summary(subset$total_CVD_model)
summary(subset$total_CVD_model_wgt)
  # comparing on same set - model may be slightly overestimating but this is before limiting to summer which makes sense

zcta_ds_all <- validation %>%
  select(ZCTA, total_CVD_direct, total_CVD_model, total_CVD_model_wgt) %>%
  pivot_longer(cols=total_CVD_direct:total_CVD_model_wgt, names_to = "group", values_to="total_CVD") %>%
  group_by(group) %>%
  summarise(min = min(total_CVD, na.rm=T),
            p25 = quantile(total_CVD, .25, na.rm=T),
            mean = mean(total_CVD, na.rm=T),
            median = median(total_CVD, na.rm=T),
            p75 = quantile(total_CVD, .75, na.rm=T),
            max = max(total_CVD, na.rm=T),
            n = sum(!is.na(total_CVD))) %>%
  mutate(comparison = "All Non-Missing")

zcta_ds_subset <- subset %>%
  select(ZCTA, total_CVD_direct, total_CVD_model, total_CVD_model_wgt) %>%
  pivot_longer(cols=total_CVD_direct:total_CVD_model_wgt, names_to = "group", values_to="total_CVD") %>%
  group_by(group) %>%
  summarise(min = min(total_CVD, na.rm=T),
            p25 = quantile(total_CVD, .25, na.rm=T),
            mean = mean(total_CVD, na.rm=T),
            median = median(total_CVD, na.rm=T),
            p75 = quantile(total_CVD, .75, na.rm=T),
            max = max(total_CVD, na.rm=T),
            n = sum(!is.na(total_CVD))) %>%
  mutate(comparison = "Same Subset")

zcta_ds <- zcta_ds_all %>% bind_rows(zcta_ds_subset)

write.csv(zcta_ds, "./Data/Results/Validation/Compare_ZCTA_Agg_vs_Actual_DescStats.csv")
