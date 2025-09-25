################################################################################
# Program Name: Tables_and_Figures.R
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
  library(colorspace)
  library(tmap)
  library(waterfalls)
  library(cowplot)
}

### Prepare Color Palettes ---

pal <- wes_palette("Zissou1", n = 5) #  Darjeeling1
pal

pal2 <- wes_palette("Darjeeling1", n = 5) #  
pal2

pal_blue <- c(lighten(pal[5],0.6), lighten(pal[5],0.3),pal[5],darken(pal[5],0.3),darken(pal[5],0.6))
pal_red <- c(lighten(pal[1],0.6), lighten(pal[1],0.3),pal[1],darken(pal[1],0.3),darken(pal[1],0.6))
pal_green <- c(lighten(pal[2],0.6), lighten(pal[2],0.3),pal[2],darken(pal[2],0.3),darken(pal[2],0.6))

pal_red_cont <- colorRampPalette(c('white', pal2[1]))
pal_green_cont <- colorRampPalette(c("#84EFD8", "#034036"))

### Figure Data ---
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### 1. BRING IN DATA ----
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Step 3a Results
cbg_an_final_25 <- read_csv("../Data/Analysis/RTP_XCESS_ATTR_Burden_CBG_Summer2018_25C_povallhh.csv") %>% # "../Data/Analysis/3_RTP_XCESS_ATTR_Burden_CBG_Summer2018_25C.csv"
  filter(!is.na(an_c_1km)) %>%
  mutate(month = month(date)) %>%
  filter(month %in% c(5,6,7,8,9)) # drop out the lag days in April and October

cbg_comp <- cbg_an_final_25 %>% 
  group_by(GEOID, pop_count) %>%
  summarise(tot_an = sum(an_c_1km),
            tot_an_tps = sum(an_c_ws_tps),
            tmean_f = mean(tmean_f_1km),
            tmean_ftps = mean(tmean_f_ws_tps),
            tmean_c = mean(tmean_c_1km),
            tmean_ctps = mean(tmean_c_ws_tps),         
            total_cvd = sum(total_CVD)) %>%
  ungroup() %>%
  mutate(GEOID = as.character(GEOID))

# TIGER/Line Shapefile, 2020, North Carolina, Block Groups: https://catalog.data.gov/dataset/tiger-line-shapefile-2020-state-north-carolina-block-groups
cbg_sf <- st_read("../Data/Original/Geography/tl_2020_37_bg/tl_2020_37_bg.shp") %>%
  select(GEOID, geometry)

cbg_comp_sf <- cbg_sf %>%
  select(GEOID, geometry) %>%
  left_join(cbg_comp, by="GEOID") %>%
  filter(!is.na(tot_an)) # drop out CBGs in NC outside of the RT focus area

# Step 3b Results
decomp <- read_csv("../Data/Analysis/RTP_XCESS_ATTR_Burden_CBG_Decomposed_povallhh.csv") %>% # "../Data/Results/3b_Attributable_Rate_Decomposition.csv"
  select(-c(pop_count)) %>%
  mutate(GEOID = as.character(GEOID))

cbg_comp_sf_plot <- cbg_comp_sf %>%
  mutate(tot_an_rate = (tot_an/pop_count)*10000,
         tot_an_rate_tps = (tot_an_tps/pop_count)*10000,
         tot_cvd_rate = (total_cvd/pop_count)*10000) %>%
  mutate(temp_diff = tmean_f - tmean_ftps,
         an_diff = tot_an - tot_an_tps,
         an_rate_diff = tot_an_rate - tot_an_rate_tps) %>%
  mutate(pct_ha = (tot_an/total_cvd)*100) %>%
  left_join(decomp, by="GEOID") %>%
  mutate(cvd_pctile_alt = round(percent_rank(tot_cvd_rate)*100, digits=0),
         temp_pctile_alt = round(percent_rank(tmean_c)*100, digits=0)) %>%
  # For subgroup burden plots
  mutate(tot_an_rate = (tot_an/pop_count)*10000) %>%
  mutate(an_qtile = round(percent_rank(tot_an)*100, digits=0),
         an_rate_qtile = round(percent_rank(tot_an_rate)*100, digits=0)) %>%
  mutate(decile = ntile(tot_an_rate,10)) %>%
  mutate(top_ten = ifelse(decile == 10,"Top 10%","Bottom 90%"))

rm(cbg_comp, cbg_sf, cbg_comp_sf, decomp, cbg_an_final_25)

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Table 3:  Multi-Level Regression Model and Attributable Burden Analysis Results
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# Note: This table combines the modeling results from Step 2 and subgroup burdens from Step 3c

step2_results <- read.csv("../Data/Figures/MRP_Model_Results.csv") %>% # "../Data/Results/2_MRP_Model_Results.csv"
  mutate(subgroup_type = case_when(grepl("age",X) ~ "age",
                                   grepl("sex",X) ~ "sex",
                                   grepl("race",X) ~ "race"),
         subgroup_value = sub(paste0(".*","_factor"),"",X),
         confint_95 = paste0("(",round(`X2.5..`,digits=2),",",round(`X97.5..`,digits=2),")")) %>%
  rename(odds_ratio = "exp.fixef.model_agg..") %>%
  select(subgroup_type, subgroup_value, odds_ratio, confint_95)

step3c_results <- read.csv("../Data/Results/3c_Subgroup_Attributable_Burdens.csv") %>%
  select(subgroup_type, subgroup_value, pop_count, tot_an, rate_an) %>%
  rename(population = pop_count,
         tot_cvd_hosp_count = tot_an,
         cvd_hosp_rate_per10k = rate_an)

table3 <- step2_results %>%
  full_join(step3c_results, by=c("subgroup_type","subgroup_value")) %>%
  arrange(subgroup_type,subgroup_value)

rm(step2_results, step3c_results)
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Figure 2: Characteristics of Census Block Groups by Heat-Attributable CVD Hospitalization Rates  ----
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

sex_totals <- census_data %>%
  group_by(GEOID, sex_factor) %>%
  summarise(pop_count = sum(pop_count)) %>% 
  pivot_wider(id_cols = GEOID, names_from = sex_factor, values_from = pop_count) %>%
  rename(pop_male = Male,
         pop_female = Female)

race_totals <- census_data %>%
  group_by(GEOID, race_factor) %>%
  summarise(pop_count = sum(pop_count)) %>%
  pivot_wider(id_cols = GEOID, names_from = race_factor, values_from = pop_count) %>%
  rename(pop_asian = Asian,
         pop_black = Black,
         pop_hispanic = Hispanic,
         pop_native = NANative,
         pop_other = Other,
         pop_white = White)

age_totals <- census_data %>%
  group_by(GEOID, age_factor) %>%
  summarise(pop_count = sum(pop_count)) %>%
  pivot_wider(id_cols = GEOID, names_from = age_factor, values_from = pop_count) %>%
  rename(pop_65_74 = `(64,74]`,
         pop_75_84 = `(74,84]`,
         pop_85_pl = `(84,110]`)


# Poverty + Education Data
acs_cbg_orig <- read_csv("../Data/Original/Census/nhgis0025_ds249_20205_blck_grp_E.csv") # ACS 2016-2020
acs_cbg <- acs_cbg_orig %>%
  dplyr::select(GEOID, AMRZE001, AMRZE022, AMRZE023, AMRZE024, AMRZE025, AMRZE025, AMR5E001, AMR5E002, AMR5E008, AMR5E014, AMR5E019, AMR5E025, AMR5E030, AMR5E037, AMR5E043, AMR5E048, AMR5E054, AMR5E059, AMR8E001) %>%
  separate(GEOID, into=c("Geography", "GEOID"), sep="US") %>%
  # Educational Attainment for the Population 25 Years and Over, ACS 2016-2020
  mutate(pop25_w_bach = AMRZE022+AMRZE023+AMRZE024+AMRZE025,
         tot_pop25 = AMRZE001,
         pct_w_bach = (AMRZE022+AMRZE023+AMRZE024+AMRZE025)/AMRZE001) %>%
  
  # Median Household Income in the Past 12 Months (in 2020 Inflation-Adjusted Dollars), ACS 2016-2020
  rename(median_hh_inc = AMR8E001,
         hh_below_pl = AMR5E002,
         tot_hh = AMR5E001) %>%
  
  # Households with Income in Past 12 Months Below Poverty Level, ACS 2016-2020
  mutate(pct_hh_below_pl = hh_below_pl/tot_hh) %>%
  # Format
  mutate(GEOID = as.numeric(GEOID)) %>%
  dplyr::select(GEOID, pop25_w_bach, tot_pop25, pct_w_bach, median_hh_inc, hh_below_pl, tot_hh, pct_hh_below_pl) %>%
  # Filter to CBGs in RDU
  filter(GEOID %in% unique(cbg_comp_sf_plot$GEOID))

cbg_pop_totals <- sex_totals %>%
  left_join(race_totals, by="GEOID") %>%
  left_join(age_totals, by="GEOID") %>%
  left_join(acs_cbg, by="GEOID") %>%
  mutate(GEOID = as.character(GEOID))

hotspots <- cbg_comp_sf_plot %>%
  left_join(cbg_pop_totals, by="GEOID") %>%
  st_drop_geometry() %>%
  group_by(top_ten) %>%
  summarise(pop_count = sum(pop_count),
            pop_male = sum(pop_male),
            pop_female = sum(pop_female),
            pop_asian = sum(pop_asian),
            pop_black = sum(pop_black),
            pop_hispanic = sum(pop_hispanic),
            pop_native = sum(pop_native),
            pop_other = sum(pop_other),
            pop_white = sum(pop_white),
            pop_65_74 = sum(pop_65_74),
            pop_75_84 = sum(pop_75_84),
            pop_85_pl = sum(pop_85_pl),
            pop25_w_bach = sum(pop25_w_bach),
            tot_pop25 = sum(tot_pop25),
            hh_below_pl = sum(hh_below_pl),
            tot_hh = sum(tot_hh)) %>%
  ungroup() %>%
  mutate(pct_male = pop_male/(pop_male+pop_female),
         pct_female = pop_female/(pop_male+pop_female),
         pct_asian = pop_asian/(pop_asian + pop_black + pop_hispanic + pop_native + pop_other + pop_white),
         pct_black = pop_black/(pop_asian + pop_black + pop_hispanic + pop_native + pop_other + pop_white),
         pct_hispanic = pop_hispanic/(pop_asian + pop_black + pop_hispanic + pop_native + pop_other + pop_white),
         pct_native = pop_native/(pop_asian + pop_black + pop_hispanic + pop_native + pop_other + pop_white),
         pct_other = pop_other/(pop_asian + pop_black + pop_hispanic + pop_native + pop_other + pop_white),
         pct_white = pop_white/(pop_asian + pop_black + pop_hispanic + pop_native + pop_other + pop_white),
         pct_65_74 = pop_65_74/(pop_65_74 + pop_75_84 + pop_85_pl),
         pct_75_84 = pop_75_84/(pop_65_74 + pop_75_84 + pop_85_pl),
         pct_85_pl = pop_85_pl/(pop_65_74 + pop_75_84 + pop_85_pl),
         pct_below_pl = hh_below_pl/tot_hh,
         pct_w_bach = pop25_w_bach/tot_pop25) %>%
  dplyr::select(top_ten, pct_male, pct_female, pct_asian, pct_black,
                pct_hispanic, pct_native, pct_other, pct_white,
                pct_65_74, pct_75_84, pct_85_pl, pct_below_pl, pct_w_bach) %>%
  pivot_longer(pct_male:pct_w_bach, names_to="var", values_to = "pct") %>%
  mutate(subgroup_type = case_when(var %in% c("pct_male", "pct_female") ~ "Sex",
                                   var %in% c("pct_asian", "pct_black", "pct_hispanic", "pct_native", "pct_other", "pct_white") ~ "Race",
                                   var %in% c("pct_65_74", "pct_75_84", "pct_85_pl") ~ "Age",
                                   var %in% c("pct_below_pl", "pct_w_bach") ~ "Poverty + Education")) %>%
  # For Labels:
  mutate(var_label = case_when(var == "pct_male" ~ "Male",
                               var == "pct_female" ~ "Female",
                               var == "pct_asian" ~ "Asian",
                               var == "pct_black" ~ "Black",
                               var == "pct_hispanic" ~ "Hisp.",
                               var == "pct_native" ~ "Native American",
                               var == "pct_other" ~ "Other",
                               var == "pct_white" ~ "White",
                               var == "pct_65_74" ~ "65-74",
                               var == "pct_75_84" ~ "75-64",
                               var == "pct_85_pl" ~ "85+",
                               var == "pct_below_pl" ~ "HHs Below PL",
                               var == "pct_w_bach" ~ "Bachelor's +")) %>%
  mutate(pct_round = round(pct*100,1))

# Sex  
sex_comp <- ggplot(hotspots[hotspots$subgroup_type=="Sex",],
                   aes(x=var_label, y=pct_round, fill=top_ten)) + # , color=top_ten
  geom_bar(position="dodge",stat="identity") + # ,color="black"
  geom_text(aes(label=pct_round, size=7), vjust=-0.5, position=position_dodge(width=0.9), size=3.75) +
  labs(y= "Percent", x = "Sex", fill = "CBG Attr. Rate Group") +
  scale_y_continuous(limits=c(0,100), n.breaks=10) +
  scale_fill_manual(values = c("Top 10%"="#034036", "Bottom 90%"="#84EFD8")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15))

# Race
race_comp <- ggplot(hotspots[hotspots$subgroup_type=="Race" & hotspots$var!="pct_native",],
                    aes(x=var_label, y=pct_round, fill=top_ten)) + # exclude native bc not in model, , color=top_ten
  geom_bar(position="dodge",stat="identity") +
  geom_text(aes(label=pct_round, size=7), vjust=-0.5, position=position_dodge(width=0.9), size=3.75) +
  labs(y= "Percent", x = "Race", fill = "CBG Attr. Rate Group") +
  scale_y_continuous(limits=c(0,100), n.breaks=10) +
  scale_fill_manual(values = c("Top 10%"="#034036", "Bottom 90%"="#84EFD8")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15))

# Age
age_comp <- ggplot(hotspots[hotspots$subgroup_type=="Age",], 
                   aes(x=var_label, y=pct_round, fill=top_ten)) + # , color=top_ten
  geom_bar(position="dodge",stat="identity") +
  geom_text(aes(label=pct_round, size=7), vjust=-0.5, position=position_dodge(width=0.9), size=3.75) +
  labs(y= "Percent", x = "Age", fill = "CBG Attr. Rate Group") +
  scale_y_continuous(limits=c(0,100), n.breaks=10) +
  scale_fill_manual(values = c("Top 10%"="#034036", "Bottom 90%"="#84EFD8")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15))

# Poverty + Education
other_comp <- ggplot(hotspots[hotspots$subgroup_type=="Poverty + Education",], 
                     aes(x=var_label, y=pct_round, fill=top_ten)) + # , color=top_ten
  geom_bar(position="dodge",stat="identity") +
  geom_text(aes(label=pct_round, size=7), vjust=-0.5, position=position_dodge(width=0.9), size=3.75) +
  labs(y= "Percent", x = "Poverty + Education", fill = "CBG Attr. Rate Group") +
  scale_y_continuous(limits=c(0,100), n.breaks=10) +
  scale_fill_manual(values = c("Top 10%"="#034036", "Bottom 90%"="#84EFD8")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15))

bar_comp <- ggarrange(sex_comp, race_comp, age_comp, other_comp, 
                      nrow=2, ncol=2, common.legend = T, legend="bottom") 

# MAP
tot_an_rate_sf <- ggplot(cbg_sf) +
  geom_sf(aes(fill=an_rate_qtile, color=top_ten),
          lwd = ifelse(cbg_sf$top_ten =="Top 10%", 1, 0.1)) +
  scale_fill_gradient(name = "CBG Attr. Rate Quantile",
                      high = "#034036", low = "#84EFD8", ) + 
  scale_color_manual(values=c("Bottom 90%"="black","Top 10%"=pal2[1]),
                     name = "CBG Attr. Rate Group") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.box.just = "left",
        axis.text = element_blank(),
        panel.grid=element_blank(),
        panel.background=element_rect(fill="white"),
        text=element_text(size=13),
        legend.text=element_text(size=12)) + #11.5 for 6/12) +
  ggtitle("Heat-Attributable CVD Hospitalization Rate per 10k")
tot_an_rate_sf

# COMBINE
ggarrange(bar_comp, tot_an_rate_sf, ncol=2, nrow=1, widths=c(1,1.15)) # 
ggsave("../Data/Figures/Fig2_CBG_Characteristics_by_AttrRate.png", height=7, width=14, units="in")

rm(hotspots, cbg_pop_totals, sex_totals, race_totals, age_totals,
   acs_cbg_orig, acs_cbg, sex_comp, race_comp, age_comp, other_comp, bar_comp)

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Figure 3: Temperature Data Comparison (Temperature and Relative Burden) ----
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Note: this code creates individual map components (final figure created in Canva)

temp_1km <- tm_shape(cbg_comp_sf_plot) +
  tm_polygons("tmean_f", 
              style="cont", 
              title="Degrees F",
              palette = pal_red_cont(100)) + # "Reds"
  tm_layout(main.title="High Resolution (1km Gridded)",
            legend.position=c("left","TOP"),
            main.title.size=1)

tmap_save(temp_1km, "../Data/Figures/temp_1km.png", height=4.5, width=4)


temp_tps <- tm_shape(cbg_comp_sf_plot) +
  tm_polygons("tmean_ftps", 
              style="cont", 
              title="Degrees F",
              palette = pal_red_cont(100)) + # "Reds"
  tm_layout(main.title="Low Resolution (WS Interpolated)",
            legend.position=c("left","TOP"),
            main.title.size=1)

tmap_save(temp_tps, "../Data/Figures/temp_tps.png", height=4.5, width=4)

tmap_save(temp_diff, "../Data/Figures/temp_diff.png", height=4.5, width=4)

anr_1km <- tm_shape(cbg_comp_sf_plot) +
  tm_polygons("tot_an_rate",
              title = "AN per 10k",
              style="quantile",
              palette = pal_green_cont(10)) +
  tm_layout(main.title="High Resolution (1km Gridded)",
            legend.position=c("left","TOP"),
            main.title.size=1) 

tmap_save(anr_1km, "../Data/Figures/anr_1km.png", height=4.5, width=4)

anr_tps <- tm_shape(cbg_comp_sf_plot) +
  tm_polygons("tot_an_rate_tps",
              title = "AN per 10k",
              style="quantile",
              palette = pal_green_cont(10)) +
  tm_layout(main.title="Low Resolution (WS Interpolated)",
            legend.position=c("left","TOP"),
            main.title.size=1) 

tmap_save(anr_tps, "../Data/Figures/anr_tps.png", height=4.5, width=4)

tmap_save(rate_diff, "../Data/Figures/rate_diff.png", height=4.5, width=4)

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Figure 4: Characteristics of CBGs by Risk Group ----
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Note: this code creates individual map components (final figure with quadrant legend created in Canva)

# Clusters 
cluster <- cbg_comp_sf_plot %>%
  mutate(cluster = case_when(an_tempdiff2_rate>0 & an_cvddiff1_rate>0 ~ "Above Avg Temp, Above Avg CVD",
                             an_tempdiff2_rate>0 & an_cvddiff1_rate<0 ~ "Above Avg Temp, Below Avg CVD",
                             an_tempdiff2_rate<0 & an_cvddiff1_rate>0 ~ "Below Avg Temp, Above Avg CVD",
                             an_tempdiff2_rate<0 & an_cvddiff1_rate<0 ~ "Below Avg Temp, Below Avg CVD")) %>%
  mutate(risk_group = case_when(an_tempdiff2_rate>0 & an_cvddiff1_rate>0 ~ "Dual Channel Risk",
                                an_tempdiff2_rate>0 & an_cvddiff1_rate<0 ~ "Heat Driven Risk",
                                an_tempdiff2_rate<0 & an_cvddiff1_rate>0 ~ "Health Driven Risk",
                                an_tempdiff2_rate<0 & an_cvddiff1_rate<0 ~ "Low Risk"))

# Cluster Demographic Data
acs_cbg_orig <- read_csv("../Data/Original/Census/nhgis0025_ds249_20205_blck_grp_E.csv") # ACS 2016-2020
acs_cbg_housing <- read_csv("../Data/Original/Census/CBG_Housing_ACS/nhgis0029_csv/nhgis0029_ds249_20205_blck_grp.csv") %>%
  select(GEOID, AMUFE001, AMUFE002, AMUFE003)

acs_cbg_ethnicity <- read_csv("../Data/Original/Census/CBG_Ethnicity_ACS/nhgis0030_csv/nhgis0030_ds249_20205_blck_grp.csv") %>%
  select(GEOID, AMP3E001, AMP3E003)

dhc_cbg_urbrur <- read_csv("../Data/Original/Census/CBG_UrbanRural_DHC/nhgis0031_csv\\nhgis0031_ds258_2020_blck_grp.csv") %>%
  separate(GEOID, into=c("Geography", "GEOID"), sep="US") %>% 
  select(-Geography) %>%
  select(GEOID, U7I001, U7I002, U7I003, U9W001,U9W002,U9W003) # using population not housing units

acs_cbg <- acs_cbg_orig %>%
  select(STATE, GEOID, AMRZE001, AMRZE022, AMRZE023, AMRZE024, AMRZE025, AMRZE025, 
         AMR5E001, AMR5E002, AMR5E008, AMR5E014, AMR5E019, AMR5E025, AMR5E030, 
         AMR5E037, AMR5E043, AMR5E048, AMR5E054, AMR5E059, AMR8E001) %>%
  filter(STATE == "North Carolina") %>%
  left_join(acs_cbg_housing, by="GEOID") %>%
  left_join(acs_cbg_ethnicity, by="GEOID") %>%
  separate(GEOID, into=c("Geography", "GEOID"), sep="US") %>%
  select(-Geography) %>%
  left_join(dhc_cbg_urbrur, by="GEOID")

rm(acs_cbg_orig, acs_cbg_ethnicity, acs_cbg_housing, dhc_cbg_urbrur)

cluster_demographics <- cluster %>%
  left_join(acs_cbg, by="GEOID") %>%
  st_drop_geometry() %>%
  group_by(risk_group) %>% # cluster
  summarise(degree_bach = sum(AMRZE022, na.rm=T),
            degree_mast = sum(AMRZE023, na.rm=T),
            degree_prof = sum(AMRZE024, na.rm=T),
            degree_doct = sum(AMRZE024, na.rm=T),
            tot_pop_25_pl = sum(AMRZE001, na.rm=T),
            hh_below_pl = sum(AMR5E002, na.rm=T),
            tot_hh = sum(AMR5E001, na.rm=T),
            tot_hu = sum(AMUFE001, na.rm=T),
            oo_hu = sum(AMUFE002, na.rm=T),
            ro_hu = sum(AMUFE003, na.rm=T),
            tot_pop = sum(AMP3E001, na.rm=T),
            white_alone = sum(AMP3E003, na.rm=T),
            tot_pop2 = sum(U7I001, na.rm=T),
            urban_pop = sum(U7I002, na.rm=T),
            rural_pop = sum(U7I003, na.rm=T)) %>%
  ungroup() %>%
  mutate(pct_w_bach = ((degree_bach+degree_mast+degree_prof+degree_doct)/tot_pop_25_pl)*100,
         pct_hh_below_pl = (hh_below_pl/tot_hh)*100,
         pct_hu_rented = (ro_hu/tot_hu)*100,
         pct_nonwhite = ((tot_pop-white_alone)/tot_pop)*100,
         pct_rural_pop = (rural_pop/tot_pop2)*100) %>%
  select(risk_group, pct_w_bach:pct_rural_pop) %>%
  pivot_longer(pct_w_bach:pct_rural_pop, names_to = "metric", values_to="percentage") %>%
  mutate(labels = case_when(metric == "pct_w_bach" ~ "Pop w/ Bachelor's Degree +",
                            metric == "pct_hu_rented" ~ "Renter Occupied Housing Units",
                            metric == "pct_nonwhite" ~ "Non-White Population",
                            metric == "pct_rural_pop" ~ "Rural Population",
                            metric == "pct_hh_below_pl" ~ "Households Below Poverty Level"))

# PLOT TOGETHER:
cluster_map_gg <- ggplot(cluster) +
  geom_sf(aes(fill=risk_group)) +
  scale_fill_manual(name = "Risk Group",
                    breaks = c("Low Risk", "Health Driven Risk", "Heat Driven Risk", "Dual Channel Risk"),
                    values = c("#8AE3FE", "#00A08A","#fdae61","#F21A00")) +
  # values=c(lighten(pal2[5],0.5), pal2[2], "#fdae61", pal[5]))) 
  theme_bw() +
  theme(axis.text = element_blank(),
        legend.position = "none",
        panel.grid=element_blank(),
        legend.text=element_text(size=16),
        plot.background = element_blank())

cluster_map_gg
ggsave("../Data/Figures/Fig4_Maps.png", height=4, width=4, units="in", dpi=300)

demo_bar <- ggplot(cluster_demographics) +
  geom_bar(aes(x=labels, y=percentage, group=risk_group, fill=risk_group), 
           position=position_dodge(), stat="identity") + # , colour="black"
  scale_fill_manual(name = "Risk Group",
                    breaks = c("Low Risk", "Health Driven Risk", "Heat Driven Risk", "Dual Channel Risk"),
                    values = c("#8AE3FE", "#00A08A","#fdae61","#F21A00")) +
  # coord_cartesian(ylim=c(74.5,76.5)) +
  labs(y= "Percent", x = "") +
  theme_bw() +
  theme(legend.position = "none",
        plot.background = element_rect(fill="white"),
        legend.text=element_text(size=16), # 16 for 6/12
        text = element_text(size = 11.5), #11.5 for 6/12
        panel.grid=element_blank()) + 
  # theme(legend.position = "bottom",legend.text = element_text(size=8)) + # axis.text.x = element_text(angle = 45, hjust = 1), +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))

demo_bar
ggsave("../Data/Figures/Fig4_Bars.png", height=4, width=4, units="in", dpi=300)


ggarrange(cluster_map_gg, demo_bar, align = "h")
combined_cluster <- ggarrange(cluster_map_gg, demo_bar, align = "h") # common.legend=T, legend="bottom", 
combined_cluster

# Estimate Population by Cluster - values in text
total_pop_by_cbg <- read_csv("../Data/Original/Census/Census_P12_Tables.csv") %>%
  select(GEOID, race_factor, `!!Total:`) %>%
  rename(total_pop = `!!Total:`) %>%
  group_by(GEOID) %>%
  summarise(total_pop = sum(total_pop))

cluster_pop <- cluster %>%
  group_by(risk_group) %>%
  summarise(count_cbg = n(),
            pop_count = sum(pop_count)) %>%
  ungroup() %>%
  mutate(total_pop = sum(pop_count),
         total_cbg = sum(count_cbg)) %>%
  mutate(pct_pop = pop_count/total_pop,
         pct_cbg = count_cbg/total_cbg)

