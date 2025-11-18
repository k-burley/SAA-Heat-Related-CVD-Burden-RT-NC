################################################################################
# Program Name: RShiny_Dashboard.R
# Program Purpose: Create R Shiny Dashboard with results from RTP CVD small area analysis

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
  library(shiny)
  library(shinydashboard)
  library(leaflet)
  library(fresh)
}


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### 1. BRING IN DATA ----
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Set up working directories, bring in necessary files. 
setwd("~/GitHub/SAA-Heat-Related-CVD-Burden-RT-NC/Code")

cbg_an_final_25 <- read_csv("../Data/Analysis/3_RTP_XCESS_ATTR_Burden_CBG_Summer2018_25C.csv") %>% 
  filter(!is.na(an_c_1km)) %>%
  mutate(LocationID = str_sub(GEOID, 1, 11)) %>% # get census tract from CBG - 459
  mutate(month = month(date)) %>%
  filter(month %in% c(5,6,7,8,9))

cbg_comp <- cbg_an_final_25 %>% 
  group_by(GEOID, LocationID, pop_count) %>%
  summarise(tot_an = sum(an_c_1km),
            tmean_f = mean(tmean_f_1km),
            tmean_c = mean(tmean_c_1km),
            total_cvd = sum(total_CVD)) %>%
  ungroup() %>%
  mutate(GEOID = as.character(GEOID))

cbg_sf <- st_read("../Data/Original/Geography/tl_2020_37_bg/tl_2020_37_bg.shp") %>%
  select(GEOID, geometry)

cbg_comp_sf <- cbg_sf %>%
  select(GEOID, geometry) %>%
  left_join(cbg_comp, by="GEOID") %>%
  filter(!is.na(tot_an)) 

# Attributable Burden decomposition
decomp <- read_csv("../Data/Results/3b_Attributable_Rate_Decomposition.csv") %>%
  select(-c(pop_count)) %>%
  mutate(GEOID = as.character(GEOID)) 

cbg_comp_sf_plot <- cbg_comp_sf %>%
  mutate(tot_an_rate = (tot_an/pop_count)*10000,
         tot_cvd_rate = (total_cvd/pop_count)*10000) %>%
  mutate(pct_ha = (tot_an/total_cvd)*100) %>%
  left_join(decomp, by="GEOID") %>%
  mutate(cvd_pctile = round(percent_rank(tot_cvd_rate)*100, digits=0),
         temp_pctile = round(percent_rank(tmean_c)*100, digits=0)) %>%
  st_transform(4326) # check if this fixes!

rm(cbg_an_final_25, cbg_comp, cbg_comp_sf, decomp)

# MAIN DF FOR THE MAP! Risk Groups and Rate
map_df <- cbg_comp_sf_plot %>%
  mutate(cluster = as.factor(case_when(an_tempdiff2_rate>0 & an_cvddiff1_rate>0 ~ "Dual Channel Risk",
                             an_tempdiff2_rate>0 & an_cvddiff1_rate<0 ~ "Heat Driven Risk",
                             an_tempdiff2_rate<0 & an_cvddiff1_rate>0 ~ "Health Driven Risk",
                             an_tempdiff2_rate<0 & an_cvddiff1_rate<0 ~ "Low Risk"))) %>%
  mutate(anrate_pctile = round(percent_rank(tot_an_rate)*100, digits=0)) %>%
  select(GEOID, geometry, cluster, an_cvddiff1_rate, an_tempdiff2_rate, 
         cvd_pctile, temp_pctile, tmean_c, tmean_f, tot_cvd_rate,
         tot_an_rate, anrate_pctile) %>%
  # Popup info
  mutate(content = paste0("CBG ID: ",GEOID)) %>%
  mutate(cvd_text = case_when((cvd_pctile %% 10 == 1 & cvd_pctile!=11) ~ paste0(round(tot_cvd_rate, 1)," hosp. per 10k"),
                                  (cvd_pctile %% 10 == 2 & cvd_pctile!=12) ~ paste0(round(tot_cvd_rate, 1)," hosp. per 10k"),
                                  (cvd_pctile %% 10 == 3 & cvd_pctile!=13) ~ paste0(round(tot_cvd_rate, 1)," hosp. per 10k"),
                                  TRUE ~ paste0(round(tot_cvd_rate, 1)," hosp. per 10k"))) %>%
  mutate(temp_text = case_when((temp_pctile %% 10 == 1 & temp_pctile!=11) ~ paste0(round(tmean_f, 1),"°F"),
                                   (temp_pctile %% 10 == 2 & temp_pctile!=12) ~ paste0(round(tmean_f, 1),"°F"),
                                   (temp_pctile %% 10 == 3 & temp_pctile!=13) ~ paste0(round(tmean_f, 1),"°F"),
                                   TRUE ~ paste0(round(tmean_f, 1),"°F"))) %>%
  mutate(anr_pct_text = case_when((anrate_pctile %% 10 == 1 & anrate_pctile!=11) ~ paste0(round(tot_an_rate, 1)," hospitalizations per 10k (",anrate_pctile,"st percentile)"),
                                  (anrate_pctile %% 10 == 2 & anrate_pctile!=12) ~ paste0(round(tot_an_rate, 1)," hospitalizations per 10k (",anrate_pctile,"nd percentile)"),
                                  (anrate_pctile %% 10 == 3 & anrate_pctile!=13) ~ paste0(round(tot_an_rate, 1)," hospitalizations per 10k (",anrate_pctile,"rd percentile)"),
                                  TRUE ~ paste0(round(tot_an_rate, 1)," hospitalizations per 10k (",anrate_pctile,"th percentile)"))) %>%
  mutate(risk_group_text1 = case_when(cluster == "Dual Channel Risk" ~ "both relatively high cardiovascular disease incidence and heat exposure",
                                     cluster == "Heat Driven Risk" ~ "relatively high heat exposure",
                                     cluster == "Health Driven Risk" ~ "relatively high cardiovascular disease incidence",
                                     cluster == "Low Risk" ~ "relatively low cardiovascular disease incidence and heat exposure")) %>%
  mutate(risk_group_text2 = case_when(cluster == "Dual Channel Risk" ~ " contribute to the heat-attributable CVD burden rate ",
                                     cluster == "Heat Driven Risk" ~ " contributes to the heat-attributable CVD burden rate ",
                                     cluster == "Health Driven Risk" ~ " contributes to the heat-attributable CVD burden rate ",
                                     cluster == "Low Risk" ~ " result in a below average heat-attributable CVD burden rate ")) %>%
  mutate(health_cont_text = case_when(an_cvddiff1_rate<0 ~ "negative",
                                      an_cvddiff1_rate>0 ~ "positive")) %>%
  mutate(health_cont_sign = case_when(an_cvddiff1_rate<0 ~ "≤",
                                      an_cvddiff1_rate>0 ~ ">")) %>%
  mutate(heat_cont_text = case_when(an_tempdiff2_rate<0 ~ "negative",
                                    an_tempdiff2_rate>0 ~ "positive")) %>%
  mutate(heat_cont_sign = case_when(an_tempdiff2_rate<0 ~ "≤",
                                    an_tempdiff2_rate>0 ~ ">")) %>%
  mutate(overall_anr_text = case_when(tot_an_rate < 5.788092 ~ "lower",
                                      tot_an_rate > 5.788092 ~ "higher"))

# Waterfall Plot Data - MAIN DF FOR WATERFALL PLOTS! 
waterfall_df <- cbg_comp_sf_plot %>% 
  select(GEOID, an_atac_rate, an_gtgc_rate, an_cvddiff1_rate, an_tempdiff2_rate) %>%
  mutate(start_1 = 0,
         start_2 = an_atac_rate,
         start_3 = an_atac_rate+an_cvddiff1_rate,
         start_4 = 0,
         amount_1 = an_atac_rate,
         amount_2 = an_cvddiff1_rate,
         amount_3 = an_tempdiff2_rate,
         amount_4 = an_gtgc_rate) %>%
  mutate(end_1 = start_1+amount_1,
         end_2 = start_2+amount_2,
         end_3 = start_3+amount_3,
         end_4 = start_4+amount_4) %>%
  pivot_longer(start_1:end_4, names_to="group", values_to="rate") %>%
  separate(group, into=c("var","group"), remove=T) %>%
  pivot_wider(names_from=var, values_from=rate) %>%
  group_by(GEOID) %>%
  mutate(id = row_number()) %>%
  ungroup() %>%
  # Fields for plotting
  mutate(group = factor(group, levels=c(1,2,3,4), 
                        labels=c("RT Summer Avg. Rate","Health Contribution","Heat Contribution", "CBG Attributable Rate"))) %>%
  mutate(label_loc = pmax(start,end)) %>%
  mutate(color_assignment = case_when(group == "RT Summer Avg. Rate" ~ "Average",
                                      group == "CBG Attributable Rate" ~ "CBG",
                                      (group %in% c("Health Contribution","Heat Contribution") & amount> 2) ~ "Large Above",
                                      (group %in% c("Health Contribution","Heat Contribution") & amount< 2 & amount>0) ~ "Small Above",
                                      (group %in% c("Health Contribution","Heat Contribution") & amount< -2) ~ "Large Below",
                                      (group %in% c("Health Contribution","Heat Contribution") & amount> -2 & amount<0) ~ "Small Below"))


### Demographics Plots
cbg_pop_totals <- read_csv("../Data/Analysis/CBG_Subgroup_Pop_Totals.csv") %>%
  mutate(GEOID = as.character(GEOID))

# Demographic statistics aggregated across the whole RTP area
demobar_allrtp <- cbg_comp_sf_plot %>%
  left_join(cbg_pop_totals) %>%
  st_drop_geometry() %>%
  mutate(area = "All Research Triangle") %>%
  group_by(area) %>%
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
  select(area, pct_male, pct_female, pct_asian, pct_black,
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
  mutate(pct_round = round(pct*100,1)) %>%
  mutate(GEOID = area)


demobar_df <- cbg_comp_sf_plot %>%
  left_join(cbg_pop_totals) %>%
  st_drop_geometry() %>%
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
  select(GEOID, pct_male, pct_female, pct_asian, pct_black,
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
  mutate(pct_round = round(pct*100,1)) %>%
  mutate(area = paste0("CBG")) %>%
  # BRING IN AGGREGATED
  bind_rows(demobar_allrtp) 

# Color Palettes
pal <- wes_palette("Zissou1", n = 5) 
pal2 <- wes_palette("Darjeeling1", n = 5) 

# green_cont_pal <- colorRampPalette(c('#84EFD8', '#034036'),100)

green_cont_pal <- colorRampPalette(c('#d2f9f1', "#00A08A", "#006E5F", '#003930'))(4) # colorRampPalette(c('#d2f9f1', "#00A08A", "#0B7E6D", "#006E5F", '#024F43', '#003930'))(6)

cont_pal <- colorNumeric(
  palette = green_cont_pal, # "inferno"
  domain = range(map_df$tot_an_rate)
)

factpal <- colorFactor(c(pal[5],pal2[2], "#fdae61", lighten(pal2[5],0.5)), map_df$cluster)

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### 2. CREATE FUNCTIONS ----
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### WATERFALL PLOT ####

# Decomposition Waterfall Selection Function
plot_waterfall <- function(data, cbg_geoid) {
  if(nrow(data)==0) {
    return(NULL)
  }
  
  ggplot() + 
    geom_rect(data=data[data$GEOID==cbg_geoid,], mapping=aes(x=group, xmin = id - 0.45, xmax = id + 0.45, ymin = end, ymax = start,
                                                                                  colour=color_assignment, fill=color_assignment), alpha=0.75) +
    geom_text(data=data[data$GEOID==cbg_geoid,], aes(x=group,y=label_loc+1.5, label=round(amount,2)),vjust=.5, size=6) +
    geom_segment(data=data[data$GEOID==cbg_geoid & data$id!=4,], 
                 aes(x=id+.45, xend=id+.55, y=end, yend=end),linetype="twodash",linewidth=0.75) +
    scale_colour_manual(values=c("Average"="gray50",
                                 "CBG" = pal2[2],
                                 "Large Above" = pal[5],
                                 "Small Above" = lighten(pal[5],0.5),
                                 "Small Below" = lighten(pal[1],0.5),
                                 "Large Below" = pal[1]),
                        aesthetics=c("colour","fill")) +
    ylim(c(0,24)) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(size=12)) +
    labs(title = "Decomposition of CBG Attributable Burden Rate",
         # subtitle = paste("CBG ",cbg_geoid),
         y="",
         x="") 
}


# Default Waterfall Plot 
ggplot() +
  geom_rect(data=waterfall_df[waterfall_df$id==1,], mapping=aes(x=group, xmin = id - 0.45, xmax = id + 0.45, ymin = end, ymax = start,
                                                           colour=color_assignment, fill=color_assignment), alpha=0.75) +
  geom_text(data=waterfall_df[waterfall_df$id==1,], aes(x=group,y=label_loc+1.5, label=round(amount,2)),vjust=.5, size=6) +
  scale_colour_manual(values=c("Average"="gray50",
                               "CBG" = pal2[2],
                               "Large Above" = pal[5],
                               "Small Above" = lighten(pal[5],0.5),
                               "Small Below" = lighten(pal[1],0.5),
                               "Large Below" = pal[1]),
                      aesthetics=c("colour","fill")) +
  ylim(c(0,24)) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size=12)) +
  labs(title = "RTP Average Heat-Attributable Hospitalizations per 10,000",
       subtitle = "Summer 2018",
       y="",
       x="")

#### RISK GROUP MAP ####

# Map Plot Function
plot_map <- function(data) {
  
  leaflet() %>%
    addTiles() %>%
    # addProviderTiles(providers$CartoDB.PositronNoLabels)  %>%
    setView(lng = -78.8, lat = 35.84, zoom = 8) %>%
    addPolygons(data=data,

                # state border stroke color
                color = 'white',

                # soften the weight of the state borders
                weight = 1,

                # values >1 simplify the polygons' lines for less detail but faster loading
                smoothFactor = .3,

                # set opacity of polygons
                fillOpacity = .65,

                # specify that the each state should be colored per paletteNum()
                fillColor = ~factpal(cluster),
                
                # Popup
                label = ~content) %>% # issue is here
    addLegend("topleft", pal=factpal, values = data$cluster,
              title = "CBG Cluster",
              # labFormat = labelFormat(prefix = "$"),
              opacity = 1)
  
}

#### RATE MAP ####
   
# FUNCTION TO PLOT RATE
# data should = map_df
plot_ratemap <- function(data) {                                
leaflet() %>%
  addTiles() %>%
  # addProviderTiles(providers$CartoDB.PositronNoLabels)  %>%
  setView(lng = -78.8, lat = 35.84, zoom = 9) %>%
  addPolygons(data=data,

              # state border stroke color
              color = 'white',

              # soften the weight of the state borders
              weight = 1,

              # values >1 simplify the polygons' lines for less detail but faster loading
              smoothFactor = .3,

              # set opacity of polygons
              fillOpacity = .75,

              # specify that the each state should be colored per paletteNum()
              fillColor = ~cont_pal(tot_an_rate),

              # popup - would like to add county and city here.
              label = ~content) %>%
  addLegend("topleft", pal=cont_pal, values = data$tot_an_rate,
            title = "CBG Heat-Attributable Hosp. per 10k",
            opacity = 1)
}

#### DEMOGRAPHICS BAR PLOT ####

plot_demobar <- function(data, cbg_geoid){
  if(nrow(data)==0) {
    return(NULL)
  }
  
  # SEX
  sex_comp <- ggplot(data[data$subgroup_type=="Sex" & data$GEOID %in% c("All Research Triangle",cbg_geoid),],
                     aes(x=var_label, y=pct_round, fill=area)) + #
    geom_bar(position="dodge",stat="identity",color="black") +
    geom_text(aes(label=format(pct_round, nsmall=1)), vjust=-0.25, position=position_dodge(width=0.9), size=4.5) +
    coord_cartesian(clip="off") +
    labs(y= "Percent", x = "Sex", fill="Geography") + # , fill = "CBG Attr. Rate Group"
    scale_y_continuous(limits=c(0,100), n.breaks=5) +
    scale_fill_manual(values = c("All Research Triangle"="gray50", "CBG"=pal2[2])) + # "Top 10%"=pal[5], "Bottom 90%"=lighten(pal2[5],0.5)
    theme_minimal() +
    theme(panel.grid = element_blank(),
          text = element_text(size = 15),
          axis.text.x = element_text(size=12))
  
  # Race
  race_comp <- ggplot(data[data$subgroup_type=="Race" & data$var!="pct_native" & data$GEOID %in% c("All Research Triangle",cbg_geoid),],
                      aes(x=var_label, y=pct_round, fill=area)) + # exclude native bc not in model
    geom_bar(position="dodge",stat="identity",color="black") +
    geom_text(aes(label=format(pct_round, nsmall=1)), vjust=-0.25, position=position_dodge(width=0.9), size=4.5) +
    coord_cartesian(clip="off") +
    labs(y= "Percent", x = "Race", fill="Geography") +
    scale_y_continuous(limits=c(0,100), n.breaks=5) +
    scale_fill_manual(values = c("All Research Triangle"="gray50", "CBG"=pal2[2])) + # "Top 10%"=pal[5], "Bottom 90%"=lighten(pal2[5],0.5)
    theme_minimal() +
    theme(panel.grid = element_blank(),
          text = element_text(size = 15),
          axis.text.x = element_text(size=12))
  
  # Age
  age_comp <- ggplot(data[data$subgroup_type=="Age" & data$GEOID %in% c("All Research Triangle",cbg_geoid),],
                     aes(x=var_label, y=pct_round, fill=area)) + 
    geom_bar(position="dodge",stat="identity",color="black") +
    geom_text(aes(label=format(pct_round, nsmall=1)), vjust=-0.25, position=position_dodge(width=0.9), size=4.5) +
    coord_cartesian(clip="off") +
    labs(y= "Percent", x = "Age", fill="Geography") +
    scale_y_continuous(limits=c(0,100), n.breaks=5) +
    scale_fill_manual(values = c("All Research Triangle"="gray50", "CBG"=pal2[2])) + # "Top 10%"=pal[5], "Bottom 90%"=lighten(pal2[5],0.5)
    theme_minimal() +
    theme(panel.grid = element_blank(),
          text = element_text(size = 15),
          axis.text.x = element_text(size=12))
  
  # Poverty + Education
  other_comp <- ggplot(data[data$subgroup_type=="Poverty + Education" & data$GEOID %in% c("All Research Triangle",cbg_geoid),],
                       aes(x=var_label, y=pct_round, fill=area)) + #
    geom_bar(position="dodge",stat="identity",color="black") +
    geom_text(aes(label=format(pct_round, nsmall=1)), vjust=-0.25, position=position_dodge(width=0.9), size=4.5) +
    coord_cartesian(clip="off") +
    labs(y= "Percent", x = "Poverty + Education", fill="Geography") +
    scale_y_continuous(limits=c(0,100), n.breaks=5) +
    scale_fill_manual(values = c("All Research Triangle"="gray50", "CBG"=pal2[2])) + # "Top 10%"=pal[5], "Bottom 90%"=lighten(pal2[5],0.5)
    theme_minimal() +
    theme(panel.grid = element_blank(),
          text = element_text(size = 15),
          axis.text.x = element_text(size=12))
  
  bar_comp <- ggarrange(sex_comp, race_comp, age_comp, other_comp,
                        nrow=2, ncol=2, common.legend = T) # sooo similar
  bar_comp
}


# Default Demographics Bar Plot
# SEX
sex_comp_default <- ggplot(demobar_allrtp[demobar_allrtp$subgroup_type=="Sex",],
                   aes(x=var_label, y=pct_round, fill=area)) + #
  geom_bar(position="dodge",stat="identity",color="black") +
  geom_text(aes(label=format(pct_round, nsmall=1)), vjust=-0.25, position=position_dodge(width=0.9), size=4.5) +
  coord_cartesian(clip="off") +
  labs(y= "Percent", x = "Sex", fill = "Geography") + # , fill = "CBG Attr. Rate Group"
  scale_y_continuous(limits=c(0,100), n.breaks=5) +
  scale_fill_manual(values = c("All Research Triangle"="gray50")) + # "Top 10%"=pal[5], "Bottom 90%"=lighten(pal2[5],0.5)
  theme_minimal() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15),
        axis.text.x = element_text(size=12))

# Race
race_comp_default <- ggplot(demobar_allrtp[demobar_allrtp$subgroup_type=="Race" & demobar_allrtp$var!="pct_native",],
                    aes(x=var_label, y=pct_round, fill=area)) + # exclude native bc not in model
  geom_bar(position="dodge",stat="identity",color="black") +
  geom_text(aes(label=format(pct_round, nsmall=1)), vjust=-0.25, position=position_dodge(width=0.9), size=4.5) +
  coord_cartesian(clip="off") +
  labs(y= "Percent", x = "Race", fill = "Geography") +
  scale_y_continuous(limits=c(0,100), n.breaks=5) +
  scale_fill_manual(values = c("All Research Triangle"="gray50")) + # "Top 10%"=pal[5], "Bottom 90%"=lighten(pal2[5],0.5)
  theme_minimal() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15),
        axis.text.x = element_text(size=12))

# Age
age_comp_default <- ggplot(demobar_allrtp[demobar_allrtp$subgroup_type=="Age",],
                   aes(x=var_label, y=pct_round, fill=area)) + #
  geom_bar(position="dodge",stat="identity",color="black") +
  geom_text(aes(label=format(pct_round, nsmall=1)), vjust=-0.25, position=position_dodge(width=0.9), size=4.5) +
  coord_cartesian(clip="off") +
  labs(y= "Percent", x = "Age", fill = "Geography") +
  scale_y_continuous(limits=c(0,100), n.breaks=5) +
  scale_fill_manual(values = c("All Research Triangle"="gray50")) + # "Top 10%"=pal[5], "Bottom 90%"=lighten(pal2[5],0.5)
  theme_minimal() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15),
        axis.text.x = element_text(size=12))

# Poverty + Education
other_comp_default <- ggplot(demobar_allrtp[demobar_allrtp$subgroup_type=="Poverty + Education",],
                     aes(x=var_label, y=pct_round, fill=area)) + #
  geom_bar(position="dodge",stat="identity",color="black") +
  geom_text(aes(label=format(pct_round, nsmall=1)), vjust=-0.25, position=position_dodge(width=0.9), size=4.5) +
  coord_cartesian(clip="off") +
  labs(y= "Percent", x = "Poverty + Education", fill = "Geography") +
  scale_y_continuous(limits=c(0,100), n.breaks=5) +
  scale_fill_manual(values = c("All Research Triangle"="gray50")) + # "Top 10%"=pal[5], "Bottom 90%"=lighten(pal2[5],0.5)
  theme_minimal() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 15),
        axis.text.x = element_text(size=12))

demobar_default <- ggarrange(sex_comp_default, race_comp_default, age_comp_default, other_comp_default,
                      nrow=2, ncol=2, common.legend = T) # sooo similar
demobar_default
rm(sex_comp_default, age_comp_default, race_comp_default, other_comp_default)

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### 3. Set Up Dashboard ----
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### UI ####

my_theme = create_theme(
  adminlte_color(
    light_blue = "#78B7C5"
  )
)

# LAYOUT
ui<- dashboardPage(
  # Header
  dashboardHeader(disable=TRUE), # width doesn't work

  # # Sidebar
  dashboardSidebar(
    disable = TRUE
    # sidebarMenu(
    #   menuItem("Research Triangle CBG Map",tabName = 'map', icon = icon("location-dot"))
    # )
  ),
  # Body
  dashboardBody(

    use_theme(my_theme), # custom skin color defined in my_theme
    # tabItems(
    ######################### Map #
    # tabItem(tabName = "map",
    h2("CBG Heat-Related CVD Burdens in North Carolina's Research Triangle Area During Summer 2018"),
    # fluidPage(
    fluidRow(
      column(width = 6,
             # TESTING - input dropdown
             selectInput(
               inputId = "select_map",
               label = "Select Map Value",
               choices = c(
                 "Risk Groups",
                 "Heat-Attributable CVD Rates"
               ),
               selected = "Counts",
               multiple = FALSE
             ),
             # Map on Upper Left
             box(width=NULL, height=830,
                 leafletOutput("mymap",height=800),
                 status='primary',
                 headerBorder =FALSE)


      ),
      column(width=6,

             # Waterfall Plot Middle Right
             box(width=NULL, height=480, # height=290
                 htmlOutput('ui_waterfall'),
                 plotOutput('waterfall',height=260),
                 status='primary',
                 headerBorder =FALSE),
             
             # Demographics Bar Plot Bottom Right
             box(width=NULL, height=410,
                 title = "Demographic Characteristics of Selected Census Block Group",
                 plotOutput('demobar',height=340),
                 status='primary',
                 headerBorder =FALSE)
      )
    )
  )
)

# fluidRow within fluidPage - adds up to 12 so pick width to add up to 12

#### SERVER ####

server <- shinyServer(function(input, output) {
  
  #### initialize reactive values ####
  rvs <- reactiveValues(map=NULL, poly_cbg = map_df)
  
  rv_shape <- reactiveVal(FALSE)
  rv_location <- reactiveValues(id=NULL,lat=NULL,lng=NULL)
  rv_location_move_old <- reactiveValues(lat=NULL,lng=NULL)
  rv_text <- reactiveValues(click='Click on the map to see census block group information.')
  
  #### output ####
  
  ## output: leaflet map ## 
  output$mymap <- renderLeaflet({
    if (input$select_map == "Heat-Attributable CVD Rates") {
      # Rate Map
      leaflet() %>%
        addTiles() %>%
        # addProviderTiles(providers$CartoDB.PositronNoLabels)  %>%
        setView(lng = -78.8, lat = 35.84, zoom = 9) %>%
        addPolygons(data=map_df,
                    
                    # state border stroke color
                    color = 'white',
                    
                    # soften the weight of the state borders
                    weight = 1,
                    
                    # values >1 simplify the polygons' lines for less detail but faster loading
                    smoothFactor = .3,
                    
                    # set opacity of polygons
                    fillOpacity = .75,
                    
                    # specify that the each state should be colored per paletteNum()
                    fillColor = ~cont_pal(tot_an_rate),
                    
                    # Specify ID for reactive values to use
                    layerId = rvs$poly_cbg$GEOID,
                    
                    # popup - would like to add county and city here.
                    label = ~content) %>%
        addLegend("topleft", pal=cont_pal, values = map_df$tot_an_rate,
                  title = "CBG Heat-Attributable Hosp. per 10k",
                  opacity = 1)
    }
    else{
    # if (input$select_map == "Risk Groups") {
      # Risk Group Map
      leaflet() %>% 
        addTiles() %>% # color background map
        # addProviderTiles(providers$CartoDB.Positron) %>% # black and white background map
        setView(lng = -78.8, lat = 35.84, zoom = 8.5) %>%
        addPolygons(data=map_df,
                    
                    # state border stroke color
                    color = 'white',
                    
                    # soften the weight of the state borders
                    weight = 1,
                    
                    # values >1 simplify the polygons' lines for less detail but faster loading
                    smoothFactor = .3,
                    
                    # set opacity of polygons
                    fillOpacity = .7,
                    
                    # specify that the each state should be colored per paletteNum()
                    fillColor = ~factpal(cluster),
                    
                    # Specify ID for reactive values to use
                    layerId = rvs$poly_cbg$GEOID,
                    
                    # Popup
                    label = ~content) %>% 
        addLegend("topleft", pal=factpal, values = map_df$cluster,
                  title = "CBG Cluster",
                  opacity = 1)
      
    }
  })
  
  ## output: Print Information on CBG ## 
  
  # TEXT ABOVE WATERFALL PLOT
  output$ui_waterfall <- renderUI({
    location_info <- reactiveValuesToList(rv_location)
    if (!all(is.null(unlist(location_info)))){ # if any entry in rv_location is not NULL
      HTML(paste(h4(strong('Results for Census Block Group:',rv_location$id)),
                 
                 h4('Overall, the attributable burden rate for this CBG is',
                    rvs$poly_cbg$anr_pct_text[rvs$poly_cbg$GEOID==rv_location$id],
                    'which is ',
                    strong(rvs$poly_cbg$overall_anr_text[rvs$poly_cbg$GEOID==rv_location$id]), 
                    ' than the Research Triangle average.'),
                 
                 h4('This census block group is in the ', 
                    strong(rvs$poly_cbg$cluster[rvs$poly_cbg$GEOID==rv_location$id]),
                    'Group, meaning that the Health Contribution is ', 
                    rvs$poly_cbg$health_cont_text[rvs$poly_cbg$GEOID==rv_location$id],
                    '(',
                    rvs$poly_cbg$health_cont_sign[rvs$poly_cbg$GEOID==rv_location$id],
                    '0) and the Heat-Contribution is ',
                    rvs$poly_cbg$heat_cont_text[rvs$poly_cbg$GEOID==rv_location$id],
                    '(',
                    rvs$poly_cbg$heat_cont_sign[rvs$poly_cbg$GEOID==rv_location$id],
                    '0).',
                    'This indicates that',
                    strong(rvs$poly_cbg$risk_group_text1[rvs$poly_cbg$GEOID==rv_location$id]),
                    rvs$poly_cbg$risk_group_text2[rvs$poly_cbg$GEOID==rv_location$id],
                    'in this CBG.')
      ))
      # HTML(paste(h4(strong('Results for Census Block Group:',rv_location$id)),
      #            
      #            h4('This census block group is in the ', 
      #            strong(rvs$poly_cbg$cluster[rvs$poly_cbg$GEOID==rv_location$id]),
      #            'Group, meaning that the contribution of local CVD incidence is positive', 
      #            strong(rvs$poly_cbg$health_cont_text[rvs$poly_cbg$GEOID==rv_location$id]), 
      #            'than the Research Triangle area average (Health Contribution ',
      #            rvs$poly_cbg$health_cont_sign[rvs$poly_cbg$GEOID==rv_location$id],
      #            ' 0) and the summer average daily temperature is ',
      #            strong(rvs$poly_cbg$heat_cont_text[rvs$poly_cbg$GEOID==rv_location$id]), 
      #            'than the Research Triangle average (Heat Contribution ',
      #            rvs$poly_cbg$heat_cont_sign[rvs$poly_cbg$GEOID==rv_location$id],
      #            ' 0).'),
      #            
      #            h4('Overall, the attributable burden rate for this CBG is',
      #               rvs$poly_cbg$anr_pct_text[rvs$poly_cbg$GEOID==rv_location$id],
      #               'which is ',
      #               strong(rvs$poly_cbg$overall_anr_text[rvs$poly_cbg$GEOID==rv_location$id]), 
      #               ' than the Research Triangle average. 
      #               At this location, the CVD incidence rate was ',
      #               rvs$poly_cbg$cvd_text[rvs$poly_cbg$GEOID==rv_location$id],
      #               'and the daily average temperature during the summer of 2018
      #               was ',
      #               rvs$poly_cbg$temp_text[rvs$poly_cbg$GEOID==rv_location$id],
      #               '.')
      # ))
    }
    else{
      HTML(paste(h4('Click on the map to see census block group information.')))
    }
  })

  ## output: Waterfall Plot ##
  output$waterfall <- renderPlot ({
    # Before CBG is selected - just show area average
    if (!rv_shape()){
      ggplot() + 
        geom_rect(data=waterfall_df[waterfall_df$id==1,], mapping=aes(x=group, xmin = id - 0.45, xmax = id + 0.45, ymin = end, ymax = start,
                                                                      colour=color_assignment, fill=color_assignment), alpha=0.75) +
        geom_text(data=waterfall_df[waterfall_df$id==1,], aes(x=group,y=label_loc+1.5, label=round(amount,2)),vjust=.5, size=6) +
        scale_colour_manual(values=c("Average"="gray50",
                                     "CBG" = pal2[2],
                                     "Large Above" = pal[5],
                                     "Small Above" = lighten(pal[5],0.5),
                                     "Small Below" = lighten(pal[1],0.5),
                                     "Large Below" = pal[1]),
                            aesthetics=c("colour","fill")) +
        ylim(c(0,24)) +
        theme_minimal() +
        theme(legend.position = "none",
              axis.text.x = element_text(size=12)) +
        labs(title = "Research Triangle Average Heat-Attributable Hospitalizations per 10,000",
             subtitle = "Summer 2018",
             y="",
             x="")
      
    # Once CBG is selected on map
    }else{
      print(rv_location$id)
      plot_waterfall(data=waterfall_df, cbg_geoid = rv_location$id)
    }
  })
  
  ## output: Demographic Bar Plot ##
  # Before CBG is selected - just show area average
  output$demobar <- renderPlot ({
    if (!rv_shape()){
      demobar_default

    # Once CBG is selected on map
    }else{
      print(rv_location$id)
      plot_demobar(data=demobar_df, cbg_geoid = rv_location$id)
    }
  })

  #### observe mouse events ####
  
  ## when any click happens
  observeEvent(input$mymap_click,{
    mymap_shape_click_info <- input$mymap_shape_click
    mymap_click_info <- input$mymap_click
    print(mymap_shape_click_info)
    rv_shape(TRUE)
    
    rv_location$id <- mymap_shape_click_info$id
    rv_location$lat <- round(mymap_click_info$lat, 4)
    rv_location$lng <- round(mymap_click_info$lng, 4)
  })

  
}
)

#### VIEW APP ####
shinyApp(ui = ui, server = server)

