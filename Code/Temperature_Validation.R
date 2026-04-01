################################################################################
# Program Name: Temperature_Validation.R
# Program Purpose: Assess similarity of multiple near-surface air temperature data sets
# Compare aggregated, model-based ZCTA level estimates vs. direct ZCTA-level data from Medicare

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
  library(sjPlot)
  library(nngeo)
  library(fields) # thin plate spline regression
}

## Set Up ----

zctas_shp <- st_read("../Data/Original/Geography/RTP_CBSAZCTA.shp")
zctacrs <- st_crs(zctas_shp)

## 2. Prep WS Data (code from Format_Heat_Data_CBG) ----

ws_files <- list.files(path="../Data/Original/Temperature/Weather_Station/GSSOD_2018", full.names = TRUE)

ws_all <- do.call(rbind,lapply(ws_files, function(i){read_csv(i)})) 

ws_clean <- ws_all %>%
  mutate(MAX = ifelse(MAX == 9999.9, NA, MAX),
         MIN = ifelse(MIN == 9999.9, NA, MIN)) %>%
  filter(TEMP_ATTRIBUTES >=12) %>% # require stations to have at least 12 hours of observations - new from 3/15 version
  # mutate(TEMP = ifelse(TEMP == 5.5, NA, TEMP)) %>%
  rename(tmean_f_ws = TEMP) %>%
  mutate(Station_ID = as.character(STATION)) %>%
  rename(date = DATE) %>%
  mutate(none_missing = TEMP_ATTRIBUTES==24) %>%
  select(Station_ID, date, tmean_f_ws, none_missing) 

# Get Station List
gssod_stations <- read_excel("../Data/Original/Temperature/Weather_Station/GSSOD_Stations_2018.xlsx") %>%
  distinct(Station_ID) %>%
  filter(!is.na(Station_ID))

# All Weather Station Data from Historical Observing Metadata Repository (https://www.ncei.noaa.gov/access/homr/reports)
isd_history <- read_excel("../Data/Original/Temperature/Weather_Station/isd_history_txt.xlsx") %>%
  mutate(USAF = str_pad(USAF, width=6, side="left", pad="0"),
         WBAN = str_pad(WBAN, width=5, side="left", pad="0")) %>%
  mutate(Station_ID = paste(USAF, WBAN, sep=""))

# Create set of viable weather stations
station_points <- isd_history %>%
  
  # Only those with GSSOD files for all of 2018
  filter(Station_ID %in% unique(gssod_stations$Station_ID)) %>%
  filter(!is.na(LON) & !is.na(LAT)) %>% # Only 1 missing in AK
  filter(END>=20190101) %>% # Horace Williams goes offline mid-2018
  
  # Create Geometry
  st_as_sf(coords = c("LON", "LAT"), crs = zctacrs) %>%
  
  # Will use row index to match back to the closest stations
  rowid_to_column("row_index") %>%
  
  # FILTER to stations in area
  filter(Station_ID %in% unique(ws_clean$Station_ID)) %>%
  select(Station_ID, geometry)

# st_write(station_points, "../Data/Original/Temperature/Weather_Station/Station_Points.shp")

stations <- unique(station_points$Station_ID)

# Find days with missing values
dates <- unique(ws_clean$date)

ws_all_dates <- expand.grid(dates, stations) %>%
  rename(date = Var1,
         Station_ID = Var2) %>%
  left_join(ws_clean, by=c("date", "Station_ID")) %>%
  mutate(flag = ifelse(is.na(tmean_f_ws), 1, 0))

area_avg <- ws_all_dates %>% 
  filter(!is.na(tmean_f_ws)) %>%
  group_by(date) %>%
  summarise(tmean_avg = mean(tmean_f_ws))

# Filled
ws_final <- ws_all_dates %>%
  left_join(area_avg, by="date") %>%
  # Replace missing values with area average for daily temp
  mutate(tmean_f_ws_nn = ifelse(flag==1, tmean_avg, tmean_f_ws)) %>% # nn is filled
  select(-c(tmean_avg))

rm(area_avg, ws_all_dates, ws_clean, ws_all, ws_files)

## 3. Prep 1km Data at WS Points (code adapted from Format_Heat_Data_CBG) ----

# Files available on Google Drive: https://drive.google.com/drive/folders/1w8cVGYO3_wturLneBdq0W-1fOvc9J2m2?usp=sharing

tmax_2018_ws <- read_csv("G:/My Drive/Global_Daily_Air_Temp_WS/tmax2018_rtp_ws.csv") %>%
  rename(system.index = `system:index`)

nmax <- max(str_count(tmax_2018_ws$system.index, "_")) + 1

tmax <- tmax_2018_ws %>%
  separate(col=system.index, into = paste0("var", seq_len(nmax)), sep="_") %>%
  mutate(year = substr(var4, 1, 4)) %>%
  mutate(date = as.Date(as.numeric(var5)-1, origin=paste0(year, "-01-01"))) %>% # it is accounting for leap years
  rename(tmax_c = b1, # , # first
         doy = var5) %>%
  select(Station_ID, year, date, doy, tmax_c)

tmin_2018_ws <- read_csv("G:/My Drive/Global_Daily_Air_Temp_WS/tmin2018_rtp_ws.csv") %>%
  rename(system.index = `system:index`)

nmax <- max(str_count(tmin_2018_ws$system.index, "_")) + 1

tmin <- tmin_2018_ws %>%
  separate(col=system.index, into = paste0("var", seq_len(nmax)), sep="_") %>%
  mutate(year = substr(var4, 1, 4)) %>%
  mutate(date = as.Date(as.numeric(var5)-1, origin=paste0(year, "-01-01"))) %>% # it is accounting for leap years
  rename(tmin_c = b1, # , first
         doy = var5) %>%
  select(Station_ID, year, date, doy, tmin_c)

ws_1km <- tmax %>%
  left_join(tmin, by=c("Station_ID", "year", "date", "doy")) %>%
  mutate(tmean_c_1km = (tmin_c + tmax_c)/2) %>%
  mutate(tmean_f_1km = ((tmean_c_1km*(9/5)+32))) %>%
  select(Station_ID, date, tmean_c_1km, tmean_f_1km) %>%
  mutate(Station_ID = as.character(Station_ID))

rm(tmax_2018_ws, tmin_2018_ws, tmax, tmin)

## 4. Combine for Validation ----

stations_all <- ws_final %>%
  left_join(ws_1km, by=c("Station_ID","date")) %>%
  mutate(diff_mean = tmean_f_ws - tmean_f_1km) %>%
  mutate(tmean_c_ws = (tmean_f_ws-32)*(5/9)) %>%
  mutate(diff_mean_c = tmean_c_1km - tmean_c_ws) %>% # 1km - WS
  mutate(summer = case_when(date >= mdy("05-01-2018") & date < mdy("10-01-2018") ~ 1,
                            TRUE ~ 0)) %>%
  mutate(all = 1) # %>%
  # filter(summer==1)
  
hist(stations_all$diff_mean_c)
hist(stations_all$diff_mean_c[stations_all$none_missing==TRUE])

summary(stations_all$diff_mean_c)
summary(stations_all$diff_mean_c[stations_all$none_missing==TRUE])
summary(stations_all$diff_mean_c[stations_all$summer==1])

# Summary Statistics
overall <- stations_all %>%
  group_by(all) %>%
  summarise(mean_diff = mean(diff_mean_c, na.rm=T),
            sd_diff = sd(diff_mean_c, na.rm=T),
            n = n()) %>%
  ungroup() %>%
  mutate(group = "all") %>%
  select(-all)

nomiss <- stations_all %>%
  filter(none_missing == TRUE) %>%
  group_by(none_missing) %>%
  summarise(mean_diff = mean(diff_mean_c, na.rm=T),
            sd_diff = sd(diff_mean_c, na.rm=T),
            n = n()) %>%
  ungroup() %>%
  mutate(group = "none_missing") %>%
  select(-none_missing)

summer <- stations_all %>%
  filter(summer == TRUE) %>%
  group_by(summer) %>%
  summarise(mean_diff = mean(diff_mean_c, na.rm=T),
            sd_diff = sd(diff_mean_c, na.rm=T),
            n = n()) %>%
  ungroup() %>%
  mutate(group = "summer") %>%
  select(-summer)

summer_nomiss <- stations_all %>%
  filter(summer == TRUE & none_missing==TRUE) %>%
  group_by(summer, none_missing) %>%
  summarise(mean_diff = mean(diff_mean_c, na.rm=T),
            sd_diff = sd(diff_mean_c, na.rm=T),
            n = n()) %>%
  ungroup() %>%
  mutate(group = "summer_nomiss") %>%
  select(-c(summer, none_missing))

summary_stats <- overall %>%
  bind_rows(nomiss) %>%
  bind_rows(summer) %>%
  bind_rows(summer_nomiss)

# Regression
reg <- lm(tmean_c_ws ~ tmean_c_1km, stations_all)
summary(reg)


# Plots
summary_by_station <- stations_all %>%
  group_by(Station_ID) %>%
  summarise(min = min(diff_mean_c, na.rm=T),
            p25 = quantile(diff_mean_c, probs = .25, na.rm=T),
            median = median(diff_mean_c, na.rm=T),
            mean = mean(diff_mean_c, na.rm=T),
            p75 = quantile(diff_mean_c, probs = .75, na.rm=T),
            max = max(diff_mean_c, na.rm=T))
# write.csv(summary_by_station, "../Data/Figures/Temperature Validation/WSvsNS_By_Station.csv")

station_plot_data <- station_points %>% 
  left_join(summary_by_station) %>%
  select(Station_ID, mean) %>%
  mutate(mean_lbl = round(mean, digits=2))

ggplot() +
  geom_sf(data=zctas_shp, color="white",fill="grey") +
  geom_sf(data=station_plot_data) +
  geom_sf_text(data=station_plot_data, aes(label=mean_lbl), nudge_y = 0.05) +
  theme_minimal() +
  ggtitle("Mean Difference By Weather Station (Gridded 1km - Weather Station)")

ggsave("../Data/Figures/Temperature Validation/WS_Diff_Map.png",dpi=300, width=8,height=8)

## 5. Validation figures  ----

tempDiff_long <- stations_all %>%
  select(Station_ID, date, tmean_c_1km, tmean_c_ws, diff_mean_c) %>%
  group_by(date) %>%
  summarise(across(c(tmean_c_1km, tmean_c_ws, diff_mean_c), mean, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(YEAR=substr(date, 1, 4)) %>%
  mutate(date = as.Date(str_remove_all(date, "-"), format="%Y%m%d")) %>%
  gather(key = "Source", value = "Temperature", -date, -YEAR) 

ggplot(tempDiff_long, aes(x = date, y = Temperature, color = Source, group=Source, linetype=Source)) +
  geom_line() +
  scale_color_manual(values=c("tmean_c_ws"="#1bc4e7", "tmean_c_1km"="#a832c5", "diff_mean_c"="#0e0d01ff"), labels=c("tmean_c_ws"="WeatherStation", "tmean_c_1km"="1km Gridded (Zhang et al. 2022)", "diff_mean_c"="Differences")) + 
  scale_linetype_manual(values=c("tmean_c_ws"="solid", "tmean_c_1km"="dashed", "diff_mean_c"="solid"), labels=c("tmean_c_ws"="WeatherStation", "tmean_c_1km"="1km Gridded (Zhang et al. 2022)", "diff_mean_c"="Differences")) +
  facet_wrap(~YEAR, scales = "free_x") +
  labs(title = "Average Temperature Comparison by Source", x = "Date", y = "Temperature (°C)") +
  theme_minimal() +
  theme(axis.text.x = element_text(size=12, angle = 90, hjust = 1),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12, face = "bold")) +
  guides(color = guide_legend(override.aes = list(linetype = c("solid", "dashed", "solid"))))

ggsave("../Data/Figures/Temperature Validation/WeatherStation_vs_NearSurface.png",dpi=300, width=12,height=8)

# All observations
ggplot(stations_all) +
  geom_histogram(aes(x = diff_mean_c), bins = 30, color="grey") +
  labs(title = "Distribution of Temperature Differences (Gridded 1km - Weather Station)", x = "Temperature Difference (°C)", y = "Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=12))

ggsave("Data/Figures/Temperature Validation/Hist_WS_vs_NS.png",dpi=300, width=12,height=8)

# WS None Missing
ggplot(stations_all[stations_all$none_missing==TRUE,]) +
  geom_histogram(aes(x = diff_mean), bins = 30, color="grey") +
  labs(title = "Distribution of Temperature Differences (WeatherStation - NearSurface)", x = "Temperature Difference (°F)", y = "Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=12))

# Both - they don't look that different
ggplot(stations_all) +
  geom_histogram(aes(x = diff_mean), bins = 30, color="grey") +
  geom_histogram(data = stations_all[stations_all$none_missing==TRUE,], aes(x = diff_mean), bins = 30, color="red") +
  
  labs(title = "Distribution of Temperature Differences (WeatherStation - NearSurface)", x = "Temperature Difference (°F)", y = "Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=12))


## 6. PRISM vs. Zhang et al. 2022  ----

# Load temperature data from near surface data - change path as needed
NearSurface <- read.csv("../Data/Original/Temperature/RTP_ZCTA_Temp_2010_2019.csv")

# Load temperature data from PRISM 800m - change path as needed
RTP_PRISM <- read.csv("../Data/Original/Temperature/RTP_PRISM_TMEAN.csv", sep="\t", header=T)

PRISM_NearSurface <- RTP_PRISM %>%
  right_join(NearSurface, by = c("GEOID10", "date")) 

PRISM_NearSurface <- PRISM_NearSurface %>%
  mutate(diff_TMEAN = tmean_c - TMEAN_PRISM) # 1km - Prism

# Summary Statistics
overall_prism <- PRISM_NearSurface %>%
  mutate(group="all") %>%
  group_by(group) %>%
  summarise(mean_diff = mean(diff_TMEAN, na.rm=T),
            sd_diff = sd(diff_TMEAN, na.rm=T),
            n = n()) 

byyear_prism <- PRISM_NearSurface %>%
  group_by(year) %>%
  summarise(mean_diff = mean(diff_TMEAN, na.rm=T),
            sd_diff = sd(diff_TMEAN, na.rm=T),
            n = n()) 

# Regression
reg_prism <- lm(TMEAN_PRISM ~ tmean_c, PRISM_NearSurface)
summary(reg_prism)

# Plots
tempDiff_long_PRISM <- PRISM_NearSurface %>%
  dplyr::select(GEOID10, date, TMEAN_PRISM, tmean_c, diff_TMEAN) %>%
  group_by(date) %>%
  summarise(across(c(TMEAN_PRISM, tmean_c, diff_TMEAN), mean, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(YEAR=substr(date, 1, 4)) %>%
  mutate(date = as.Date(str_remove_all(date, "-"), format="%Y%m%d")) %>%
  gather(key = "Source", value = "Temperature", -date, -YEAR) 

ggplot(tempDiff_long_PRISM, aes(x = date, y = Temperature, color = Source, group=Source, linetype=Source)) +
  geom_line() +
  scale_color_manual(values=c("TMEAN_PRISM"="#1bc4e7", "tmean_c"="#a832c5", "diff_TMEAN"="#0e0d01ff"), labels=c("TMEAN_PRISM"="PRISM", "tmean_c"="Zhang et al. 2022", "diff_TMEAN"="Differences")) + 
  scale_linetype_manual(values=c("TMEAN_PRISM"="solid", "tmean_c"="dashed", "diff_TMEAN"="solid"), labels=c("TMEAN_PRISM"="PRISM", "tmean_c"="Zhang et al. 2022", "diff_TMEAN"="Differences")) +
  facet_wrap(~YEAR, scales = "free_x") +
  labs(title = "Average Temperature Comparison by Source", x = "Date", y = "Temperature (°C)") +
  theme_minimal() +
  theme(axis.text.x = element_text(size=12, angle = 90, hjust = 1),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12, face = "bold")) +
  guides(color = guide_legend(override.aes = list(linetype = c("solid", "dashed", "solid"))))

ggsave("../Data/Figures/Temperature Validation/Prism_vs_NearSurface.png",dpi=300, height=8, width=12)

ggplot(PRISM_NearSurface) +
  geom_histogram(aes(x = diff_TMEAN), bins = 30, color="grey") +
  labs(title = "Distribution of Temperature Differences (Gridded 1km - PRISM 800m)", x = "Temperature Difference (°C)", y = "Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=12))

ggsave("../Data/Figures/Temperature Validation/Hist_PRISM_vs_NS.png",dpi=300, width=12,height=8)

reg_prism <- lm(tmean_c ~ TMEAN_PRISM, PRISM_NearSurface)
summary(reg_prism)

# Map
tempDiff_2018_PRISM <- PRISM_NearSurface %>%
  filter(year == 2018) %>%
  group_by(GEOID10) %>%
  summarise(diff_TMEAN = mean(diff_TMEAN, na.rm=TRUE)) %>%
  ungroup() %>%
  select(GEOID10, diff_TMEAN) %>%
  mutate(GEOID10 = as.character(GEOID10))

zctas_prismvs1km <- zctas_shp %>%
  left_join(tempDiff_2018_PRISM, by="GEOID10") %>%
  mutate(mean_lbl = round(diff_TMEAN, digits=2))

ggplot() +
  geom_sf(data=zctas_prismvs1km, aes(fill=diff_TMEAN), color="white") +
  geom_sf_text(data=zctas_prismvs1km, aes(label=mean_lbl), color="black", size = 2.5) + # , nudge_y = 0.05
  scale_fill_gradient2(midpoint=0) +
  theme_minimal()  +
  ggtitle("2018 Mean Difference By ZCTA (Gridded 1km - Weather Station)")

ggsave("../Data/Figures/Temperature Validation/PRISM_Diff_Map.png",dpi=300, width=8,height=8)
