################################################################################
# Program Name: Format_Heat_Data_CBG.R
# Program Purpose: Format 2 soures of daily air temp data by CBG for 2018 only
# This data will be used in Step 3

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

# Files too large for GitHub, available on Google Drive: https://drive.google.com/drive/folders/12SBHV9RcsHHl3Pr5CyY28woV4Wm52dZu?usp=drive_link

################################################################################
# Format 1km Zhang et al. data ######

tmax_2018_cbg <- read_csv("G:/My Drive/Global_Daily_Air_Temp_CBG/tmax2018_rtp_cbg.csv") %>%
  rename(system.index = `system:index`)

nmax <- max(str_count(tmax_2018_cbg$system.index, "_")) + 1

tmax <- tmax_2018_cbg %>%
  separate(col=system.index, into = paste0("var", seq_len(nmax)), sep="_") %>%
  mutate(year = substr(var4, 1, 4)) %>%
  mutate(date = as.Date(as.numeric(var5)-1, origin=paste0(year, "-01-01"))) %>% # it is accounting for leap years
  rename(tmax_c = mean,
         doy = var5) %>%
  select(GEOID, year, date, doy, tmax_c)

tmin_2018_cbg <- read_csv("G:/My Drive/Global_Daily_Air_Temp_CBG/tmin2018_rtp_cbg.csv") %>%
  rename(system.index = `system:index`)

nmax <- max(str_count(tmin_2018_cbg$system.index, "_")) + 1

tmin <- tmin_2018_cbg %>%
  separate(col=system.index, into = paste0("var", seq_len(nmax)), sep="_") %>%
  mutate(year = substr(var4, 1, 4)) %>%
  mutate(date = as.Date(as.numeric(var5)-1, origin=paste0(year, "-01-01"))) %>% # it is accounting for leap years
  rename(tmin_c = mean,
         doy = var5) %>%
  select(GEOID, year, date, doy, tmin_c)

rm(tmax_2018_cbg, tmin_2018_cbg)

cbg_1km <- tmax %>%
  left_join(tmin, by=c("GEOID", "year", "date", "doy")) %>%
  mutate(tmean_c_1km = (tmin_c + tmax_c)/2) %>%
  mutate(tmean_f_1km = ((tmean_c_1km*(9/5)+32)))

rm(tmax, tmin)
################################################################################
# IDENTIFY WEATHER STATION MATCHES FOR EACH CBG #######

## 1. Bring in Data ----

# Census Block Groups - get centroids
cbgs <- st_read("Data/RTP_Boundary/RTP_CBG2020.shp")

cbg_centroid <- st_centroid(cbgs)

cbgcrs <- st_crs(cbg_centroid)

# # Weather Station Data Available from Global Surface Summary of the Day for 2018
# gssod_stations <- read_excel("Data/Original/Temperature/Weather_Station/GSSOD_Stations_2018.xlsx") %>%
#   distinct(Station_ID) %>%
#   filter(!is.na(Station_ID))
#   
# # All Weather Station Data from Historical Observing Metadata Repository (https://www.ncei.noaa.gov/access/homr/reports)
# isd_history <- read_excel("Data/Original/Temperature/Weather_Station/isd_history_txt.xlsx") %>%
#   mutate(USAF = str_pad(USAF, width=6, side="left", pad="0"),
#          WBAN = str_pad(WBAN, width=5, side="left", pad="0")) %>%
#   mutate(Station_ID = paste(USAF, WBAN, sep=""))
# 
# # Create set of viable weather stations
# station_points <- isd_history %>%
#   
#   # Only those with GSSOD files for all of 2018
#   filter(Station_ID %in% unique(gssod_stations$Station_ID)) %>%
#   filter(!is.na(LON) & !is.na(LAT)) %>% # Only 1 missing in AK
#   filter(END>=20190101) %>% # Horace Williams goes offline mid-2018
#   
#   # Create Geometry 
#   st_as_sf(coords = c("LON", "LAT"), crs = cbgcrs) %>%
#   
#   # Will use row index to match back to the closest stations
#   rowid_to_column("row_index")
#   
# # Find the closest weather station to each of the CBGs based on centroid
# cbg_station_join <- st_nn(cbg_centroid, station_points, k = 1, returnDist = T)
# 
# station_info <- station_points %>%
#   st_drop_geometry() %>%
#   select(`row_index`, `STATION NAME`, `ST CALL`, `Station_ID`)
# 
# # Create Final Set of Matches
# cbg_station_match <- cbgs %>%
#   mutate(row_index = as.integer(cbg_station_join[[1]]),
#          distance_m = as.numeric(cbg_station_join[[2]])) %>%
#   left_join(station_info, by="row_index") %>%
#   st_drop_geometry()
# 
# write.csv(cbg_station_match, "Data/Original/Temperature/Weather_Station/CBG_WeatherStation_Match.csv", row.names=F)

cbg_station_match <- read_csv("Data/Original/Temperature/Weather_Station/CBG_WeatherStation_Match.csv")

# List of Files to Retrieve
stations <- unique(cbg_station_match$Station_ID)

rm(cbg_station_join, cbgcrs, gssod_stations, isd_history, station_info, station_points)

################################################################################
# CLEAN WEATHER STATION DATA #######

ws_files <- list.files(path="Data/Weather_Station/GSSOD_2018", full.names = TRUE)

ws_all <- do.call(rbind,lapply(ws_files, function(i){read_csv(i)}))

ws_clean <- ws_all %>%
  mutate(MAX = ifelse(MAX == 9999.9, NA, MAX),
         MIN = ifelse(MIN == 9999.9, NA, MIN)) %>%
  filter(TEMP_ATTRIBUTES >=12) %>% # require stations to have at least 12 hours of observations - new from 3/15 version
  # mutate(TEMP = ifelse(TEMP == 5.5, NA, TEMP)) %>%
  rename(tmean_f_ws = TEMP) %>%
  mutate(Station_ID = STATION) %>%
  rename(date = DATE) %>%
  select(Station_ID, date, tmean_f_ws)

# Find days with missing values
dates <- unique(cbg_1km$date)

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
  mutate(tmean_f_ws_nn = ifelse(flag==1, tmean_avg, tmean_f_ws)) %>%
  select(-c(tmean_avg, flag))

cbg_ws <- cbg_station_match %>%
  select(GEOID, Station_ID) %>%
  mutate(GEOID = as.numeric(GEOID)) %>%
  full_join(ws_final, by="Station_ID") 

################################################################################
# Alternative - Thin Plate Spline Regressions #######

# Weather Station Data Available from Global Surface Summary of the Day for 2018
gssod_stations <- read_excel("Data/Weather_Station/GSSOD_Stations_2018.xlsx") %>%
  distinct(Station_ID) %>%
  filter(!is.na(Station_ID))

# All Weather Station Data from Historical Observing Metadata Repository (https://www.ncei.noaa.gov/access/homr/reports)
isd_history <- read_excel("Data/Weather_Station/isd_history_txt.xlsx") %>%
  mutate(USAF = str_pad(USAF, width=6, side="left", pad="0"),
         WBAN = str_pad(WBAN, width=5, side="left", pad="0")) %>%
  mutate(Station_ID = paste(USAF, WBAN, sep=""))

# Create set of viable weather stations
station_points <- isd_history %>%
  
  # Only those with GSSOD files for all of 2018
  filter(Station_ID %in% unique(gssod_stations$Station_ID)) %>%
  filter(!is.na(LON) & !is.na(LAT)) %>% # Only 1 missing in AK
  filter((LON<=-77 & LON>=-80 & LAT>=34 & LAT<=37)) %>% # 
  filter(BEGIN<= 20180101 & END>=20190101)
 
ws_files <- list.files(path="Data/Weather_Station/TPS_Stations", full.names = TRUE)

ws_all <- do.call(rbind,lapply(ws_files, function(i){read_csv(i)}))

rm(gssod_stations, isd_history)

# Get CBG points for predictions
cbg_points <- cbg_centroid %>%
  mutate(LON = st_coordinates(geometry)[,1],
         LAT = st_coordinates(geometry)[,2]) %>%
  st_drop_geometry() %>%
  select(GEOID, LAT, LON)

geoids <- unique(cbgs$GEOID)
dates <- unique(cbg_1km$date)

# Create grid of all dates for all CBGs
cbg_all_dates <- expand.grid(dates, geoids) %>%
  rename(date = Var1,
         GEOID = Var2) %>%
  left_join(cbg_points, by=c("GEOID")) %>%
  mutate(ID = match(date, unique(date)))

date_ids <- cbg_all_dates %>%
  distinct(date, ID)

ws_clean <- ws_all %>%
  filter(TEMP_ATTRIBUTES>=12) %>% # change from 24 to 12 to be consistent with NN - new from 3/15 version
  select(STATION, DATE, TEMP, LATITUDE, LONGITUDE) %>%
  rename(Station_ID = STATION,
         date = DATE,
         LAT = LATITUDE,
         LON = LONGITUDE) %>%
  left_join(date_ids, by="date") %>%
  filter(!is.na(ID))

# List of days to loop through
Ndays <- max(cbg_all_dates$ID)

# Empty DF for Predictions
cbg_tps <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(cbg_tps) <- c("date", "GEOID", "LAT", "LON", "ID", "TEMP_pred")

# Loop through dates to predict temperature
for(day0 in 1:Ndays){
  sub.temp <- ws_clean[ws_clean$ID == day0,]
  print(day0)
  
  tick <- proc.time()
  daily_predictions <- cbg_all_dates[cbg_all_dates$ID == day0,]
  
  xLocs <- sub.temp[,c('LON','LAT')]
  
  daily_predictions$tmean_f_ws_tps <- predict(Tps(xLocs,sub.temp[,'TEMP'],lon.lat=T),x=daily_predictions[,c('LON','LAT')])
  
  cbg_tps <- rbind(cbg_tps, daily_predictions)

  tock <- proc.time()
  print(tick-tock)

}

cbg_tps <- cbg_tps %>% 
  mutate(GEOID = as.numeric(GEOID))

################################################################################

# Combine temp data from 1km grid and weather stations 
cbg_final <- cbg_1km %>% # Zhang et al. 2021, 1km gridded 
  left_join(cbg_ws, by=c("GEOID", "date")) %>% # Interpolated Weather Station, Nearest Neighbor
  left_join(cbg_tps, by=c("GEOID", "date")) %>% # Interpolated Weather Station, Thin Plate Spline Regression
  select(GEOID, year, date, tmean_f_1km, tmean_f_ws_nn, tmean_f_ws_tps)

write.csv(cbg_final, "Data/Original/Temperature/RTP_CBG_Temp_2018.csv")
