################################################################################
# Program Name: Format_Heat_Data.R
# Program Purpose: Format air temperature data for ZCTA

# Heat Data processed in Google Earth Engine: https://code.earthengine.google.com/6f5fff216b6ed42c7f2a27126f00a410
# and ArcGIS for 3 days that are missing from GEE, https://iastate.figshare.com/collections/A_global_seamless_1_km_resolution_daily_land_surface_temperature_dataset_2003_2020_/5078492 
# Export by year and metric (mean of tmax, mean of tmin)

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
}

## 1. Set Up ----

# Files too large for GitHub, available on Google Drive: https://drive.google.com/drive/folders/1RFte5bKwYLHzHVt2So5bZWuHbAaW4sri?usp=drive_link

# 2. Format + Combine ZCTA Level Data from Google Earth Engine

"G:\My Drive\Global_Daily_Air_Temp\tmin2019_rtp_zcta.csv"

# TMAX
tmax_files <- list.files(path="G:/My Drive/Global_Daily_Air_Temp", pattern = "tmax", full.names = TRUE)

tmax_all_orig <- do.call(rbind,lapply(tmax_files, function(i){read.csv(i, encoding="UTF-8")})) #tmax obs for all years

nmax <- max(str_count(tmax_all_orig$system.index, "_")) + 1

tmax_all <- tmax_all_orig %>%
  separate(col=system.index, into = paste0("var", seq_len(nmax)), sep="_") %>%
  mutate(year = substr(var4, 1, 4)) %>%
  mutate(date = as.Date(as.numeric(var5)-1, origin=paste0(year, "-01-01"))) %>% # it is accounting for leap years
  rename(tmax_c = mean,
         doy = var5) %>%
  select(GEOID10, year, date, doy, tmax_c)

rm(tmax_all_orig)

# TMIN
tmin_files <- list.files(path="G:/My Drive/Global_Daily_Air_Temp", pattern = "tmin", full.names = TRUE)

tmin_all_orig <- do.call(rbind,lapply(tmin_files, function(i){read.csv(i, encoding="UTF-8")})) #tmax obs for all years

tmin_all <- tmin_all_orig %>%
  separate(col=system.index, into = paste0("var", seq_len(nmax)), sep="_") %>%
  mutate(year = substr(var4, 1, 4)) %>%
  mutate(date = as.Date(as.numeric(var5)-1, origin=paste0(year, "-01-01"))) %>% # it is accounting for leap years
  rename(tmin_c = mean,
         doy = var5) %>%
  select(GEOID10, year, date, doy, tmin_c)

rm(tmin_all_orig)

# Bring in the separately processed missing data
# ARCGIS Processing:
# 1. Bring in ZCTA data (RTP_CBSAZCTA with 86 ZCTA)
# 2. Bring in rasters individually from the Z drive, downloaded from https://iastate.figshare.com/collections/A_global_seamless_1_km_resolution_daily_land_surface_temperature_dataset_2003_2020_/5078492 
# 3. Calculate "zonal statistics as table" mean value using ZCTA5CE10
# 4. Inspect attribute table and export as CSV

# NOTE: In 2010 ZCTAs, one ZCTA (27568) is too small for ArcGIS zonal stats... manually replaced farther down
 
tmin_2012_135 <- read_csv("G:/My Drive/Global_Daily_Air_Temp/Additional_Days/zonal_statistics_tmin_2012_135.csv") %>%
  mutate(year = "2012",
         doy = "135") %>%
  mutate(tmin_c_new = MEAN*0.1) %>% # orig is in 0.1 celsius
  select(c("ZCTA5", "year", "doy", "tmin_c_new")) %>%
  rename(GEOID10 = ZCTA5)

tmax_2012_135 <- read_csv("G:/My Drive/Global_Daily_Air_Temp/Additional_Days/zonal_statistics_tmax_2012_135.csv") %>%
  mutate(year = "2012",
         doy = "135") %>%
  mutate(tmax_c_new = MEAN*0.1) %>% # orig is in 0.1 celsius
  select(c("ZCTA5", "year", "doy", "tmax_c_new")) %>%
rename(GEOID10 = ZCTA5)

t_2012_135 <- tmin_2012_135 %>%
  left_join(tmax_2012_135, by=c("GEOID10", "year", "doy"))

tmin_2014_257 <- read_csv("G:/My Drive/Global_Daily_Air_Temp/Additional_Days/zonal_statistics_tmin_2014_257.csv") %>%
  mutate(year = "2014",
         doy = "257") %>%
  mutate(tmin_c_new = MEAN*0.1) %>% # orig is in 0.1 celsius
  select(c("ZCTA5", "year", "doy", "tmin_c_new")) %>%
  rename(GEOID10 = ZCTA5)

tmax_2014_257 <- read_csv("G:/My Drive/Global_Daily_Air_Temp/Additional_Days/zonal_statistics_tmax_2014_257.csv") %>%
  mutate(year = "2014",
         doy = "257") %>%
  mutate(tmax_c_new = MEAN*0.1) %>% # orig is in 0.1 celsius
  select(c("ZCTA5", "year", "doy", "tmax_c_new")) %>%
  rename(GEOID10 = ZCTA5)

t_2014_257 <- tmin_2014_257 %>%
  left_join(tmax_2014_257, by=c("GEOID10", "year", "doy"))

tmin_2019_146 <- read_csv("G:/My Drive/Global_Daily_Air_Temp/Additional_Days/zonal_statistics_tmin_2019_146.csv") %>%
  mutate(year = "2019",
         doy = "146") %>%
  mutate(tmin_c_new = MEAN*0.1) %>% # orig is in 0.1 celsius
  select(c("ZCTA5", "year", "doy", "tmin_c_new")) %>%
  rename(GEOID10 = ZCTA5)

tmax_2019_146 <- read_csv("G:/My Drive/Global_Daily_Air_Temp/Additional_Days/zonal_statistics_tmax_2019_146.csv") %>%
  mutate(year = "2019",
         doy = "146") %>%
  mutate(tmax_c_new = MEAN*0.1) %>% # orig is in 0.1 celsius
  select(c("ZCTA5", "year", "doy", "tmax_c_new")) %>%
  rename(GEOID10 = ZCTA5)

t_2019_146 <- tmin_2019_146 %>%
  left_join(tmax_2019_146, by=c("GEOID10", "year", "doy"))

rm(tmin_2012_135, tmax_2012_135, tmin_2014_257, tmax_2014_257, tmin_2019_146, tmax_2019_146)

missing <- t_2012_135 %>%
  bind_rows(t_2014_257) %>%
  bind_rows(t_2019_146)

rm(t_2012_135, t_2014_257, t_2019_146)

# COMBINE
zcta_final <- tmax_all %>%
  left_join(tmin_all, by=c("GEOID10", "year", "date", "doy")) %>%
  # Fill in missing data
  left_join(missing, by=c("GEOID10", "year", "doy"))

# CHECK: How closely do ArcGIS and GEE estimates align?
check <- zcta_final[is.na(zcta_final$tmin_c),]
cor.test(check$tmax_c, check$tmax_c_new, method = "pearson") # 0.9999863 
ggplot(data=check, aes(x=tmax_c, y=tmax_c_new)) + geom_point()

check2 <- zcta_final[is.na(zcta_final$tmax_c),]
cor.test(check2$tmin_c, check2$tmin_c_new, method = "pearson") # 0.9999414 
ggplot(data=check2, aes(x=tmin_c, y=tmin_c_new)) + geom_point()

# REPLACE MISSING
zcta_final <- zcta_final %>%
  mutate(tmin_c = ifelse(is.na(tmin_c), tmin_c_new, tmin_c),
         tmax_c = ifelse(is.na(tmax_c), tmax_c_new, tmax_c)) %>%
  # Manual check + fixes on one tiny ZCTA that is too small for zonal statistics in ArcGIS
  mutate(tmin_c = ifelse(GEOID10 == 27568 & year == 2012 & doy == 135, 15.4, tmin_c),
         tmax_c = ifelse(GEOID10 == 27568 & year == 2014 & doy == 257, 25.8, tmax_c),
         tmin_c = ifelse(GEOID10 == 27568 & year == 2019 & doy == 146, 20.8, tmin_c)) %>%
  mutate(tmean_c = (tmin_c + tmax_c)/2) %>%
  select(-c(tmin_c_new, tmax_c_new))

table(is.na(zcta_final$tmax_c))
table(is.na(zcta_final$tmin_c))
table(is.na(zcta_final$tmean_c))

ggplot(zcta_final, aes(x=tmax_c)) + 
  geom_histogram(color="black", fill="white")

ggplot(zcta_final, aes(x=tmin_c)) + 
  geom_histogram(color="black", fill="white")

ggplot(zcta_final, aes(x=tmean_c)) + 
  geom_histogram(color="black", fill="white")

write.csv(zcta_final, "Data/Analysis/RTP_ZCTA_Temp_2010_2019.csv", row.names=F)




