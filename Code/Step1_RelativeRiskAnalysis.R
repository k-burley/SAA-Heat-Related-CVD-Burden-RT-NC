###############################################################################
# PURPOSE: SUBGROUP ANALYSIS. This program is used to prepare data and run conditional Poisson Models 
#          that examine the association between ZCTA level temperature and cardiovascular
#          hospitalizations among Medicare cohort (65+). 
#          
# 1. Temperature exposure metrics:
#       -)  Mean Temperature
#
# 2. LAGS: 
#       -) 0 to 6 days

# 3. Meteorological control: 
#       -) ambient cubic mean RH (lag 0)
#
# 4. Time-varying control
#      -) Federal holidays (factors), summer time spline, with 4 df by day of year ns(DOY,df=4)
#       
# conditional Poisson Formula Based on Armstrong et al., 2014: http://www.ag-myresearch.com/2014_armstrong_bmcmrm.html
# gnm(formula, eliminate = NULL, ofInterest = NULL, constrain = numeric(0),
#       constrainTo = numeric(length(constrain)), family = gaussian,
#       data = NULL, subset, weights, na.action, method = "gnmFit",
#       checkLinear = TRUE, offset, start = NULL, etastart = NULL,
#       mustart = NULL, tolerance = 1e-06, iterStart = 2, iterMax = 500,
#       trace = FALSE, verbose = TRUE, model = TRUE, x = TRUE,
#       termPredictors = FALSE, ridge = 1e-08, ...)
#

# Updated on: December 26,2024 - Melissa
# Note: Some of the data used in this step cannot be publicly disclosed by the terms of the Data Use Agreement with US Centers for Medicare and Medicaid Services. 
#       These specific datasets are not included on the public GitHub and the filepaths for them were removed from the script. 

#################################################################################

####################### SECTION 1 #######################
#         LOAD LIBRARIES
#         SET WORKING DRIVE
#         Set input and output directories
#########################################################

# This section loads the needed libraries and sets the working drive where relevant datasets are stored

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
library(readxl)

# Filepaths removed for in_dir, out_dir, med_in_dir

##################################################################################
##################### SECTION 3 - Model Parameters ##############################
#
# ESTABLISH MODEL PARAMETERS
# Set User-defined functions
# Create empty dataframes and lists to hold output of interest
# Define all of the formulas to be used in the modeling loop

# ###############################################################################
zip_zcta_crosswalk<-read.csv("RTP_ZipcodeZCTA_crosswalk2009to2020.csv", sep="\t", header=T) # filepath removed
avg_num_bene_df<-fread(paste0(in_dir,"AVGNUMBENE_20092019_120CBSA_ZIP.csv"), drop = 1)
avg_num_bene_df$ZIP<-str_pad(avg_num_bene_df$ZIP,5,"0",side = "left")
#
#
rtp_temp_df<-read_csv(paste0(in_dir,"RTP_ZCTA_Temp_2010_2019.csv"))
str(rtp_temp_df)
zcta_list<-unique(rtp_temp_df$GEOID10)

rtp_temp_df<- rtp_temp_df %>%
  rename("ZCTA"="GEOID10") %>%
  mutate(ZCTA = as.character(ZCTA)) %>%                            
  dplyr::select(ZCTA,date,tmean_c) %>%
  filter(ZCTA %in% zcta_list)

zcta_merged_list<-unique(rtp_temp_df$ZCTA)

######################################################################################################################
## Join Health Data
## import ZIPCODE data, but only keep zipcodes that match ZCTAs
## Commented out lines were requred to create dataset, but it is not necessary to run them each time. The dataset that 
## Was created is saved and imported. 
######################################################################################################################

# #### Import Health Data at ZIPCODE Level
# diagno <- "CVD"
# year <- "2010"
# 
# yrs_list<-as.character(seq(2011,2019))
# med_df <- read.csv(paste0(med_in_dir,diagno,"_ZIPCODE_CBSA120_Year",year,".csv"))
# names(med_df)
# 
# med_df <- med_df %>%
#   rename("ZIP" = "BENE_ZIP_CD") %>%
#   mutate(ZIP = str_pad(ZIP,5,"0",side = "left"),
#          DATE = as.Date(paste0(substr(ADMSN_DT,1,4),"-",substr(ADMSN_DT,5,6),"-",substr(ADMSN_DT,7,8)),format="%Y-%m-%d")) %>%
#   mutate_if(is.numeric, ~replace_na(.,0)) %>%
#   dplyr::select(ZIP,DATE,CVD_TOTAL,CVD_MALE,CVD_FEMALE,CVD_UNK,CVD_WHT,CVD_BLK,CVD_OTH,CVD_NTV,CVD_85Plus) %>%
#   filter((ZIP %in% zip_zcta_list))
# 
# med_df_sum <- med_df %>%
#   group_by(ZIP) %>%
#   summarise(sum(CVD_TOTAL, na.rm = T))
# 
# ## 109 ZIPs with counts
# 
# for(i in seq(yrs_list)) {
#   year <- yrs_list[i]
#   med_yr <- read.csv(paste0(med_in_dir,diagno,"_ZIPCODE_CBSA120_Year",year,".csv"))
#   
#   med_yr <- med_yr %>%
#     rename("ZIP" = "BENE_ZIP_CD") %>%
#     mutate(ZIP = str_pad(ZIP,5,"0",side = "left"),
#            DATE = as.Date(paste0(substr(ADMSN_DT,1,4),"-",substr(ADMSN_DT,5,6),"-",substr(ADMSN_DT,7,8)),format="%Y-%m-%d")) %>%
#     mutate_if(is.numeric, ~replace_na(.,0)) %>%
#     dplyr::select(ZIP,DATE,CVD_TOTAL,CVD_MALE,CVD_FEMALE,CVD_UNK,CVD_WHT,CVD_BLK,CVD_OTH,CVD_NTV,CVD_85Plus) %>%
#     filter((ZIP %in% zip_zcta_list))
#   
#   med_df<-rbind(med_df,med_yr)
# }
# rm(med_yr)
# 
# med_df_sum <- med_df %>%
#   group_by(ZIP) %>%
#   summarise(sum(CVD_TOTAL, na.rm = T))
# 
# med_df$ZCTA<-med_df$ZIP
# noaa_df<-read.csv("NOAA_DailyWeatherData_2009-2021.csv") # filepath removed
# noaa_df<- noaa_df %>%
#   dplyr::select(ZCTA,DATE,RH_prct) %>%
#   mutate(Year = year(DATE),
#          Month = month(DATE),
#          ZCTA = as.character(ZCTA),
#          DATE = as.Date(DATE)) %>%
#   filter(Year<2020 &
#            Month<=10 & Month >=4) %>%
#   filter(ZCTA %in% zip_zcta_list)
# 
# 
# 
# rtp_df<-left_join(rtp_temp_df,med_df, by = c("ZCTA","date"="DATE")) %>%
#   mutate_if(is.numeric, ~replace_na(.,0)) %>%
#   mutate(Year = year(date)) %>%
#   left_join(noaa_df, by = c("ZCTA","date"="DATE", "Year")) %>%
#   mutate(ZIP = ifelse(is.na(ZIP), ZCTA, ZIP)) %>%
#   left_join(avg_num_bene_df, by = c("ZIP", "Year"="YEAR"))
# 
# 
# length(unique(rtp_df$ZCTA))   ## 84 ZCTA
# 
# rm(med_df, rtp_temp_df,noaa_df)
# write_feather(rtp_df,paste0(in_dir,"RTP_ZIP_ZCTA_CVD_EXPO_20102019.feather"))

#-------------------------------------------------------------------------------------------------------------------
rtp_df<-read_feather(paste0(in_dir,"RTP_ZIP_ZCTA_CVD_EXPO_20102019.feather"))
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
# ------- Formula ------------
#-------------------------------------------------------------------------------------

formula.i <- "cb_rtp.ns2knots.lag6 + ns(RH_prct,df=4) + holiday_f + ns(DOY,4)+offset(log(TotalBENE_TOTAL))"

#-------------------------------------------------------------------------------------------------------------------

# Define function to calculate Q-AIC - from: https://github.com/gasparrini/2013_gasparrini_BMCmrm_Rcodedata
fqaic <- function(model) {
  loglik <- sum(dpois(model$y,model$fitted.values,log=TRUE))
  phi <- summary(model)$dispersion
  qaic <- -2*loglik + 2*summary(model)$df[3]*phi
  return(qaic)
}

# INITIALISE DATA FRAME TO HOLD MODEL FIT STATISTICS TO FACILATE COMPARISON AND DEFINE FUNCTION TO EXTRACT THEM
setNames(modelfitstats <- data.frame(matrix(NA,0,7)),c("model","df.residual", "dfmodel","dispersion","deviance","AIC","QAIC"))

#function to extract model fit statistics
extractmodelfit <- function(fit,name) {
  modelfit <- data.frame(name, 
                         summary(fit)$df.residual,
                         summary(fit)$df[1],
                         summary(fit)$dispersion,       
                         summary(fit)$deviance,
                         summary(fit)$deviance + 2*summary(fit)$df[1]*summary(fit)$dispersion,
                         fqaic(fit))                                  #my_BIC  <- dev +  p * log(n)
  names(modelfit) <- c("model","df.residual","dfmodel","dispersion","deviance","AIC","QAIC")
  return(modelfit)
}

all_RRs_bylag<-NULL
all_cumulRRs_bylag<-NULL


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

##################################################################################
##################### SECTION 4 - RUN THE MODELS ##############################
# Run 4 loop models that examine mean temperature - resp associations in each CBSA for each resp outcome (n=17)
# Total # of models 17*120 = 2280
###############################################################################

########################################################################################################################################################
## Warm Season
## outcome = "CVD_TOTAL"
########################################################################################################################################################

cbsa_name <- "RTP, NC"
#---------------------------------------------------------------------------------------------------------------

subgroup_name<-"CVD_RTP_2575.nsarglag_ZIPZCTA_6DayLag_YD"

######################################################################################################################################################
## Create Data and Variables for model grouping and strata
######################################################################################################################################################

rtp_df<-as.data.frame(rtp_df)


rtp_df<-rtp_df[order(rtp_df$ZCTA,rtp_df$date), ]

# rtp_df<- rtp_df %>%
#   group_by(ZCTA) %>%
#   mutate(Avg.Temp = mean(tmean_c, na.rm = T)) %>%
#   ungroup() 

rtp_df$year_f   <- as.factor(format(rtp_df$date, format="%Y"))

rtp_df$DOY <- as.numeric(strftime(rtp_df$date,format="%j"))

rtp_df$month_f  <- as.factor(months(rtp_df$date))

rtp_df$dow_f    <- as.factor(weekdays(rtp_df$date))
print(table(rtp_df$dow_f))

rtp_df$zcta_f    <-as.factor(rtp_df$ZCTA)

rtp_df <- rtp_df %>%
  mutate(DOW_NUM=wday(rtp_df$date),
         incre_mult=rep(c(0,7,14), len = nrow(rtp_df), each=7)) %>%
  mutate(dow_num_f=as.factor(DOW_NUM+incre_mult))

#need to have stratum variable for the modeling
rtp_df$stratum  <- as.factor(rtp_df$year_f:rtp_df$dow_num_f)

# Crossbasis grouping variable
rtp_df$zcta_yr_f  <- as.factor(paste0(rtp_df$year_f,":",rtp_df$zcta_f))


# Federal holidays: NewYears, Independence, Veterans, Christmas, MLK, Presidents, Memorial, Labor, Columbus, and Thanksgiving
hld <- tis::federalHolidays(2010:2019, board = T, businessOnly = T) 
holid <- as.Date(as.character(hld), format = "%Y%m%d")
rtp_df$holiday <- ifelse(rtp_df$date %in% holid, 1, 0)
rtp_df$holiday_f<-as.factor(rtp_df$holiday)

CBSA.Title<-"RTP, NC"

#Grab the name of the exposure used
exposure<-"tmean_c"

print(exposure)

#Grab the name of the formula used
formula_name<-substr(formula.i,4,21)
print(formula_name)

tspline="DOY4"


basis<-substr(formula.i,1,21) #grab name of the crossbasis used in the formula. all crossbasis names have the same length

print(basis)
#Get name for saving output and identifying which model was run  

name=paste0(subgroup_name)
print(name)


#Get percentile values for centering and prediction
#cen_50 is the 50th percentile and the a priori centering value 
cen_50<-round(quantile(eval(parse(text=paste0("rtp_df$",exposure))),0.5,na.rm=T),2) 
print(cen_50)

#Generate cumulative RRs comparing 95th to 50th percentile,
#95th percentile value
at95=round(quantile(eval(parse(text=paste0("rtp_df$",exposure))),0.95,na.rm=T),2)
at99=round(quantile(eval(parse(text=paste0("rtp_df$",exposure))),0.99,na.rm=T),2)

#################################################################################################################################   
####### CROSSBASIS FOR MEAN TEMPERATURE #########################

cb_rtp.ns2knots.lag6 <- crossbasis(rtp_df$tmean_c, lag = c(0,6), #constrain the lag effects to account for correlation across lag days
                                   argvar =  list(fun = "ns", knots=quantile(rtp_df$tmean_c,c(0.25,0.75),na.rm=T)),
                                   arglag = list(fun="ns",knots=logknots(c(0,6),df=5)), group = rtp_df$zcta_yr_f)


#stratum based on year:dow

# #EXCLUDE EMPTY STRATA, OTHERWISE BIAS IN gnm WITH quasipoisson - This is commonly done, but often not explicitly stated in the methods.
# For the best example, see the supplement of https://www.ahajournals.org/doi/10.1161/CIRCULATIONAHA.122.061832
# https://github.com/gasparrini/CTS-smallarea/blob/main/03.mainmod.R
#	IMPORTANT FOR CONDITIONAL POISSON. IF YOU ACCIDENTLY INCLUDE THE EMPTY STRATA, THESE STRATA STILL WON'T BE INFORMATIVE TO THE ESTIMATE, BUT YOUR N WILL BE ARTIFICIALLY TOO BIG, SO YOUR CONFIDENCE INTERVALS WILL BE SMALLER THAN THEY SHOULD BE

#keeping empty strata make standard errors smaller than they should be
#to remove empty strata do subset=keep in the formula


rtp_df<-data.table(rtp_df)
rtp_df[,keep:=sum(CVD_TOTAL)>0, by=stratum]


# y is now your dependent variable (i.e., outcome, or the number of events)
y = "CVD_TOTAL"
z= c(formula.i) #z is now the the specific formula string in the loop
model = reformulate(z, response = y) #this writes out the model formula with the correct specification
print(model)

## Warm Season
fit <-try(gnm(model,data =rtp_df, family = "quasipoisson",
              eliminate = factor(stratum), na.action = "na.exclude",subset=keep))

print(class(fit))

#Saves the model coeficients
coef <- coef(fit)

#Saves the variance-covariance matrix
vcov <- vcov(fit)

#This creates a dataframe to add all of the model fit statistics for each run of the loop
modelfitstats <- rbind(modelfitstats,extractmodelfit(fit,name))
#
# #Handy code to capture model summary for record keeping purposes
o<-capture.output(summary(fit))

writeLines(o, paste0(out_dir,name,"_summary_fit.txt"))
# Over dispersion can be detected by dividing the residual deviance by the degrees of freedom (residual?).
# If this quotient is much greater than one, the model is overdispersed and quasipoisson is required


################### defining variables for prediction and plotting ###############
# cen_50 is the 50th percentile of the CBSA-specific curve (defined above, early in the loop)
# at95 is the 95th percentile of the CBSA-specific curve (defined above, early in the loop)

lagRRat95<-as.data.frame(crosspred(cb_rtp.ns2knots.lag6, fit,cen=cen_50, at=at95)[c("matRRfit","matRRlow","matRRhigh")]) #RRs for each lag comparing 50th to 95th percentile
cumlagRRat95<-as.data.frame(crosspred(cb_rtp.ns2knots.lag6, fit,cen=cen_50, at=at95,cumul=T)[c("cumRRfit","cumRRlow","cumRRhigh")])
#-------------------------------------------------------------------------------------------------------------------------------------------------------------

overall.plot<- crosspred(cb_rtp.ns2knots.lag6, fit,cen=cen_50)
lagslices.plot<-crosspred(cb_rtp.ns2knots.lag6, fit, at=at95,cen=cen_50) #plot of RR at each lag
lagcumul.plot<-crosspred(cb_rtp.ns2knots.lag6, fit, cen=cen_50, cumul=TRUE, bylag=1, at=at95)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------
lagRRat95$ID<-1

long <- melt(setDT(lagRRat95), id.vars = 22)
long$lagday<-rep(c("lag0","lag1","lag2","lag3","lag4","lag5","lag6"),3)
long$RRtype<-rep(c("RRfit","RRlow","RRhigh"),each = 7, len = 21)

long$RRtype<-as.factor(long$RRtype)
lagRR.95<-dcast(long, lagday ~ RRtype)
lagRR.95$LAGnum<-as.numeric(gsub("lag","",lagRR.95$lagday))
lagRR.95$CBSA.ID<-cbsa_name
lagRR.95$CBSA.Name<-CBSA.Title
lagRR.95$cen_per<-"50th percentile"
lagRR.95$contrast_per<-"95th percentile"
lagRR.95$model<-name
lagRR.95<-lagRR.95[order(lagRR.95$LAGnum),]
rownames(lagRR.95)<-lagRR.95$LAGnum
lagRR.95$LAGnum<-NULL

#-------------------------------------------------------------------------------------------------------------------------------------------------------------
# #Generate cumulative RRs at each lag comparing 95th (flipped percentiles) to 88th percentile, create a dataframe that holds the values and adds rows for each element in the loop
cumlagRRat95$ID<-1
# 
long <- melt(setDT(cumlagRRat95), id.vars = 22)
long$lagday<-rep(c("lag0","lag1","lag2","lag3","lag4","lag5","lag6"),3)
long$RRtype<-rep(c("cumulRRfit","cumulRRlow","cumulRRhigh"),each = 7, len = 21)

long$RRtype<-as.factor(long$RRtype)
cumullagRR.95<-dcast(long, lagday ~ RRtype)
cumullagRR.95$CBSA.ID<-cbsa_name
cumullagRR.95$CBSA.Name<-CBSA.Title
cumullagRR.95$cen_per<-"50th percentile"
cumullagRR.95$contrast_per<-"95th percentile"
cumullagRR.95$model<-name
cumullagRR.95$LAGnum<-as.numeric(gsub("lag","",cumullagRR.95$lagday))
cumullagRR.95<-cumullagRR.95[order(cumullagRR.95$LAGnum),]
rownames(cumullagRR.95)<-cumullagRR.95$LAGnum
cumullagRR.95$LAGnum<-NULL


#combine all RRs for each lag day together and save
all_RRs_bylag<-rbind(all_RRs_bylag,lagRR.95)

#combine all cumulative RRs for each lag day together and save
all_cumulRRs_bylag<-rbind(all_cumulRRs_bylag,cumullagRR.95) #sneaky way to concatenate new rows created by each loop


#### PERFORM THE CROSSREDUCE for each CBSA ####x#####
#https://github.com/gasparrini/2015_gasparrini_EHP_Rcodedata

# CROSSREDUCE
#Reduces model parameters to overall cumulative function. The crossreduce here will have 3 variables representing temperature splines
#Reduced variables are centered at each CBSAs' 50th percentile.
#crossreduce info: https://www.rdocumentation.org/packages/dlnm/versions/2.4.7/topics/crossreduce

reduced<- crossreduce(cb_rtp.ns2knots.lag6,type="overall",fit,
                      cen=cen_50,model.link="log")

overall.coef<-coef(reduced) #store the coef 
overall.vcov<-vcov(reduced) # store the vcov 
############# PERFORM CROSSREDUCE TO GET LAG EFFECT - OVERALL REPORTING
#NOTE: This is tricky. Centering value and prediction value is super important to specify here. You must use the centering and prediction values you want to report
#type="var" gives the lag-response association for a particular contrast value (in this case, value is 5th percentile) for all lags used in prediction
# Below, we center at the cbsa specific 50th percentile (cen=cen_50), and we compare to the cbas-specific 95th percentile ("value=at5)
# This cross-reduce extracts lag-day specific associations and not cumulative lag-day associations

lag_red<- crossreduce(cb_rtp.ns2knots.lag6,type="var",
                      value=at95,
                      fit,cen=cen_50,model.link="log")

lag.coef_red<-coef(lag_red)
lag.vcov_red<-vcov(lag_red)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plots
#file name for saving the plots
graph.file.out <- paste0(out_dir, name,".pdf")
print(graph.file.out)

pdf(file = paste0(graph.file.out),   # The directory you want to save the file in
    width = 9, # The width of the plot in inches
    height = 12) # The height of the plot in inches
layout(matrix(c(1,2,3,4,5,6),3,2,byrow=T),respect=T)
par(mfrow=c(3,2))

par(mar=c(4,4.5,10,2.4))

plot(overall.plot, "overall", xlab="Temperature", ylab="RR" ,
     main=paste0("Cumulative RR by mean temperature centered at 50th percentile \n in RTP"))
abline(v = at95, col = "red", lty = 2)
abline(v = cen_50, col = "black", lty = 2)

####### PLOT THE ESTIMATED LAG STRUCTURE for each lag comparing the 95th percentile to the Median
## Histogram of temperatures
hist(rtp_df$tmean_c, breaks = 20)

plot(lagslices.plot, "slices", var=at95, ci="bars", type="p", col=2, pch=19,
     xlab="Lag Day", ylab="RR", ylim=c(.90,1.1),main=paste0("Lag-response at 95th percentile \ncompared to 50th percentile in ",CBSA.Title))


#this plots the cumulative risk across lag days
plot(lagcumul.plot, "slices", var=at95, col=2, cumul=TRUE, ylab="Cumulative RR, 95th percentile",
     ylim=c(.9,1.1),main=paste0("Cumulative Lag-response at 95th percentile \ncompared to 50th percentile in ",CBSA.Title))


dev.off()
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------     

gc()


model_info_out=paste0(out_dir,"allmodeloutput_",subgroup_name,".RData")
#Saves the coef and vcov info into an R list of lists. Should be able to use this for 2nd stage meta analysis
save(coef,vcov,file=model_info_out)

#Save the model fit statistics
save(modelfitstats,file=paste0(out_dir,subgroup_name,"modelfitstats.RData")) #save in O drive

RR_lag_file_out <- paste0(out_dir, subgroup_name,"all_RRs_bylag.RDS")
saveRDS(all_RRs_bylag,file=RR_lag_file_out)

#Save all cumulative RRs for each lag day
cumRR_bylag.file_out <- paste0(out_dir,subgroup_name,"all_cumulRRs_bylag.RDS")
saveRDS(all_cumulRRs_bylag,file=cumRR_bylag.file_out)    

#file name for reduced coef and reduced vcov
red_info_out=paste0(out_dir,"reduced_coef_vcov_",subgroup_name,".RData")
#Saves the coef and vcoc from reduced models in an r list of lists. 
save(overall.coef,overall.vcov,file=red_info_out)

#file name for reduced coef and reduced vcov for LAGGED EFFECT
lag_red_info_out=paste0(out_dir,"LAG_reduced_coef_vcov_",subgroup_name,".RData")
#Saves the coef and vcoc from reduced models in an r list of lists. 
#saves lag-response 95th to 50th percentile
save(lag.coef_red,lag.vcov_red,file=lag_red_info_out)






