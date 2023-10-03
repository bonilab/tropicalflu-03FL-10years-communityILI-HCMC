rm(list = ls())

source("functions/Detrend and Smooth.R")
source('functions/plotting.R')
library(ggplot2)
library(dplyr)
library(reshape2)

########## Zeta score of ILI percent ##########

### load ILIPERCENT ###
load("Rdata/clean ILIPERCENT.RData")

### Get the zeta score for each clinic
ZETA_PERC = apply(PERCENT,2,zeta_score,window_width = 365,min_available_days = 150)

colnames(ZETA_PERC) = colnames(PERCENT)

##### Average zeta scores across all the clinics
zeta_perc_avg = rowMeans(ZETA_PERC,na.rm = T)

rm(PERCENT,ZETA_PERC)

##### 7-day averaged zeta score of ILI percent and stack by years
zeta_perc_avg_smth7 = moving_average(zeta_perc_avg,window_width = 7)
zeta_perc_avg_df = data.frame(Date = seq(as.Date("2010-01-01",tz = "UTC"),as.Date("2019-12-31",tz = "UTC"),by = "1 day"),
                              zeta_score_avg = zeta_perc_avg,
                              smoothed_zeta_score = zeta_perc_avg_smth7)

#### base plot version ####

#### Entire time series ###
plot(x = zeta_perc_avg_df$Date,y = zeta_perc_avg_df$zeta_score_avg,main = "Zeta Score During 2010-01-01 to 2019-12-31",type = "l",xlab = "Date",ylab = "",
     xaxs = "i")
annual = seq(as.Date("2010-01-01"),as.Date("2020-01-01"),"1 year")
abline(v = annual, lty = "dashed",col = "grey")
lines(x = zeta_perc_avg_df$Date,y = zeta_perc_avg_df$smoothed_zeta_score,col = "red",lwd = 1.5)
legend("topright",legend = c("Daily Zeta Score","7 day-window Smoothed Zeta Score"),col = c("black","red"),
       lty = c("solid","solid"),lwd = c(1,1.5))

##### zeta score by year ###
plot_data_by_year(data_df = zeta_perc_avg_df,data_name = 'zeta_score_avg',
                  time_name = 'Date',ymin = 0.5,ymax = 1.6,
                  smoothed_data_name = 'smoothed_zeta_score',Nrow = 5,Ncol = 2,
                  tag = '',
                  save = F,
                  width = 12,
                  height = 10)




save(zeta_perc_avg_df,file = "Rdata/7d.smth.zeta.ilipercent.Rdata")

rm(list = ls())
