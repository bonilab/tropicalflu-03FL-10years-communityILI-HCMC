### Calculate ILI+ from aggregated pcr fraction and daily zeta score ###
# Note: pcr reporting is incomplete after 2019-06-26, use 20201109 version to calculate ili+

setwd('~/Dropbox/Influenza and Respiratory Disease in Vietnam/code')
rm(list = ls())
library(dplyr)
load('Rdata/PCR data aggregated by 21d_holiday_zeroed_20220210.RData')
load('Rdata/7d.smth.zeta.ilipercent.Rdata')
source('functions/Detrend and Smooth.R')
source('functions/iliplus data analysis.R')

#### Overall ILI+ from daily zeta and 21-day aggregated PCR ####
#### ILI+ = ZETA * positive percent in the 21 day-interval ####

#### Use the latest file ###
ili_pcr = read.csv("datasets/AllPCRandSubtypeResults.03FLand17FL.20220210.csv",
                   header = T,stringsAsFactors = F,
                   sep = ',')

ili_pcr$SampleDate = paste0(ili_pcr$SampleYear,"-",ili_pcr$SampleMonth,"-",ili_pcr$SampleDay)
ili_pcr$SampleDate = as.Date(ili_pcr$SampleDate,tz = "UTC")

timeperiod = seq(min(ili_pcr$SampleDate),max(ili_pcr$SampleDate),'1 day')
timeperiod = zeta_perc_avg_df$Date[zeta_perc_avg_df$Date %in% timeperiod]

#### Number of nasal swab from May 23 2012 to Dec 31 2019 ###
nrow(ili_pcr[ili_pcr$SampleDate %in% timeperiod,])

all_flu = ili_pcr[ili_pcr$FluA == 1 | ili_pcr$FluB == 1,]

## Number of samples positive for influenza ##
## 545(20.9%) patients are positive for flu 
## 6 (0.2%) patients are positive for flu A and flu B
## 163 (6.3%) for H1N1 
## 169 (6.5%) for H3N2
## 209 (8%) for B
ili_pcr %>%
    filter(SampleDate %in% timeperiod) %>%
    summarise(flu = sum(FluA == 1 | FluB == 1)/2604,
              both = sum(FluA == 1 & FluB == 1),
              H1 = sum(H1 == 1),
              H3 = sum(H3 == 1),
              B = sum(FluB == 1)/2604,
              A = sum(FluA == 1))

aggregated_pcr_21d = aggregated_pcr_21d[aggregated_pcr_21d$firstday %in% timeperiod,]

##### Subtype ILI+ from daily zeta and 21-day subtype pos pcr rate #####
all_daily_subtype_ilip = data.frame()
ili_plus_21d_df = data.frame()

for(i in 1:nrow(aggregated_pcr_21d)){
    
    if(i == nrow(aggregated_pcr_21d)) {
        
        period_this_interval = seq(aggregated_pcr_21d$firstday[i],max(timeperiod),"1 day")
        
    } else {
        
        period_this_interval = seq(aggregated_pcr_21d$firstday[i],(aggregated_pcr_21d$firstday[i+1] -1),"1 day")
    }
    
    zeta_this_interval = zeta_perc_avg_df$zeta_score_avg[which(zeta_perc_avg_df$Date %in% period_this_interval)]
    
    daily_ilip = aggregated_pcr_21d$pos_perc_np_swab[i]* zeta_this_interval/100
    daily_H1_ilip = round(aggregated_pcr_21d$pos_perc_H1_np_swab[i] * zeta_this_interval/100,2)
    daily_H3_ilip = round(aggregated_pcr_21d$pos_perc_H3_np_swab[i] * zeta_this_interval/100,2)
    daily_B_ilip = round(aggregated_pcr_21d$pos_perc_B_np_swab[i] * zeta_this_interval/100,2)
    
    all_daily_subtype_ilip = rbind(all_daily_subtype_ilip,
                                   data.frame(date = period_this_interval,
                                              daily_H1_ilip = daily_H1_ilip,
                                              daily_H3_ilip = daily_H3_ilip,
                                              daily_B_ilip = daily_B_ilip))
    ili_plus_21d_df = rbind(ili_plus_21d_df,
                          data.frame(date = period_this_interval,
                                     daily_ilip = daily_ilip))
}

ili_plus_21d_df$ili_plus_smth_7d = moving_average(ili_plus_21d_df$daily_ilip,window_width = 7)

#### NA at 2019-12-31 in iliplus is because zeta is NA on that day
all_subtype_ilip = list(H1_ilip = data.frame(date = all_daily_subtype_ilip$date,
                                             daily_ilip = all_daily_subtype_ilip$daily_H1_ilip,
                                             ilip_smth_7d = round(moving_average(all_daily_subtype_ilip$daily_H1_ilip,7),2)),
                        H3_ilip = data.frame(date = all_daily_subtype_ilip$date,
                                             daily_ilip = all_daily_subtype_ilip$daily_H3_ilip,
                                             ilip_smth_7d = round(moving_average(all_daily_subtype_ilip$daily_H3_ilip,7),2)),
                        B_ilip = data.frame(date = all_daily_subtype_ilip$date,
                                            daily_ilip = all_daily_subtype_ilip$daily_B_ilip,
                                            ilip_smth_7d = round(moving_average(all_daily_subtype_ilip$daily_B_ilip,7),2)))

save(ili_plus_21d_df,file = 'Rdata/Overall.ILIplus.21-day.aggregated.pcr.Rdata')
save(all_subtype_ilip,file = "Rdata/Subtype.iliplus.21d.aggregated.pcr.Rdata")

#### Plot stacked iliplus ####
all_subtype_smth_df = data.frame(date = all_subtype_ilip$H1_ilip$date,
                                 H1plus = all_subtype_ilip$H1_ilip$ilip_smth_7d,
                                 H3plus = all_subtype_ilip$H3_ilip$ilip_smth_7d,
                                 Bplus = all_subtype_ilip$B_ilip$ilip_smth_7d)

plot_stacked_iliplus_by_year(iliplus_data = all_subtype_smth_df,
                             Bplus_label = "Influenza B",
                             H1plus_label = "Influenza H1N1",
                             H3plus_label = "Influenza H3N2",
                             tag = 'A',
                             change_tag_position = T,
                             margins = c(t = 40,b = 0,l = 0, r = 0),
                             tag_coordiantes = c(0.19,1.03))
                             


