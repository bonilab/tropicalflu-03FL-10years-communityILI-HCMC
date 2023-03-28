#### ACF of ILI plus ####
rm(list = ls())
setwd("~/Dropbox/Influenza and Respiratory Disease in Vietnam/code")
load('Rdata/Overall.ILIplus.21-day.aggregated.pcr_holiday0.RData')

## ILI plus of ILI percent
annual = seq(365,1095,365)
# jpeg(filename = '../plots/ACF/ACF_annotated_iliplus_holiday_zeroed_minAD150_2020_updated.jpg',
#      width = 500,height = 300)
acf_ili_perc = acf(ili_plus_21d_df$ili_plus_smth_7d,
                   lag.max = 365*3,
                   na.action = na.pass,
                   col = "grey",ylim = c(-0.2,0.5),
                   main = '',
                   xaxt = 'n',yaxt = 'n',
                   xlab = '',ylab = '',
)

max_lag = acf_ili_perc$lag[which.max(acf_ili_perc$acf[190:365])] + 190
max_cycle = seq(max_lag,1095,max_lag)
axis(1,at = max_cycle, labels = max_cycle,cex.axis = 1.5)
axis(2, at = seq(-0.2,0.5,0.2),labels = seq(-0.2,0.5,0.2),cex.axis = 1.5)
abline(v = max_cycle,lty = "dashed")
points(x = annual, y = acf_ili_perc$acf[annual],pch = 20,cex = 1.5)
mtext("B",side = 3,outer = F,line = 1,adj = 0,cex = 2)
# dev.off()


####### peak lag in every time series ################
rm(list = ls())
load('../Rdata/20220325_ILI+/Overall.ILIplus.21-day.aggregated.pcr.RData')
peak_day_ts = data.frame()

for (i in 2012:2015){
    for (j in 2016:2019){
        
        if(i == 2012){
            start = ili_plus_21d_df$date[1]
        } else {
            start = as.Date(paste0(i,"-01-01"),tz = "UTC") 
        }
        
        end = as.Date(paste0(j,'-12-31'),tz = 'UTC')
        
        if(j - i <4) {
            
            next()
        } else {
            
            timeseries = ili_plus_21d_df$ili_plus_smth_7d[which(ili_plus_21d_df$date >= start & ili_plus_21d_df$date <= end)]
            
            acf_temp = acf(timeseries,lag.max = 365*3,
                           na.action = na.pass,
                           plot = F)
            info_this_ts = data.frame(peak_lag = which.max(acf_temp$acf[150:450]) + 150,
                                      peak_acf = max(acf_temp$acf[150:450]),
                                      end_year = j,
                                      year_number = j - i)
            peak_day_ts = rbind(peak_day_ts,info_this_ts)
            
        }
        
    }
}

peak_day_ts$year_number = as.character(peak_day_ts$year_number)
peak_day_ts

library(ggplot2)

ggplot(peak_day_ts) +
    geom_point(aes(x = end_year,y = peak_lag,color = year_number,size = peak_acf)) +
    scale_color_manual(values = c('4' = 'orange','5' = 'firebrick','6' = 'forestgreen','7' = 'navy')) +
    scale_x_continuous(breaks = unique(peak_day_ts$end_year)) +
    labs(tag = 'C') +
    xlab("End Year")+
    ylab("Peak Day")+ 
    theme_bw() +
    theme(text = element_text(size = 20),
          axis.text = element_text(size = 20),
          plot.margin = margin(t = 55,r = 0,b = 0, l = 0),
          plot.tag = element_text(size = 30),
          plot.tag.position = c(0.13,1.1)) +
    guides(color = guide_legend(title = "Number of Years\nof Data Included",
                                title.vjust = 1,
                                override.aes = list(size = 5)),
           size = guide_legend(title = "ACF Value"))

################## Periodogram ####################

### pad the time series to be exact divided by 365 ###
# jpeg(filename = '../plots/Periodogram/Periodogram_annotated_iliplus_holiday_zeroed_minAD150_2020_updated.jpg',
#      width = 500,height = 300)
ilip_spec = spectrum(ili_plus_21d_df$ili_plus_smth_7d,
                     plot = F,log = 'no',
                     pad = 0.051,
                     fast = F)

cycle_spec = data.frame(cycle = 1/ilip_spec$freq,
                        spec = ilip_spec$spec)
cycle_spec = cycle_spec[order(cycle_spec$spec,decreasing = T),]

plot(cycle_spec$cycle,cycle_spec$spec,type = 'h',
     xlim = c(0,500),
     xlab = "",ylab = "",main = "",
     lwd = 2,mgp = c(2,1,0),cex.axis = 1.5)
mtext("D",side = 3,outer = F,line = 1,adj = 0,cex = 2)
points(x = cycle_spec$cycle[1],y = cycle_spec$spec[1],
       pch = 20,cex = 1.5)
# dev.off()



