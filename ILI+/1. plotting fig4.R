library(ggpubr)
library(cowplot)
library(gridGraphics)
library(ggplot2)
library(RColorBrewer)

rm(list = ls())
setwd("~/Dropbox/Influenza and Respiratory Disease in Vietnam/code")
source('functions/iliplus data analysis.R')
load('Rdata/Overall.ILIplus.21-day.aggregated.pcr_holiday0.RData')
load('Rdata/Subtype.iliplus.21d.aggregated.pcr.Rdata')

all_subtype_smth_df = data.frame(date = all_subtype_ilip$H1_ilip$date,
                                 H1plus = all_subtype_ilip$H1_ilip$ilip_smth_7d,
                                 H3plus = all_subtype_ilip$H3_ilip$ilip_smth_7d,
                                 Bplus = all_subtype_ilip$B_ilip$ilip_smth_7d)

ilip_p = plot_stacked_iliplus_by_year(iliplus_data = all_subtype_smth_df,
                             Bplus_label = "Influenza B",
                             H1plus_label = "Influenza H1N1",
                             H3plus_label = "Influenza H3N2",
                             tag = '',
                             change_tag_position = T,
                             margins = c(t = 40,b = 0,l = 0, r = 0),
                             tag_coordiantes = c(0.19,1.03))
ggsave(filename = 'plots/fig4a.jpeg',plot = ilip_p,
       width = 6,height = 10,units = 'in')


#### ACF of ILI plus ####

## ILI plus of ILI percent
annual = seq(365,1095,365)

acf_ili_plus = acf(ili_plus_21d_df$ili_plus_smth_7d,
                   lag.max = 365*3,
                   na.action = na.pass,
                   col = "grey",ylim = c(-0.2,0.5),
                   main = '',
                   xaxt = 'n',yaxt = 'n',
                   xlab = '',ylab = '',
                   plot = F)

max_lag = acf_ili_plus$lag[which.max(acf_ili_plus$acf[190:365])] + 190
max_cycle = seq(max_lag,1095,max_lag)

acf_df = data.frame(lag = acf_ili_plus$lag,
                    acf = acf_ili_plus$acf)
acf_p = 
    ggplot(acf_df) + 
    geom_col(aes(x = lag, y = acf),fill = 'grey') + 
    geom_point(aes(x = annual, y = acf),
               data = data.frame(annual = seq(365,1095,365), 
                                 acf = acf_df$acf[which(acf_df$lag %in% seq(365,1095,365))]),
               size = 2) + 
    geom_vline(xintercept = max_cycle,linetype = 'dashed') +
    geom_hline(yintercept = c(-1.96/sqrt(nrow(acf_df)),
                              1.96/sqrt(nrow(acf_df))),
               linetype = 'dashed',
               color = 'blue') + 
    scale_x_continuous(breaks = max_cycle) + 
    scale_y_continuous(limits = c(-0.3,0.5)) + 
    xlab('Days') + 
    ylab('ACF') + 
    theme_bw() + 
    theme(text = element_text(size = 15))

acf_p

ggsave(filename = 'plots/fig4b.jpeg',plot = acf_p,
       width = 6,height = 3,units = 'in')
    

 ####### peak lag in every time series ################
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

n.cols = 7
palettes = colorRampPalette(brewer.pal(n.cols,'YlOrBr'))
colors = palettes(n.cols)

peak_day_p = 
    ggplot(peak_day_ts) +
    geom_point(aes(x = end_year,y = peak_lag,color = year_number,size = peak_acf)) +
    # scale_color_manual(values = c('4' = 'orange','5' = 'firebrick','6' = 'forestgreen','7' = 'navy')) +
    scale_color_manual(values = colors[3:7]) + 
    scale_x_continuous(breaks = unique(peak_day_ts$end_year)) +
    xlab("End Year")+
    ylab("Peak Day")+ 
    theme_bw() +
    theme(text = element_text(size = 15),
          axis.text = element_text(size = 15),
          # plot.tag = element_text(size = 30),
          # plot.tag.position = c(0.13,1.1),
          ) +
    guides(color = guide_legend(title = "Number of Years\nof Data Included",
                                title.vjust = 1,
                                override.aes = list(size = 5)),
           size = guide_legend(title = "ACF Value"))
ggsave(filename = 'plots/fig4c.jpeg',plot = peak_day_p,
       width = 6,height = 3,units = 'in')

peak_day_p


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

spec_p = 
    ggplot(data = cycle_spec) + 
    geom_col(aes(x = cycle, y = spec), width = 2) +
    geom_point(data = cycle_spec[which.max(cycle_spec$spec),],
               aes(x = cycle, y = spec),
               size = 2) + 
    scale_x_continuous(limits = c(0,500)) + 
    xlab('Days') + 
    ylab('Spectral Density') + 
    theme_bw() + 
    theme(text = element_text(size = 15),
          axis.text = element_text(size = 15),
          # plot.tag = element_text(size = 30),
          # plot.tag.position = c(0.13,1.1),
    )
ggsave(filename = 'plots/fig4d.jpeg',plot = spec_p,
       width = 6,height = 3,units = 'in')
    


ps = ggarrange(ilip_p,
          ggarrange(acf_p,peak_day_p,spec_p,
                    nrow = 3,
                    labels = c('B','C','D')),
          ncol = 2,
          labels = 'A')

# ps

ggsave(filename = 'plots/Figure4.jpeg',
       plot = ps,
       width = 12,
       height = 10,
       units = 'in')



# cycle_spec$cycle = as.character(cycle_spec$cycle)
# cycle_spec = cycle_spec[order(cycle_spec$spec,decreasing = T),]




