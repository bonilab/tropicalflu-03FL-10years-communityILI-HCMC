library(ggpubr)
library(cowplot)
library(gridGraphics)
library(ggplot2)
library(RColorBrewer)

rm(list = ls())
load("Rdata/7d.smth.zeta.ilipercent.Rdata")
source('functions/Detrend and Smooth.R')

########### ACF of 7-day smoothed zeta score #################

acf_smoothed_zeta_score = acf(zeta_perc_avg_df$smoothed_zeta_score,lag.max = 1095,col = "grey",
                              na.action = na.pass,
                              ylim = c(-0.2,0.5),
                              xaxt = 'n',yaxt = 'n',
                              xlab = '',ylab = '',
                              main = "",
                              cex.lab = 1.2)
title(xlab = 'Lag',
      ylab = 'ACF',
      line = 2,cex.main = 1.5)

### The lag of acf starts from 0
acf_lag = cbind(acf_smoothed_zeta_score$acf,acf_smoothed_zeta_score$lag)
max_lag = acf_smoothed_zeta_score$lag[which.max(acf_smoothed_zeta_score$acf[150:400])] + 149
max_lag
max_cycle = seq(max_lag,1095,max_lag)
max_cycle
axis(1,at = max_cycle,labels = max_cycle,cex.axis = 1)
axis(2, at = seq(-0.2,0.5,0.2),labels = seq(-0.2,0.5,0.2),cex.axis = 1)
abline(v = max_cycle,lty = "dashed",col = "black")
points(x = seq(365,1095,365),
       y = c(acf_smoothed_zeta_score$acf[365],acf_smoothed_zeta_score$acf[730],acf_smoothed_zeta_score$acf[1095]),
       pch = 20,cex = 2)

########### Multi-panel graph of ACF of zeta score of ILI PERCENT #############

annual = seq(1,365*4,365)
acf_ps = list()
k = 1

for (i in 2010:2014){
    for (j in 2015:2019){
        
        
        if (j - i < 5) {
            
            acf_ps[[k]] = ggplot() + theme_void()
            k = k + 1
            
        } else {
            start = as.Date(paste0(i,"-01-01"),tz = "UTC")
            end = as.Date(paste0(j,"-12-31"),tz = "UTC")
            
            # timeseries = zeta_perc_avg_df$smoothed_zeta_score[which(zeta_perc_avg_df$Date >= start & zeta_perc_avg_df$Date < end)]
            timeseries = zeta_perc_avg_df$smoothed_zeta_score[which(zeta_perc_avg_df$Date >= start & zeta_perc_avg_df$Date < end)]
            
            acf_temp = acf(timeseries,lag.max = 365*3,col = "grey",
                           na.action = na.pass,
                           ylim = c(-0.2,0.5),
                           xlab = '',ylab = '',
                           xaxt = 'n',yaxt = 'n',
                           main = '',
                           plot = F)
            
            max_lag = acf_temp$lag[which.max(acf_temp$acf[150:450])] + 150
            max_cycle = seq(max_lag,1095,max_lag)
            
            acf_df = data.frame(lag = acf_temp$lag,
                                acf = acf_temp$acf)
            
            acf_ps[[k]] = 
                ggplot(acf_df) + 
                geom_area(aes(x = lag, y = acf),fill = 'grey') + 
                geom_point(aes(x = annual, y = acf),
                           data = data.frame(annual = seq(365,1095,365), 
                                             acf = acf_df$acf[which(acf_df$lag %in% seq(365,1095,365))]),
                           size = 2) + 
                geom_vline(xintercept = max_cycle,linetype = 'dashed') +
                geom_hline(yintercept = c(-1.96/sqrt(length(timeseries)),
                                          1.96/sqrt(length(timeseries))),
                           linetype = 'dashed',
                           color = 'blue') + 
                scale_x_continuous(breaks = max_cycle) + 
                scale_y_continuous(limits = c(-0.3,0.5)) + 
                xlab('') + 
                ylab('') + 
                annotate('text', x = 1030,y = 0.4,label = paste0(i,'-\n',j),
                         size = 5) + 
                theme_bw() + 
                theme(text = element_text(size = 15))
            
            # if(i == 2010 & j == 2015) {
            #     
            #     acf_ps[[k]] + 
            #         labs(tag = 'A') +
            #         theme(plot.tag = element_text(size = 20),
            #               plot.tag.position = c(0.16,1.1))
            # }
            
            k = k + 1
        }
    }
}


all_acf_ps = 
    plot_grid(acf_ps[[1]],  acf_ps[[2]],  acf_ps[[3]],  acf_ps[[4]],  acf_ps[[5]],  acf_ps[[6]], 
              acf_ps[[7]],  acf_ps[[8]],  acf_ps[[9]],  acf_ps[[10]], acf_ps[[11]], acf_ps[[12]],
              acf_ps[[13]], acf_ps[[14]], acf_ps[[15]], acf_ps[[16]], acf_ps[[17]], acf_ps[[18]],
              acf_ps[[19]], acf_ps[[20]], acf_ps[[21]], acf_ps[[22]], acf_ps[[23]], acf_ps[[24]],
              acf_ps[[25]],
              nrow = 5, ncol = 5)

all_acf_ps = plot_grid(all_acf_ps, labels = 'A')

all_acf_ps


####### Separate time series as 5 - 10 years, and calculate their ACF ##########
multiple_acf_list = list()
j = 1

for(ts_length in 5:9){
    
    all_ts_acf = c()
    for(start_year in 2010:(2019- ts_length)){
        
        target_ts = seq(as.Date(paste0(start_year,"-01-01"),tz = "UTC"),
                        as.Date(paste0(start_year + ts_length,"-12-31"),tz = "UTC"),"1 day")
        target_zeta = zeta_perc_avg_df$smoothed_zeta_score[which(zeta_perc_avg_df$Date %in% target_ts)]
        
        target_ts_acf = acf(target_zeta,lag.max = 1095,plot = F,
                            na.action = na.pass)$acf
        all_ts_acf = cbind(all_ts_acf,target_ts_acf)
        
    }
    
    colnames(all_ts_acf) = mapply(paste0,2010:(2019-ts_length),"-",2010:(2019-ts_length)+ts_length)
    multiple_acf_list[[j]] = all_ts_acf
    j = j + 1
}

# rm(list = setdiff(ls(),c("multiple_acf_list","multiple_zeta_ts",
#                          "zeta_perc_avg_df")))


############# Get the peak ACF for annual and non-annual cycle #################
multiple_acf_mat = do.call(cbind, multiple_acf_list)
former_year = as.numeric(gsub("^(\\d{4})-(\\d{4})","\\1",dimnames(multiple_acf_mat)[[2]]))
followed_year = as.numeric(gsub("^(\\d{4})-(\\d{4})","\\2",dimnames(multiple_acf_mat)[[2]]))


##### Extract ACF when lag = 365
acf_annual = data.frame(ts = dimnames(multiple_acf_mat)[[2]],
                        startyear = former_year,
                        endyear = followed_year,
                        years = as.character(followed_year - former_year + 1),
                        acf_365 = multiple_acf_mat[365 + 1,],
                        cycle = "annual cycle")

#### Extract the peak ACF around 200 days
peak_acf_200 = sapply(1:ncol(multiple_acf_mat), function(x) {
    
    max(multiple_acf_mat[180:220,x])
    
})

#### Extract the peak lag around 200 days
peak_day_200 = sapply(1:ncol(multiple_acf_mat),function(x) {
    
    which.max(multiple_acf_mat[190:220,x]) + 190 - 1
})

peak_day_200

#### Summarize the ACF and peak lag around 200 days
acf_non_annual = data.frame(ts = dimnames(multiple_acf_mat)[[2]],
                            startyear = former_year,
                            endyear = followed_year,
                            years = as.character(followed_year - former_year + 1),
                            acf_200 = peak_acf_200,
                            cycle = "non-annual cycle")

acf_annual$years = factor(acf_annual$years,levels = c('6','7','8','9','10'))
acf_non_annual$years = factor(acf_non_annual$years,levels = c('6','7','8','9','10'))

acf_annual$years = as.numeric(acf_annual$years)
acf_non_annual$years = as.numeric(acf_non_annual$years)

n.cols = 7
palettes = colorRampPalette(brewer.pal(n.cols,'YlOrBr'))
colors = palettes(n.cols)

acf_annual$years = as.character(acf_annual$years)
acf_non_annual$years = as.character(acf_non_annual$years)

p = ggplot() +
    geom_point(aes(x = endyear,y = acf_365,color = years,shape = cycle),
               data = acf_annual,size = 5,alpha = 0.7)+ 
    geom_point(aes(x = endyear,y = acf_200,color = years,shape = cycle),
               data = acf_non_annual,size = 5,alpha = 0.7) +
    geom_hline(yintercept = c(-0.04,0.04),linetype = "dashed",
               color = "navy")+
    # scale_color_manual(values = c('6' = 'firebrick','7' = 'forestgreen','8' = 'navy',
    #                               '9' = 'orange','10' = 'black')) +
    # scale_color_brewer(palette = 'YlOrBr',values = 0.2) + 
    scale_color_manual(values = colors[3:7]) + 
    scale_shape_manual(values = c("annual cycle" = 16,"non-annual cycle" = 17),
                       labels = c('Annual Cycle','Non-annual Cycle')) +
    xlab("End Year")+
    ylab("ACF")+ 
    theme_bw() +
    theme(text = element_text(size = 15),
          axis.text = element_text(size = 15),
          plot.margin = margin(t = 55,r = 10,b = 10, l = 10),
          legend.position = 'bottom') +
    guides(color = guide_legend(title = "Number of Years of\nData Included",
                                title.vjust = 1,
                                nrow = 2),
           shape = guide_legend(title = "Cycle",
                                nrow = 2)) +
    coord_cartesian(clip = 'off')


p = plot_grid(p, labels = 'B')
p

# ggsave(filename = 'plots/shift_cycle_acf_zeta.jpg',
#        plot = p, 
#        width = 7,
#        height = 6,
#        units = 'in')



############# Periodogram ##############

### choose the periodogram without smoothing
# jpeg(filename = 'plots/SuppFig4.jpg',width = 6,height = 4,units = 'in',
#      res = 300)
zeta_perc_spec = spectrum(na.omit(zeta_perc_avg_df$smoothed_zeta_score),
                          ylab = "Amplitude",xlab = "Frequency",
                          main = "Smoothed periodogram of Zeta Score",
                          fast = F,log = 'no')
which.max(zeta_perc_spec$spec)
1/zeta_perc_spec$freq[which.max(zeta_perc_spec$spec)]

cycle_spec = cbind(1/zeta_perc_spec$freq,zeta_perc_spec$spec)
cycle_spec = cycle_spec[order(cycle_spec[,2],decreasing = T),]

# par(oma = c(0.5,0.5,1,0.5))
plot(x = 1/zeta_perc_spec$freq ,y = zeta_perc_spec$spec ,type = "h",xlim = c(0,500),
     xlab = "",ylab = "",main = "",
     lwd = 2,mgp = c(2,1,0),cex.axis = 1.5)
# mtext("A",side = 3,outer = F,line = 2,adj = 0,cex = 2)

#### The two cycles with the highest spectral density ##
points(x = cycle_spec[1,1], y = cycle_spec[1,2],pch = 20,cex = 1.5)
points(x = cycle_spec[2,1], y = cycle_spec[2,2],pch = 20,cex = 1.5)
# dev.off()


