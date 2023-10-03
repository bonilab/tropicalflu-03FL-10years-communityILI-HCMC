rm(list = ls())
library(reshape2)
library(ggplot2)
setwd("~/Dropbox/Influenza and Respiratory Disease in Vietnam/code")
#### Correlation of zeta score from each year,the last day in leap year is removed

load('Rdata/7d.smth.zeta.ilipercent.Rdata')

total_series = c()
for(i in 2010:2019){
    
    series = zeta_perc_avg_df$smoothed_zeta_score[zeta_perc_avg_df$Date >= as.Date(paste0(i,"-01-01"),tz = UTC) &
                                                   zeta_perc_avg_df$Date <= as.Date(paste0(i,"-12-31"),tz = UTC)]
    total_series = cbind(total_series,series)
    
}

colnames(total_series) = as.character(seq(2010,2019,1))

cor(total_series[,"2018"],total_series[,"2019"])

corr_zeta = matrix(NA,nrow = ncol(total_series),ncol = ncol(total_series))
for(i in 1:ncol(total_series)){
    for(j in 1:ncol(total_series)){
        
        if(j > i) {
            next()
        } else {
    corr_zeta[i,j] = round(cor(total_series[,i],total_series[,j]),digits = 3)
        }
    }
}

rownames(corr_zeta) = as.character(seq(2010,2019,1))
colnames(corr_zeta) = as.character(seq(2010,2019,1))

image(corr_zeta,axes = FALSE)
axis(side = 1,at = seq(0,1,length.out = 10),labels = seq(2010,2019,1))

molten_corr_zeta = melt(corr_zeta)

ggplot(molten_corr_zeta) +
    geom_raster(aes(x = Var1,y = Var2,fill = value)) +
    geom_text(aes(x = Var1,y = Var2,label = value),size = 3.5) +
    scale_x_continuous(breaks = seq(2010,2019,1),expand = c(0,0)) +
    scale_fill_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0,
                         na.value = 'white')+
    #scale_fill_viridis_c(option = 'D') +
    scale_y_continuous(breaks = seq(2010,2019,1),expand = c(0,0),
                       position = 'right') +
    ylab("")+
    xlab("")+
    # ggtitle("Correlation map of zeta score of ILI percent")+
    theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = '../paper/Figures/SuppFig3.jpeg',
       width = 8,height = 7,units = 'in')


