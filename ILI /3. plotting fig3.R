rm(list = ls())
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(cowplot)
library(ggpubr)

setwd("~/Dropbox/Influenza and Respiratory Disease in Vietnam/code")
source('functions/Multiple_stepfunc_v5.R')
load('Rdata/7d.smth.zeta.ilipercent.Rdata')

############## plot zeta step function grid search ##################

####### Read all the fits #######
all_fits_pt1 = list()
for(i in 1:100) {
    
    load(paste0('Rdata/Results from Roar/Roar_zeta_minAvailDay_150_sensitivity_analysis_1_100/best_fit_',i,'.RData'))
    all_fits_pt1[[i]] = best_fit
}

all_fits_pt2 = list()
for(i in 1:95) {
    
    load(paste0('Rdata/Results from Roar/Roar_zeta_minAvailDay_150_sensitivity_analysis_101_195/best_fit_',i,'.RData'))
    all_fits_pt2[[i]] = best_fit
}

all_fits_pt3 = list()
for(i in 1:40) {
    
    load(paste0('Rdata/Results from Roar/Roar_zeta_minAvailDay_150_sensitivity_analysis_cyc150_189/best_fit_',i,'.RData'))
    all_fits_pt3[[i]] = best_fit
}

all_fits_pt4 = list()
for(i in 1:70) {
    
    load(paste0('Rdata/Results from Roar/Roar_zeta_minAvailDay_150_sensitivity_analysis_cyc385_450/best_fit_',i,'.RData'))
    all_fits_pt4[[i]] = best_fit
}

all_fits = c(all_fits_pt1,all_fits_pt2,all_fits_pt3,all_fits_pt4)

all_N = sapply(1:length(all_fits),function(x){all_fits[[x]]$N})
all_aic = sapply(1:length(all_fits),function(x){all_fits[[x]]$aic})
all_cycle = sapply(1:length(all_fits), function(x){all_fits[[x]]$cycle})
cycles = unique(all_cycle)

all_results = data.frame(cycle = all_cycle, N = all_N,aic = all_aic)
all_results[which.min(all_results$aic),]
all_results = all_results[order(all_results$aic),]


selected_results = all_results[all_results$aic < 0,]
p1 = 
    ggplot(selected_results) +
    geom_raster(aes(x = cycle,y = N,fill = aic)) +
    scale_fill_viridis_c(option = "H") +
    scale_x_continuous(breaks = seq(min(cycles),max(cycles),15)) +
    labs(fill = 'AIC') +
    xlab('Cycle Length') +
    ylab('Number of Steps') +
    theme_bw() +
    theme(text = element_text(size  = 15),
          axis.text = element_text(size = 13),
          # plot.margin = margin(t = 55,b = 10,r = 5,l = 0),
          # plot.tag = element_text(size = 35),
          # plot.tag.position = c(0.05,1.1)
          )
ggsave(filename = 'plots/fig3a.jpeg',
       plot = p1, 
       width = 9,
       height = 5,
       units = 'in')
       
          


############# plot diff AIC of zeta ########

source('functions/get normal model.R')
source('functions/predictor importance.R')
load('Rdata/Dataset_zeta_minAvailDay_150_scaled_weather_human_AR7.RData')

##### Null model #####
null_model = lm(logzeta ~ 1,data = new_dataset)
AIC(null_model)


### Add cycle to the data ###
load('Rdata/zeta_best_stepfunc_for_sensitivity_analysis_cyc150_450.Rdata')
load('Rdata/7d.smth.zeta.ilipercent.Rdata')
sim_ts_365 = all_zeta_sim_ts[which(rownames(all_zeta_sim_ts) == '365'),]
sim_ts_210 = all_zeta_sim_ts[which(rownames(all_zeta_sim_ts) == '210'),]

new_dataset$sim_ts_cyc365 = sim_ts_365[which(zeta_perc_avg_df$Date %in% new_dataset$date)]
new_dataset$sim_ts_cyc210 = sim_ts_210[which(zeta_perc_avg_df$Date %in% new_dataset$date)]

#### full model ####
full_model = get_normal_model(response_name = 'logzeta',
                              data = new_dataset[,-1])

summary(full_model)
confint(full_model)

#### Get the importance of each predictor in the full model ####
all_preds = names(coef(full_model))[-1]

pred_impo = data.frame(predictor = all_preds,
                       aic_diff = NA)

for(i in 1:nrow(pred_impo)) {
    
    pred_impo$aic_diff[i] = diff_aic_backward(predictor_to_assess = all_preds[i],
                                              response_name = 'logzeta',
                                              full_model = full_model,
                                              data = new_dataset[,-1])
}


pred_impo = pred_impo[order(pred_impo$aic_diff),]
### improvement should be positive ##
pred_impo$aic_diff = -pred_impo$aic_diff

#### Add the total diff-AIC from all the annual covariates ###
# pred_impo = rbind(pred_impo,
#                   data.frame(predictor = 'annual_vars',
#                              aic_diff = sum(pred_impo$aic_diff[-which(pred_impo$predictor %in% c('zeta_lag7d','sim_ts_cyc210'))])))

#### change the predictor name ###
pred_impo$predictor_name = c('7-day lagged zeta','Non-annual cycle', 'Annual cycle',
                             'AH','14-day lagged AH',
                             '7-day lagged RF', '21-day lagged RF',
                             '21-day lagged temp','21-day lagged AH')

# pred_impo = pred_impo[order(pred_impo$aic_diff,decreasing = T),]

pred_position = pred_impo$predictor_name

#### generate figures ###
p2 = 
    ggplot(pred_impo) +
    geom_col(aes(x = predictor_name, y = aic_diff)) +
    scale_x_discrete(limits = pred_position,
                     labels = function(x){str_wrap(x,width = 2)}) +
    xlab('') +
    ylab('∆AIC') +
    # labs(tag = 'B') +
    theme_bw() + 
    theme(text = element_text(size  = 15),
          axis.text = element_text(size  = 15),
          # plot.margin = margin(t = 65,b = 10,r = 5,l = 0),
          # plot.tag = element_text(size = 35),
          # plot.tag.position = c(0.1,1.1)
    )

p2

ggsave(filename = 'plots/fig3b.jpeg',
       plot = p2, 
       width = 9,
       height = 5,
       units = 'in')



##### plot full model with zeta ###
####### Plot Full model only ############
new_dataset$fitted = exp(fitted(full_model))

### 95% prediction interval ###
pred_int = predict(full_model,newdata = new_dataset,
                   interval = 'prediction')

pred_int = exp(pred_int)

pred_int = data.frame(pred_int)
pred_int$date = new_dataset$date

p3 = 
    ggplot() + 
    geom_ribbon(aes(x = date, ymin = lwr, ymax = upr), 
                data = pred_int,
                fill = 'orange',alpha = 0.5) + 
    geom_line(aes(x = date, y = exp(logzeta)),color = 'darkgrey',data = new_dataset) +
    geom_line(aes(x = date, y = fitted),color = 'red',alpha = 0.6,data = new_dataset) +
   
    # geom_vline(xintercept = seq(as.Date('2010-01-01'),as.Date('2020-01-01'),"1 year"),
    #            linetype = 'dashed',color = 'grey') +
    scale_x_date(limits = c(as.Date('2010-01-01'),as.Date('2019-12-31')),
                 date_breaks = '1 year',
                 date_labels = '%Y',
                 expand = c(0,0)) +
    xlab('Date') + 
    ylab('ILI zeta score') + 
    theme_bw() + 
    theme(text = element_text(size  = 20))
    
p3
ggsave(filename = 'plots/fig3c.jpeg',
       plot = p3, 
       width = 18,
       height = 5,
       units = 'in')


ps = ggarrange(ggarrange(p1,p2, ncol = 2,labels = c("","B")),
          p3,
          nrow = 2,
          labels = c('A','C'))
ps
ggsave(filename = 'plots/Figure3.jpeg',
       plot = ps, 
       width = 18,
       height = 10,
       units = 'in')
