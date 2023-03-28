######### Reading N-stepfunc fits of overall ILI+ from Roar #####
library(ggplot2)
library(ggpubr)
rm(list = ls())
setwd("~/Dropbox/Influenza and Respiratory Disease in Vietnam/code")
source('functions/Multiple_stepfunc_v5.R')
load('Rdata/Overall.ILIplus.21-day.aggregated.pcr_holidayNA.RData')

all_fits_pt1 = list()
for(i in 1:100) {
    
    load(paste0('Rdata/Results from Roar/Roar_updated_iliplus_minAD150_sensitivity_analysis_pt1/best_fit_',i,'.RData'))
    all_fits_pt1[[i]] = best_fit
}

all_fits_pt2 = list()
for(i in 1:100) {
    
    load(paste0('Rdata/Results from Roar/Roar_updated_iliplus_minAD150_sensitivity_analysis_pt2/best_fit_',i,'.RData'))
    all_fits_pt2[[i]] = best_fit
}

all_fits_pt3 = list()
for(i in 1:105) {
    
    load(paste0('Rdata/Results from Roar/Roar_updated_iliplus_minAD150_sensitivity_analysis_pt3/best_fit_',i,'.RData'))
    all_fits_pt3[[i]] = best_fit
}

all_fits = c(all_fits_pt1,all_fits_pt2,all_fits_pt3)

all_N = sapply(1:length(all_fits),function(x){all_fits[[x]]$N})
all_aic = sapply(1:length(all_fits),function(x){all_fits[[x]]$aic})
all_cycle = sapply(1:length(all_fits), function(x){all_fits[[x]]$cycle})
cycles = unique(all_cycle)

all_results = data.frame(cycle = all_cycle, N = all_N,aic = all_aic)
all_results = all_results[order(all_results$aic),]
all_results[which.min(all_results$aic),]

valid_results = all_results[all_results$aic < 0,]

library(ggplot2)
p1 = ggplot(valid_results) +
    geom_raster(aes(x = cycle,y = N,fill = aic)) +
    scale_fill_viridis_c(option = "H") +
    scale_x_continuous(breaks = seq(min(cycles),max(cycles),15)) +
    labs(fill = 'AIC') +
    xlab('Cycle Length') +
    ylab('Number of Steps') +
    theme_bw() +
    theme(text = element_text(size = 15),
          axis.text = element_text(size = 13))
p1

ggsave(filename = 'plots/fig5a.jpeg',
       plot = p1,
       width = 9, 
       height = 5, 
       units = 'in')



########## plot predictor importance #############
load('Rdata/Overall.ILIplus.21-day.aggregated.pcr_holidayNA.RData')
source('functions/predictor importance.R')
source('functions/get gamma hurdle model.R')
library(ggplot2)
library(cowplot)

#### cyc330 has the lowest AIC, load N-stepfunc from cyc330 cyc385 ###
load('Rdata/ILI+_holidayNA_best_stepfunc_for_sensitivity_analysis_cyc150_450.Rdata')
load('Rdata/Data of overall iliplus_holidayNA.UV_WS_absent_scaled_weather.AR_21.human for gamma hurdle model.Rdata')
#ili_plus_21d_df$sim_ts_385 = all_iliplus_sim_ts[which(rownames(all_iliplus_sim_ts) == '385'),]
ili_plus_21d_df$sim_ts_330 = all_iliplus_sim_ts[which(rownames(all_iliplus_sim_ts) == '330'),]
#ili_plus_21d_df$sim_ts_365 = all_iliplus_sim_ts[which(rownames(all_iliplus_sim_ts) == '365'),]
ili_plus_21d_df$sim_ts_385 = all_iliplus_sim_ts[which(rownames(all_iliplus_sim_ts) == '385'),]


data_for_fit$sim_ts_330 = ili_plus_21d_df$sim_ts_330[which(ili_plus_21d_df$date %in% data_for_fit$date)]
#data_for_fit$sim_ts_365 = ili_plus_21d_df$sim_ts_365[which(ili_plus_21d_df$date %in% data_for_fit$date)]
data_for_fit$sim_ts_385 = ili_plus_21d_df$sim_ts_385[which(ili_plus_21d_df$date %in% data_for_fit$date)]

#### weather human AR21 cyc330 and cyc365 cyc385 ###
full_gamma_hurdle = fit_gamma_hurdle_model(data = data_for_fit[,-which(colnames(data_for_fit) %in% c('date'))],
                                           response_name = 'smth_7d_iliplus',
                                           interaction = F)
summary(full_gamma_hurdle)
AIC(full_gamma_hurdle)

#### predictor importance of logistic regression #####

### need to rebuild logistic model based on full gamma hurdle model ###
data_for_logi = data_for_fit

data_for_logi$present_flu = ifelse(data_for_logi$smth_7d_iliplus !=0,1,0)

data_for_logi$smth_7d_iliplus = NULL

call_logi = as.character(full_gamma_hurdle$call[5])

null_logi_model = glm(present_flu~1,family = 'binomial',
                      data = data_for_logi)
AIC(null_logi_model)

full_logi_model = glm(paste0('present_flu',call_logi),
                      family = 'binomial',
                      data = data_for_logi)

summary(full_logi_model)
AIC(full_logi_model)

#### predictor importance of logistic model #####

#### adding non-linear transformation,sqrt(AR21) ####
logi_pred_impo = data.frame(preds = names(coef(full_logi_model))[-1],
                            preds_name = c('21-day lagged ILI+','Temp',
                                           '21-day lagged RF','385-day cycle',
                                           'School term', '7-day lagged AH',
                                           'RF','330-day cycle',
                                           '7-day lagged RF'),
                            aic_diff = NA)

for(i in 1:nrow(logi_pred_impo)) {
    
    logi_pred_impo$aic_diff[i] = diff_aic_backward_glm(predictor_to_assess = logi_pred_impo$preds[i],
                                                       response_name = 'present_flu',
                                                       full_model = full_logi_model,
                                                       family = 'binomial',
                                                       data = data_for_logi)
}

logi_pred_impo = logi_pred_impo[order(logi_pred_impo$aic_diff),]
logi_pred_impo$aic_diff = -logi_pred_impo$aic_diff
logi_pred_impo$model_name = 'Logistic model'
logi_pred_impo

preds_position_logi = logi_pred_impo$preds_name

library(ggplot2)
library(stringr)
library(ggpubr)

ggplot(logi_pred_impo) +
    facet_grid(.~ model_name) +
    geom_col(aes(x = preds,y = aic_diff))

p_logi = ggplot(logi_pred_impo) +
    facet_grid(.~ model_name) +
    geom_col(aes(x = preds_name, y = aic_diff)) +
    scale_x_discrete(limits = preds_position_logi,
                     labels = function(x){str_wrap(x,width = 2)}) +
    xlab('') +
    ylab('∆AIC') +
    # labs(tag = 'B') +
    theme_bw() + 
    theme(text = element_text(size = 15),
          axis.text = element_text(size = 15))

# plot.margin = margin(t = 50,b = 2,r = 2,l = 2),
# plot.tag = element_text(size = 40),
# plot.tag.position = c(0.09,1.18))

p_logi
ggsave(filename = 'plots/fig5b1.jpeg',
       plot = p_logi,
       width = 9, height =2.5, units = 'in')

###### predictor importance of gamma model #######
data_for_gamma = data_for_fit[-which(data_for_fit$smth_7d_iliplus == 0),]

gamma_call = as.character(full_gamma_hurdle$call[2])
gamma_full_model = glm(gamma_call,family = Gamma('log'),
                       data = data_for_gamma)
summary(gamma_full_model)

# preds_impo_gamma = data.frame(preds = names(coef(gamma_full_model))[-1],
#                               pred_names = c('21-day lagged ILI+','21-day lagged Temp','330-day cycle',
#                                              '5th root(21-day lagged RF)','385-day cycle',
#                                              '5th root(7-day lagged AH)','21-day lagged AH',
#                                              'School term','7-day lagged RF',
#                                              'RF'),
#                               aic_diff = NA)

preds_impo_gamma = data.frame(preds = names(coef(gamma_full_model))[-1],
                              pred_names = c('330-day cycle','385-day cycle',
                                             '21-day lagged ILI+','21-day lagged temp',
                                             'RF','School term','AH','7-day lagged RF',
                                             '21-day lagged RF','Temp'),
                              aic_diff = NA)

for(i in 1:nrow(preds_impo_gamma)) {
    
    preds_impo_gamma$aic_diff[i] = diff_aic_backward_glm(predictor_to_assess = preds_impo_gamma$preds[i],
                                                         response_name = 'smth_7d_iliplus',
                                                         full_model = gamma_full_model,
                                                         family = Gamma('log'),
                                                         data = data_for_gamma)
    
}

preds_impo_gamma = preds_impo_gamma[order(preds_impo_gamma$aic_diff),]
preds_impo_gamma$aic_diff = -preds_impo_gamma$aic_diff
preds_impo_gamma$model_name = 'Gamma model'
preds_impo_gamma

preds_position_gamma = preds_impo_gamma$pred_names

p_gamma = ggplot(preds_impo_gamma) +
    facet_grid(.~ model_name) +
    geom_col(aes(x = pred_names, y = aic_diff)) +
    scale_x_discrete(limits = preds_position_gamma,
                     labels = function(x){str_wrap(x,width = 2)}) +
    xlab('') +
    ylab('∆AIC') +
    theme_bw() + 
    theme(text = element_text(size = 15),
          axis.text = element_text(size = 15))
# plot.margin = margin(t = 2,b = 2,r = 2,l = 2))

p_gamma
ggsave(filename = 'plots/fig5b2.jpeg',
       plot = p_gamma,
       width = 9, height =2.5, units = 'in')


##### plot full model only with data ######
load('Rdata/ILI+_param_boostrap.Rdata')
### ggplot version ###
p_fit = ggplot() + 
    geom_ribbon(aes(x = date, ymin = lower, ymax = upper),
                data = fitted_data,
                fill = 'orange',alpha = 0.5) +
    geom_line(aes(x = date, y = smth_7d_iliplus),color = 'grey',data = data_for_fit) + 
    geom_line(aes(x = date, y = fitted),color = 'red',alpha = 0.6,
              data = data.frame(date = data_for_fit$date,
                                fitted = fitted(full_gamma_hurdle))) + 
    scale_x_date(date_breaks = '1 year',date_labels = '%Y') +
    theme_bw() + 
    xlab('Date') + 
    ylab('ILI+') + 
    theme(text = element_text(size = 15))


p_fit
ggsave(filename = 'plots/fig5c.jpeg',
       plot = p_fit,
       width = 18,height = 5,units = 'in')


ps = ggarrange(ggarrange(p1,labels = 'A',
                    ggarrange(p_logi,p_gamma,nrow = 2,labels = c('B','')),
                    ncol = 2),
          p_fit,
          nrow = 2,
          labels = c('','C'))

ps

ggsave(filename = 'plots/Figure5.jpeg',
       plot = ps, 
       width = 18,
       height = 10,
       units = 'in')

