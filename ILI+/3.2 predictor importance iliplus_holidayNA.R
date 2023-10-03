#### predictor importance of ILI+ regression ######
rm(list = ls())
load('Rdata/Overall.ILIplus.21-day.aggregated.pcr.RData')
source('functions/predictor importance.R')
source('functions/get gamma hurdle model.R')
library(ggplot2)
library(cowplot)

#### cyc330 has the lowest AIC, load N-stepfunc from cyc330 cyc385 ###
load('Rdata/ILI+_best_stepfunc_for_sensitivity_analysis_cyc150_450.Rdata')
load('Rdata/Data of overall iliplus_UV_WS_absent_scaled_weather.AR_21.human for gamma hurdle model.Rdata')
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
confint(full_gamma_hurdle)
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
                            preds_name = c('21-day lagged ILI+','Temp','21-day lagged RF',
                                           '330-day cycle',
                                           '14-day lagged AH',
                                           '7-day lagged RF',
                                           'School term','385-day cycle',
                                           'RF'),
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
                              pred_names = c('21-day lagged ILI+','21-day lagged Temp','330-day cycle',
                                             '385-day cycle','7-day lagged AH','21-day lagged RF',
                                             '7 day lagged RF','School term','RF','21-day lagged AH'),
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

### ggplot version ###
p_fit = ggplot() + 
    geom_line(aes(x = date, y = smth_7d_iliplus),color = 'grey',data = data_for_fit) + 
    geom_line(aes(x = date, y = fitted),color = 'red',data = data.frame(date = data_for_fit$date,
                                                                        fitted = fitted(full_gamma_hurdle))) + 
    scale_x_date(date_breaks = '1 year',date_labels = '%Y') +
    theme_bw() + 
    xlab('') + 
    ylab('') + 
    theme(text = element_text(size = 20))


p_fit
ggsave(filename = 'plots/fig5c.jpeg',
       plot = p_fit,
       width = 18,height = 5,units = 'in')







### baseplot version ###
jpeg(filename = 'plots/ILIplus_fitted_full_model.jpg',
     width = 2000, height = 500, units = 'px')
plot(data_for_fit$date,
     data_for_fit$smth_7d_iliplus,
     type = "l",lwd = 1.5,
     xlim = c(as.Date('2012-01-01'),as.Date('2019-12-31')),
     ylim = c(0,1),
     main = '',xlab = '',ylab = '',
     xaxs = "i",
     xaxt = "n",
     col = 'grey',
     cex.axis = 2)
axis(1,at = seq(as.Date('2012-01-01'),as.Date('2020-01-01'),"1 year"),
     labels = format(seq(as.Date('2012-01-01'),as.Date('2020-01-01'),"1 year"),"%Y"),
     cex.axis = 2)
abline(v = seq(as.Date('2012-01-01'),as.Date('2020-01-01'),"1 year"),
       lty = "dashed",
       col = "grey")
lines(data_for_fit$date,fitted(full_gamma_hurdle),col = 'red',lwd = 1.5)
legend('topright',
       legend = c('Observed ILI+','Predicted ILI+'),
       lty = c('solid','solid'),
       col = c('grey','red'),
       bty = 'n',
       cex = 1.7)

# mtext('C',side = 3,adj = 0,outer = F,line = 2,cex = 3)
dev.off()
