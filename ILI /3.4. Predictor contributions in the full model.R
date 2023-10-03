####### Predictor importance on lognormal model ##########
rm(list = ls())
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

#### change the predictor name ###
pred_impo$predictor_name = c('7-day lagged zeta','Non-annual Cycle', 'Annual cycle',
                             'AH','14-day lagged AH',
                             '7-day lagged RF', '21-day lagged RF',
                             '21-day lagged temp','21-day lagged AH')

#### Add the total diff-AIC from all the annual covariates ###
# pred_impo = rbind(pred_impo,
#                   data.frame(predictor = 'annual_vars',
#                              aic_diff = sum(pred_impo$aic_diff[-which(pred_impo$predictor %in% c('zeta_lag7d','sim_ts_cyc210'))]),
#                              predictor_name = 'Annual\ncovariates'))

pred_position = pred_impo$predictor_name

library(ggplot2)
library(stringr)

ggplot(pred_impo) +
        geom_col(aes(x = predictor_name, y = aic_diff)) +
        scale_x_discrete(limits = pred_position,
                         labels = function(x){str_wrap(x,width = 2)}) +
    xlab('') +
    ylab('∆AIC') +
    labs(tag = 'B') +
    theme(text = element_text(size = 25),
          axis.text = element_text(size = 20),
          plot.margin = margin(t = 65,b = 10,r = 5,l = 0),
          plot.tag = element_text(size = 35),
          plot.tag.position = c(0.1,1.1))
# ggsave(filename = '../plots/weather data analysis/annotated_predictor_importance_full_model_zeta_sigma_unfixed.jpg',
#        width = 3000, height = 2200, units = c('px'))


### relative importance ###
library(relaimpo)

relaimpo_zeta = calc.relimp(full_model)
relaimpo_zeta@lmg
barplot(relaimpo_zeta@lmg)

### hierarchical partitioning ###
library(hier.part)
rela_impo_zeta = hier.part(y = new_dataset$logzeta,
                           xcan = new_dataset[,which(colnames(new_dataset) %in% names(coef(full_model))[-1])],
                           family = 'gaussian',
                           gof = 'logLik')

I.perc_relaimpo = rela_impo_zeta$I.perc
I.perc_relaimpo_df = data.frame(pred = rownames(I.perc_relaimpo),
                                I.perc = I.perc_relaimpo)

I.perc_relaimpo_df = I.perc_relaimpo_df[order(I.perc_relaimpo_df$ind.exp.var,
                                              decreasing = T),]

barplot(I.perc_relaimpo_df$ind.exp.var,
        names.arg = I.perc_relaimpo_df$pred)


##### plot full model with zeta ###
####### Plot Full model only ############
jpeg(filename = '../plots/weather data analysis/7d.smth_zeta_with_full_model.jpg',
     width = 1500, height = 500, units = 'px')
plot(new_dataset$date,
     exp(new_dataset$logzeta),
     type = "l",lwd = 1.5,
     xlim = c(as.Date('2010-01-01'),as.Date('2019-12-31')),
     ylim = c(0.5,1.55),
     main = '',xlab = '',ylab = '',
     xaxs = "i",
     xaxt = "n",
     cex.axis = 1.7)
axis(1,at = seq(as.Date('2010-01-01'),as.Date('2020-01-01'),"1 year"),
     labels = format(seq(as.Date('2010-01-01'),as.Date('2020-01-01'),"1 year"),"%Y"),
     cex.axis = 1.7)
abline(v = seq(as.Date('2010-01-01'),as.Date('2020-01-01'),"1 year"),
       lty = "dashed",
       col = "grey")
lines(new_dataset$date,exp(fitted(full_model)),col = 'red',lwd = 1.5)
mtext('C',side = 3,adj = 0,outer = F,line = 1,cex = 2)
dev.off()
