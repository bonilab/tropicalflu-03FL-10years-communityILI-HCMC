###### plot the models with/without the cycles #####
rm(list = ls())
source('functions/get normal model.R')
load('Rdata/Dataset_zeta_minAvailDay_150_scaled_weather_human_AR7.RData')

# ### weather fit ###
weather_predictor_names = c(paste('temperature_2m',c('lag_7d','lag_14d','lag_21d'),sep = '_'),
                            paste('AH',c('lag_7d','lag_14d','lag_21d'),sep = '_'),
                            paste('RF',c('lag_7d','lag_14d','lag_21d'),sep = '_'),
                            paste('WS2M',c('lag_7d','lag_14d','lag_21d'),sep = '_'))
# 
fit_weather = get_normal_model(response_name = 'logzeta',
                               data = new_dataset[,which(colnames(new_dataset) %in% c(weather_predictor_names,'logzeta'))])

##### Null model #####
null_model = lm(logzeta ~ 1,data = new_dataset)
AIC(null_model)

### The model is different than the one evaluating contribution of cycles ##
## AIC-based model selection is added to get the best model ##
fit_weather_human_ar7 = get_normal_model(response_name = 'logzeta',
                                         data = new_dataset[,-1])
summary(fit_weather_human_ar7)
AIC(fit_weather_human_ar7)

### Add cycle to the data ###
load('Rdata/zeta_best_stepfunc_for_sensitivity_analysis_cyc150_450.Rdata')
load('Rdata/7d.smth.zeta.ilipercent.Rdata')
sim_ts_365 = all_zeta_sim_ts[which(rownames(all_zeta_sim_ts) == '365'),]
sim_ts_210 = all_zeta_sim_ts[which(rownames(all_zeta_sim_ts) == '210'),]

new_dataset$sim_ts_cyc365 = sim_ts_365[which(zeta_perc_avg_df$Date %in% new_dataset$date)]
new_dataset$sim_ts_cyc210 = sim_ts_210[which(zeta_perc_avg_df$Date %in% new_dataset$date)]


#### Model with cycle ###
formula_base_model = paste(names(coef(fit_weather_human_ar7))[-1],
                           collapse = '+')

### non-annual cycle only ###
fit_weather_human_ar7_cyc210 = lm(formula = paste('logzeta ~ ',formula_base_model,'+sim_ts_cyc210'),
                                  data = new_dataset)
AIC(fit_weather_human_ar7_cyc210)

fit_weather_human_ar7_cyc365 =  lm(formula = paste('logzeta ~ ',formula_base_model,'+sim_ts_cyc365'),
                                   data = new_dataset)
AIC(fit_weather_human_ar7_cyc365)

full_model = lm(formula = paste('logzeta ~ ',formula_base_model,'+sim_ts_cyc365+sim_ts_cyc210'),
                data = new_dataset)
summary(full_model)
AIC(full_model)

#### Add cycles and AIC-based selection ###
fit_210 = get_normal_model(response_name = 'logzeta',
                           data = new_dataset[,-which(colnames(new_dataset) %in% c('date','sim_ts_cyc365'))])
summary(fit_210)
AIC(fit_210)


fit_365 = get_normal_model(response_name = 'logzeta',
                           data = new_dataset[,-which(colnames(new_dataset) %in% c('date','sim_ts_cyc210'))])
summary(fit_365)
AIC(fit_365)


fit_365_210 = get_normal_model(response_name = 'logzeta',
                               data = new_dataset[,-which(colnames(new_dataset) %in% c('date'))])
summary(fit_365_210)
AIC(fit_365_210)

#### entire ts ######
models_list = list(fit_weather_human_ar7,fit_weather_human_ar7_cyc210,
                   fit_weather_human_ar7_cyc365,full_model)
colors = c('navy','forestgreen','orange','red')
model_legends = c('No Cycle','Non-annual Cycle',
                  'Annual Cycle','Both Cycles')

par(oma = c(1,1,3,1))

model_legends = c('No Cycle','Non-annual\nCycle',
                  'Annual Cycle','Both Cycles')
layout(mat = matrix(c(1,2,3,4,5),nrow = 5),
       heights = c(rep(1,4),0.5))

for(i in 1:length(models_list)) {
    par(mar = c(1,2.5,2,2))
    plot(new_dataset$date,exp(new_dataset$logzeta),
         type = "l",lwd = 1.5,
         xlim = c(as.Date('2010-01-01'),as.Date('2019-12-31')),
         ylim = c(0.5,1.6),
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
    lines(new_dataset$date,exp(fitted(models_list[[i]])),col = colors[i],lwd = 1.5)
    
}
par(mar = c(0,0,0,0))
plot(1,type = 'n',axes = F, xlab = '',ylab = '')
legend('center',
       legend = c('ILI Zeta Score',model_legends),
       col = c('black',colors),
       lty = 'solid',
       lwd = 3,
       # ncol = 2,
       horiz = T,
       cex = 1.5,
       inset = 0,
       text.width = 0.08,
       xpd = T,
       x.intersp = 0.5,
       bty = 'n')

mtext("A",side = 3,outer = T,adj = 0,cex = 2)

####### Plot Full model only ############
jpeg(filename = '../_ANALYSIS/plots/weather data analysis/annotated_full_model_fitted_zeta.jpg',
     width = 2000, height = 500, units = 'px')
plot(new_dataset$date,
     exp(new_dataset$logzeta),
     type = "l",lwd = 1.5,
     xlim = c(as.Date('2010-01-01'),as.Date('2019-12-31')),
     ylim = c(0.5,1.55),
     main = '',xlab = '',ylab = '',
     xaxs = "i",
     xaxt = "n",
     col = 'grey',
     cex.axis = 1.7)
axis(1,at = seq(as.Date('2010-01-01'),as.Date('2020-01-01'),"1 year"),
     labels = format(seq(as.Date('2010-01-01'),as.Date('2020-01-01'),"1 year"),"%Y"),
     cex.axis = 1.7)
abline(v = seq(as.Date('2010-01-01'),as.Date('2020-01-01'),"1 year"),
       lty = "dashed",
       col = "grey")
lines(new_dataset$date,exp(fitted(full_model)),col = 'red',lwd = 1.5)
mtext('C',side = 3,adj = 0,outer = F,line = 1,cex = 3)
dev.off()
