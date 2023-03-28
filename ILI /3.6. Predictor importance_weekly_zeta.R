####### predictor importance in the regression of weekly zeta score ############
source('functions/pred_impo_pipeline functions.R')
library(relaimpo)
load('Rdata/Weekly zeta from averaging daily zeta in a week.Rdata')

### get weekly weather data in HCMC ###
weekly_weather = get_weather_data(weather_csv_location = 'datasets/POWER_Point_Daily_20100101_20191231_010d7714N_106d6895E_LST_VN.csv',
                                  plot_raw = T,
                                  plot_weekly = T,
                                  rows_to_skip = 22,
                                  ncol = 16,
                                  save = F)

### get fitted stepfunc for the weekly zeta ###
weekly_stepfunc_info = get_stepfunc(results_location = 'Rdata/Results from Roar/Roar_weekly_zeta_hcmc/best_fit_',
                                    N_fits = 88,
                                    data_to_fit = zeta_weekly_avg_df$weekly_zeta,
                                    save = F)
all_sim_ts = weekly_stepfunc_info$all_sim_ts

annual_cycle = all_sim_ts[which(rownames(all_sim_ts) == '52'),]
non_annual_cycle = all_sim_ts[which(rownames(all_sim_ts) == '30'),]

### visual check ###
plot(zeta_weekly_avg_df$date, zeta_weekly_avg_df$weekly_zeta,type = 'l',main = 'HCMC')
lines(zeta_weekly_avg_df$date,annual_cycle,col = 'red')
lines(zeta_weekly_avg_df$date,non_annual_cycle,col = 'blue')

#### get regression dataset ###
weekly_zeta_reg_data = get_regression_data(weather_data = weekly_weather,
                                           ili_zeta_df = zeta_weekly_avg_df,
                                           time_name = 'date',
                                           response_name = 'weekly_zeta',
                                           annual_cycle = annual_cycle,
                                           non_annual_cycle = non_annual_cycle,
                                           save = F)

#### get the full model ###
full_model = get_normal_model(response_name = 'ili_zeta',
                              data = weekly_zeta_reg_data)
summary(full_model)

### plot data with model ###
plot(weekly_zeta_reg_data$date,weekly_zeta_reg_data$ili_zeta,type = 'l',main = 'HCMC')
lines(weekly_zeta_reg_data$date,fitted(full_model),col = 'red')


#### predictor importance in the full model ###
preds_aic = get_pred_impo_full_model(model = full_model,
                                     response_name = 'ili_zeta',
                                     family = 'gaussian',
                                     data = weekly_zeta_reg_data)
preds_aic$diff_aic = -preds_aic$diff_aic
preds_aic = preds_aic[order(preds_aic$diff_aic,decreasing = T),]
barplot(preds_aic$diff_aic,names.arg = preds_aic$preds,main = 'HCMC')


zeta_relaimpo = calc.relimp(full_model)
save(zeta_relaimpo,file = 'Rdata/rela_impo_zeta_reg_HCMC.Rdata')

zeta_lmg = zeta_relaimpo@lmg[order(zeta_relaimpo@lmg,decreasing = T)]

barplot(zeta_lmg,main = 'HCMC')

