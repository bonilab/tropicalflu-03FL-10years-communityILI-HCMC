######## Estimate 2-step function for 7d.smth.ili+.21d-aggregated.pcr ##########
### Given cycle = 365 and non-annual ###
rm(list = ls())
load('Rdata/Overall.ILIplus.21-day.aggregated.pcr.Rdata')
# source('03. Fit iliplus with step function/Multiple_stepfunc_v5.R')
source('functions/Multiple_stepfunc_v5.R')
source('functions/plotting.R')

fits_iliplus_2N_365 = fit_multi_stepfunc(data = ili_plus_21d_df$ili_plus_smth_7d,
                                         N = 2,
                                         cycle = 365,
                                         bps_lower_limit = 55,
                                         bps_upper_limit = 300,
                                         vals_lower_limit = 0.1,
                                         vals_upper_limit = 0.5,
                                         sigma_lower_limit = 0.166,
                                         sigma_upper_limit = 0.166,
                                         iters = 100)

save(fits_iliplus_2N_365,file = 'Rdata/2_stepfunc_cyc365_iliplus_20230125.Rdata')

load('Rdata/2_stepfunc_cyc365_iliplus_holidayNA_20230125.Rdata')

best_fit_ilip_2N_365 = select_best_fit(all_fits = fits_iliplus_2N_365,
                                       N = 2,
                                       cycle = 365,
                                       iters = 100)

best_fit_ilip_2N_365
# bp1 = 54 bp2 = 324 val1 = 0.291 val2 = 0.185 sigma = 0.166
# AIC = -2020

multi_stepfunc_aic(data = ili_plus_21d_df$ili_plus_smth_7d,
                   cycle = 365,
                   c(bp1 = 55,bp2 = 314, val1 = 0.287, val2 = 0.182,sigma = 0.166))
## AIC = -2024

# start day in a yearly scale: 
# format(as.Date('2012-05-23'),'%j')
# 144 + 55 = 199
# 144 + 314 - 365 = 103

plot_data_with_stepfunc_by_year_fixed_cycle(data_df = ili_plus_21d_df,
                                            data_name = 'ili_plus_smth_7d',
                                            time_name = 'date',
                                            # stepfunc_params = best_fit_ilip_2N_365$all_results,
                                            stepfunc_params = c(bp1 = 55,bp2 = 314, val1 = 0.287, val2 = 0.182,sigma = 0.166),
                                            ymin = 0, 
                                            ymax = 1,
                                            cycle = 365,
                                            Nrow = 4,
                                            Ncol = 2)

# High period is during Beginning of March to mid-July, missed peak in 2017, 2018.

############# Cycle = 330 ##############
# fits_iliplus_2N_330 = fit_multi_stepfunc(data = ili_plus_21d_df$ili_plus_smth_7d,
#                                          N = 2,
#                                          cycle = 330,
#                                          bps_lower_limit = 1,
#                                          bps_upper_limit = 330,
#                                          vals_lower_limit = 0.1,
#                                          vals_upper_limit = 0.3,
#                                          sigma_lower_limit = 0.16,
#                                          sigma_upper_limit = 0.17,
#                                          iters = 150)
# save(fits_iliplus_2N_330,file = 'Rdata/2_stepfunc_cyc330_iliplus.Rdata')

load('Rdata/2_stepfunc_cyc330_iliplus_holidayNA.Rdata')
best_fit_ilip_2N_330 = select_best_fit(all_fits = fits_iliplus_2N_330,
                                       N = 2,
                                       cycle = 330,
                                       iters = 150)
best_fit_ilip_2N_330
# bp1 = 1 bp2 = 85 val1 = 0.166 val2 = 0.341 sigma = 0.154 
# AIC = -2414

plot_data_with_stepfunc_by_year_fixed_cycle(data_df = ili_plus_21d_df,
                                            data_name = 'ili_plus_smth_7d',
                                            time_name = 'date',
                                            stepfunc_params = best_fit_ilip_2N_330$all_results,
                                            cycle = 330,
                                            ymin = 0,ymax = 0.8,
                                            Nrow = 4,Ncol = 2)

######## 385-day cycle #########
# fits_iliplus_2N_385 = fit_multi_stepfunc(data = ili_plus_21d_df$ili_plus_smth_7d,
#                                          N = 2,
#                                          cycle = 385,
#                                          bps_lower_limit = 1,
#                                          bps_upper_limit = 385,
#                                          vals_lower_limit = 0.1,
#                                          vals_upper_limit = 0.3,
#                                          sigma_lower_limit = 0.16,
#                                          sigma_upper_limit = 0.17,
#                                          iters = 150)
# save(fits_iliplus_2N_385,file = 'Rdata/2_stepfunc_cyc385_iliplus.Rdata')

load('Rdata/2_stepfunc_cyc385_iliplus_holidayNA.Rdata')
best_fit_ilip_2N_385 = select_best_fit(all_fits = fits_iliplus_2N_385,
                                       N = 2,
                                       cycle = 385,
                                       iters = 150)
best_fit_ilip_2N_385
# bp1 = 47 bp2 = 277 val1 = 0.301 val2 = 0.146 sigma = 0.154
# AIC = -2405

plot_data_with_stepfunc_by_year_fixed_cycle(data_df = ili_plus_21d_df,
                                            data_name = 'ili_plus_smth_7d',
                                            time_name = 'date',
                                            stepfunc_params = best_fit_ilip_2N_385$all_results,
                                            cycle = 385,
                                            ymin = 0,ymax = 0.8,
                                            Nrow = 4,Ncol = 2)
