#### Estimate 2-stepfunc of 7d-smoothed zeta score ####
rm(list = ls())
load("Rdata/7d.smth.zeta.ilipercent.Rdata")
source("functions/Multiple_stepfunc_v5.R")
source("functions/Detrend and Smooth.R")
source('functions/plotting.R')

####### Try 1000 fits, to calculate how many repeats are needed to get global minimum, taking account of NM degeneracy issue

### sigma needs to change to be close to the SD of 7d-smoothed zeta

#### Cycle = 365 ###
# Based on grid search, best start value: bp1 = 55,bp2 = 135, val1 = 1.02, val2 = 0.92, sigma = 0.122

# fits_2N_365_7d_zeta = fit_multi_stepfunc(data = zeta_perc_avg_df$smoothed_zeta_score,
#                                          N = 2,
#                                          cycle = 365,
#                                          bps_lower_limit = 55,
#                                          bps_upper_limit = 135,
#                                          vals_lower_limit = 0.92,
#                                          vals_upper_limit = 1.02,
#                                          sigma_lower_limit = 0.11,
#                                          sigma_upper_limit = 0.13,
#                                          iters = 100)
# 
# 
# save(fits_2N_365_7d_zeta,file = "Rdata/1000 fits of 2-stepfunc for 7d_smth_zeta_cycle_365.Rdata")

load("Rdata/1000 fits of 2-stepfunc for 7d_smth_zeta_cycle_365.Rdata")
best_fit_365 = select_best_fit(all_fits = fits_2N_365_7d_zeta,
                               N = 2,
                               cycle = 365,
                               iters = 1000)

best_fit_365
# params: bp1 = 58 bp2 = 136 val1 = 1.015 val2 = 0.915 sigma = 0.122 
# AIC = -5016


#### plot with zeta ###
plot_data_with_stepfunc_by_year_fixed_cycle(data_df = zeta_perc_avg_df,
                                            data_name = 'smoothed_zeta_score',
                                            time_name = 'Date',
                                            stepfunc_params = best_fit_365$all_results,
                                            cycle = 365,
                                            Nrow = 5,
                                            Ncol = 2)

### Cycle is non-annual ###
# Based on grid search, best start value: bp1 = 70, bp2 = 140, val1 = 1.02, val2 = 0.94, sigma = 0.122, cycle = 210

# fits_2N_210_7d_zeta = fit_multi_stepfunc(data = zeta_perc_avg_df$smoothed_zeta_score,
#                                          N = 2,
#                                          cycle = 210,
#                                          bps_lower_limit = 70,
#                                          bps_upper_limit = 140,
#                                          vals_lower_limit = 0.94,
#                                          vals_upper_limit = 1.04,
#                                          sigma_lower_limit = 0.11,
#                                          sigma_upper_limit = 0.13,
#                                          iters = 1000)
# save(fits_2N_210_7d_zeta,file = "Rdata/1000 fits of 2-stepfunc for 7d_smth_zeta_cycle_210.Rdata")
load("Rdata/1000 fits of 2-stepfunc for 7d_smth_zeta_cycle_210.Rdata")
    
best_fit_210 = select_best_fit(all_fits = fits_2N_210_7d_zeta,
                               N = 2,
                               cycle = 210,
                               iters = 1000)
best_fit_210
# params: bp1 = 69 bp2 = 175 val1 = 1.027 val2 = 0.954 sigma = 0.122 
# AIC = -4931

##### results: annual cycle fits better than non-annual cycle #######


### Plot with 7d.smth.zeta ###

#### plot with zeta ###
plot_data_with_stepfunc_by_year_fixed_cycle(data_df = zeta_perc_avg_df,
                                            data_name = 'smoothed_zeta_score',
                                            time_name = 'Date',
                                            stepfunc_params = best_fit_210$all_results,
                                            cycle = 210,
                                            Nrow = 5,
                                            Ncol = 2)















