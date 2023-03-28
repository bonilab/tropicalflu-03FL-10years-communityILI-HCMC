#### Likelihood ratio test for overall ILI+ ###
rm(list = ls())
source('functions/likelihood function of flat line.R')
source('functions/Multiple_stepfunc_v5.R')
load('Rdata/Overall.ILIplus.21-day.aggregated.pcr.RData')
load('Rdata/2_stepfunc_cyc365_iliplus_holiday0.Rdata')
load('Rdata/2_stepfunc_cyc385_iliplus_holiday0.Rdata')

best_fit_ilip_2N_365 = select_best_fit(all_fits = fits_iliplus_2N_365,
                                       N = 2,
                                       cycle = 365,
                                       iters = 100)

best_fit_ilip_2N_365

best_ll_365 = (2*5 - best_fit_ilip_2N_365$aic)/2 # 1045

best_fit_ilip_2N_385 = select_best_fit(all_fits = fits_iliplus_2N_385,
                                       N = 2,
                                       cycle = 385,
                                       iters = 100)

best_fit_ilip_2N_385

best_ll_385 = (2*5 - best_fit_ilip_2N_385$aic)/2 # 1271


#### log likelihood given a flat line ###
flat_line_fit = optim(c(val = mean(ili_plus_21d_df$ili_plus_smth_7d,na.rm = T),
                        sigma = 0.1),
                      flat_line_ll,
                      data = ili_plus_21d_df$ili_plus_smth_7d,
                      control = list(fnscale = -1))
flat_line_fit

### 95% CI from Wilks theorem ###
upper_ll = flat_line_fit$value + qchisq(p = 0.999,df = 3)/2
upper_ll # 974

## AIC of null model 
2*(2 - flat_line_fit$value) # - 1928
