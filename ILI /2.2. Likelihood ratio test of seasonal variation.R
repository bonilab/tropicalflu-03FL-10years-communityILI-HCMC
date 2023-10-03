###### 95% CI of likelihood interval based on Wilks' Theorem ####
### If the likelihood of a step function falls into the 95% CI 
## given flat line(null hypothesis: no seasonal variation)
# There is no seasonal variation
rm(list = ls())
source('functions/likelihood function of flat line.R')
source("functions/Multiple_stepfunc_v5.R")
load('Rdata/7d.smth.zeta.ilipercent.Rdata')
load('Rdata/1000 fits of 2-stepfunc for 7d_smth_zeta_cycle_365.Rdata')
load("Rdata/1000 fits of 2-stepfunc for 7d_smth_zeta_cycle_210.Rdata")

### The best 2-step function given ts is annual 
best_fit_365 = select_best_fit(all_fits = fits_2N_365_7d_zeta,
                               N = 2,
                               cycle = 365,
                               iters = 1000)

best_fit_365

### Get log-likelihood from AIC
best_ll_365 = (2*5 - best_fit_365$aic)/2 # 2513


### The best 2-step function given zeta is nonannual 
best_fit_210 = select_best_fit(all_fits = fits_2N_210_7d_zeta,
                               N = 2,
                               cycle = 210,
                               iters = 1000)
best_fit_210

#### Get log-likelihood from AIC 
best_ll_210 = (2*5 - best_fit_210$aic)/2 # 2471


### The best loglik fitting a straight line ###
flat_line_fit = optim(c(val = 1,sigma = 0.12),
                      flat_line_ll,
                      data = zeta_perc_avg_df$smoothed_zeta_score,
                      control = list(fnscale = -1))


flat_line_fit
# val1 = 0.994 sigma = 0.129 
# ll = 2310.726

# AIC
2*(2 - flat_line_fit$value) # -4617

### The -2*(log(L0/La)) ~ Chi-squared(difference of df)
# construct the 95% CI of null hypothesis
upper_ll = flat_line_fit$value + qchisq(p = 0.999,df = 3)/2
upper_ll # 2319


#### both annual and non-annual likelihood is higher than the critical region if there is no seasonal variation
# This indicates zeta shows siginificant seasonal variation
