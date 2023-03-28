###### All Functions #####

params_initialize = function(N,
                             cycle_lower_limit,
                             cycle_upper_limit,
                             bps_lower_limit,
                             bps_upper_limit,
                             vals_lower_limit,
                             vals_upper_limit,
                             sigma_lower_limit, 
                             sigma_upper_limit, 
                             random = T,
                             whole_cycle_search = T,
                             bps = NULL,
                             vals = NULL,
                             sigma = NULL){
    
    if(bps_lower_limit <= 0){
        
        stop("bps_lower_limit is out of range")
    }
    
    if(random == T){
        
        ### start value is randomly drawn from a given range
        cycle = sample(c(cycle_lower_limit,cycle_upper_limit),1)
        
        ### If search through the entire cycle is needed, the upper limit of bp is the sampled cycle ##
        if(whole_cycle_search == T){
            
            bps = sort(sample(c(bps_lower_limit:cycle),N))
            
        } else if(whole_cycle_search == F){
            
            if(bps_upper_limit >cycle_lower_limit){
                stop("bps_upper_limit is out of range")
            }
            bps = sort(sample(c(bps_lower_limit:bps_upper_limit),N))
            
        } else {
            
            stop("The input of whole_cycle_search is not permitted")
        }
        
        vals = runif(N, min = vals_lower_limit, max = vals_upper_limit)
        sigma = runif(1, min = sigma_lower_limit, max = sigma_upper_limit)
        
    } else if(random == F) {
        
        #### bps must be in a sorted order and consistent with vals ####
        
        bps = bps
        vals = vals
        sigma = sigma
        cycle = cycle
        
    } else {
        stop(paste("The input of random is not permitted"))
    }
    
    ##### Name the parameters, each parameter must be named in order ####
    names(bps) = paste0("bp",1:length(bps))
    names(vals) = paste0("val",1:length(vals))
    names(sigma) = 'sigma'
    names(cycle) = "cycle"
    
    
    return(c(bps, vals, sigma,cycle))
}


list_to_vector <- function(l, prefix) {
    return(unlist(l[grep(prefix,names(l))]))
}


multi_stepfunc = function(t,...){
    
    #### Capture varying length of bps and vals
    
    # browser()
    ellipsis_list = as.list(...)
    
    # print(ellipsis_list)
    
    bps = list_to_vector(ellipsis_list, prefix = 'bp')
    vals = list_to_vector(ellipsis_list, prefix = 'val')
    cycle = list_to_vector(ellipsis_list,prefix = "cycle")
    
    t = t %% cycle
    
    ### The first step starts with the first day of the cycle
    ### The last step ends with the last day in the cycle
    
    y = vals[1]
    
    for (i in 1:(length(bps)-1)) {
        y = y + 
            ((vals[i+1] - vals[1]) * (t > bps[i] & t <= bps[i+1]))
    }
    
    return(y)
    
}


##### Function that calculate AIC ###

multi_stepfunc_aic = function(data,bps_diff, ...){
    
    ellipsis_list = as.list(...)
    
    bps = list_to_vector(ellipsis_list, prefix = 'bp')
    vals = list_to_vector(ellipsis_list, prefix = 'val')
    sigma = list_to_vector(ellipsis_list, prefix = 'sigma')
    cycle = list_to_vector(ellipsis_list,prefix = "cycle")
    
    # cat('\n bps: ', bps, '\n vals: ', vals, '\n sigma: ', sigma, '\n\n')
    
    #### bps cannot be the start or the end of a cycle 
    ## bps must be sorted
    ## every bps is at least 15 days distant from its neighbor values
    
    #### Simulate time series
    # sim_ts = sapply(1:length(data), multi_stepfunc,c(bps,vals,cycle))
    # 
    # log.lik = sum(dnorm(data, mean = sim_ts, sd = sigma, log = TRUE),na.rm = T)
    # sse = sum((sim_ts - data)^2)
    

    if(bps[1] <= 0 || bps[length(bps)] >= cycle || prod(bps == sort(bps)) != 1 ||
       prod(sapply(1:(length(bps)-1),function(x){bps[x+1] - bps[x] >=bps_diff})) != 1 ||
       bps[1] + cycle - bps[length(bps)] <bps_diff ||
       sigma <= 0){

        log.lik = -1e32
        # sse = 1e32

    } else {

        #### Simulate time series
        sim_ts = sapply(1:length(data), multi_stepfunc,c(bps,vals,cycle))

        log.lik = sum(dnorm(data, mean = sim_ts, sd = sigma, log = TRUE),na.rm = T)
        # sse = sum((sim_ts - data)^2)

    }

    k = (length(sigma) + length(bps) + length(vals)) 
    
    # cat('\n k = ', k, ', log.lik = ', log.lik, '\n\n')
    
    aic = (2 * k) - (2 * log.lik)
    
    return(aic)
    
}

#start = params_initialize(N = 2, cycle = 365, bps_lower_limit = 5, bps_upper_limit = 300, 
#	vals_lower_limit = 0.2, vals_upper_limit = 2, sigma_lower = 0.5, sigma_upper = 2)
# zz = multi_stepfunc_aic(data = dat, cycle = 365, z2$best_par)

library(pbmcapply)
#### Main function ####
fit_multi_stepfunc = function(data,N,
                              bps_diff,
                              cycle_lower_limit,
                              cycle_upper_limit,
                              bps_lower_limit,
                              bps_upper_limit,
                              vals_lower_limit,
                              vals_upper_limit,
                              sigma_lower_limit, 
                              sigma_upper_limit, 
                              whole_cycle_search = T,
                              iters = 8){
    
    #### Optimize from random start parameters
    # browser()
    all_fits = list()
    
    all_fits = pbmclapply(1:iters, function(x){
        
        start_value = params_initialize(N = N, 
                                        cycle_lower_limit = cycle_lower_limit,
                                        cycle_upper_limit = cycle_upper_limit,
                                        bps_lower_limit = bps_lower_limit,
                                        bps_upper_limit = bps_upper_limit,
                                        vals_lower_limit = vals_lower_limit,
                                        vals_upper_limit = vals_upper_limit, 
                                        sigma_lower_limit = sigma_lower_limit, 
                                        sigma_upper_limit = sigma_upper_limit,
                                        whole_cycle_search = whole_cycle_search)
        #### Use NM in optim 
        fit = optim(par = start_value,
                    multi_stepfunc_aic,
                    data = data,
                    bps_diff = bps_diff,
                    method = "Nelder-Mead",
                    control = list(maxit = 30000))
        
        return(list(start_value = start_value,fit = fit))
        
    },mc.cores = 8)
    
}

select_best_fit = function(all_fits, N,iters){
    
    all_converge = sapply(1:length(all_fits),function(x){all_fits[[x]]$fit$convergence})
    
    converged_fits = all_fits[which(all_converge == 0)]
    
    if(length(converged_fits) == 0){
        best_fit = list()
    } else {
        
        all_aic = sapply(1:length(converged_fits),function(x){converged_fits[[x]]$fit$value})
        all_params = lapply(1:length(converged_fits),function(x){converged_fits[[x]]$fit$par})
        all_params = do.call(rbind,all_params)
        
        all_results = cbind(all_params,all_aic)
        
        best_fit = list(all_results = all_results[which.min(all_results[,"all_aic"]),],
                        N = N,
                        converged_rate = length(converged_fits)/iters)
    }
    
    return(best_fit)
}


##### Test multi-stepfunc-cycle-as-param ###
# rm(list = ls())
# source("03. Fit iliplus with step function/Multiple_stepfunc_cycle_as_param.R")
# 
# ### Use simulated data
# sim_data = rep(c(rep(10,20),c(10:60),rep(60,30)),6)
# plot(sim_data)
# bps = c(bp1 = 10,bp2 = 30)
# vals = c(val1 = 50,val2 = 80)
# cycle = 80
# sigma = 1
# 
# multi_stepfunc_aic(data = sim_data,
#                    c(bps,vals,cycle = cycle,sigma = sigma))
# 
# fit = optim(par = c(bps,vals,cycle = cycle,sigma = sigma),
#             multi_stepfunc_aic,
#             data = sim_data,
#             method = "Nelder-Mead",
#             control = list(maxit = 10000))
# fit
# 
# plot(sim_data,type = "l")      
# lines(sapply(1:length(sim_data),multi_stepfunc,fit$par),col = "red")
# 
# ### Use zeta ###
# ### Grid search on bp2 and cycle ###
# load("../Rdata/zeta score 1 of ILI percent.RData")
# bp1_can = seq(40,60,5)
# bp2_can = seq(130,150,5)
# cycle_can = seq(200,400,5)
# 
# VAL1 = 1
# VAL2 = 0.9
# 
# SIGMA = 0.09
# 
# args_comb = expand.grid(bp1_can,bp2_can,cycle_can)
# colnames(args_comb) = c("bp1","bp2","cycle")
# 
# library(pbmcapply)
# aic = pbmclapply(1:nrow(args_comb),function(x){
#     multi_stepfunc_aic(data = zeta_perc_avg_df$smoothed_zeta_score,
#                        c(bp1 = args_comb[x,"bp1"],
#                          bp2 = args_comb[x,"bp2"],
#                          val1 = VAL1,
#                          val2 = VAL2,
#                          sigma = SIGMA,
#                          cycle = args_comb[x,"cycle"]))},
#     mc.cores = 5)
# 
# args_comb$aic = do.call(rbind,aic)
# 
# args_comb[which.min(args_comb$aic),]
# 
# plot_args_comb = args_comb[args_comb$bp1 == 55,]
# 
# library(ggplot2)
# ggplot(plot_args_comb) +
#     geom_raster(aes(x = cycle,y = bp2,fill = aic)) +
#     scale_fill_viridis_c(option = "A")
# 
# plot(x = cycle_can,y = plot_args_comb[plot_args_comb$bp2 == 130,"aic"],type = "l")
# 
# fit_zeta = optim(par = c(bp1 = 55,bp2 = 130,val1 = 1,val2 = 0.9,sigma = 0.09,cycle = 350),
#                  multi_stepfunc_aic,
#                  data = zeta_perc_avg_df$smoothed_zeta_score,
#                  method = "Nelder-Mead",
#                  control = list(maxit = 10000))
# fit_zeta

# Sys.time()
# all_fits = fit_multi_stepfunc(data = zeta_perc_avg_df$smoothed_zeta_score,
#                               N = 2,
#                               cycle_lower_limit = 190,
#                               cycle_upper_limit = 380,
#                               bps_lower_limit = 1,
#                               bps_upper_limit = 190,
#                               vals_lower_limit = 0.9,
#                               vals_upper_limit = 1,
#                               sigma_lower_limit = 0.05,
#                               sigma_upper_limit = 0.15,
#                               whole_cycle_search = T,
#                               iters = 200)
# Sys.time()
# save(all_fits,file = "../Rdata/200_fits_zeta_2_step_cycle_unfixed_whole_cycle_search.Rdata")
# 
# best_fit = select_best_fit(all_fits = all_fits,N = 2,iters = 200)
# 
# best_fit

# all_converge = sapply(1:length(all_fits),function(x){all_fits[[x]]$fit$convergence})
# 
# converged_fits = all_fits[which(all_converge == 0)]
# all_aic = sapply(1:length(converged_fits),function(x){converged_fits[[x]]$fit$value})
# all_params = lapply(1:length(converged_fits),function(x){converged_fits[[x]]$fit$par})
# all_params = do.call(rbind,all_params)
# 
# all_results = cbind(all_params,all_aic)
# best_results = all_results[all_results[,"all_aic"]<quantile(all_results[,"all_aic"],0.1),]
# 
# hist(best_results[,"cycle"])
# hist(best_results[,"all_aic"])
