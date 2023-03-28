###### All Functions #####

params_initialize = function(N,cycle,
                             bps_lower_limit,
                             bps_upper_limit,
                             vals_lower_limit,
                             vals_upper_limit,
                             sigma_lower, 
                             sigma_upper, 
                             random = T,
                             bps = NULL,
                             vals = NULL, 
                             sigma = NULL){
    
    if(random == T){
        
        ##### Initialize value #####
        
        
        
        if(bps_upper_limit > cycle | bps_lower_limit <=0){
            
            stop("bps out of range")
            
        } else {
            
            bps = sort(sample(c(bps_lower_limit:bps_upper_limit),N))
            
        }
        
        vals = runif(N, min = vals_lower_limit, max = vals_upper_limit)
        sigma = runif(1, min = sigma_lower, max = sigma_upper)
        
    } else if(random == F) {
        
        #### bps must be in a sorted order and consistent with vals ####
        
        bps = bps
        vals = vals
        sigma = sigma
        
    } else {
        stop(paste("The input of random is not permitted"))
    }
    
    ##### Name the parameters, each parameter must be named in order ####
    names(bps) = paste0("bp",1:length(bps))
    names(vals) = paste0("val",1:length(vals))
    names(sigma) = 'sigma'
    
    
    return(c(bps, vals, sigma))
}


list_to_vector <- function(l, prefix) {
    return(unlist(l[grep(prefix,names(l))]))
}


multi_stepfunc = function(t, cycle, ...){
    
    #### Capture varying length of bps and vals
    
    # browser()
    ellipsis_list = as.list(...)
    
    # print(ellipsis_list)
    
    bps = list_to_vector(ellipsis_list, prefix = 'bp')
    vals = list_to_vector(ellipsis_list, prefix = 'val')
    
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

multi_stepfunc_aic = function(data, cycle,...){
    
    ellipsis_list = as.list(...)
    
    bps = list_to_vector(ellipsis_list, prefix = 'bp')
    vals = list_to_vector(ellipsis_list, prefix = 'val')
    sigma = list_to_vector(ellipsis_list, prefix = 'sigma')
    
    # cat('\n bps: ', bps, '\n vals: ', vals, '\n sigma: ', sigma, '\n\n')
    
    #### bps cannot be the start or the end of a cycle 
    ## bps must be sorted
    ## every bps is at least 15 days distant from its neighbor values
    
    if(bps[1] <= 0 || bps[length(bps)] >= cycle || prod(bps == sort(bps)) != 1 ||
       prod(sapply(1:(length(bps)-1),function(x){bps[x+1] - bps[x] >=2})) != 1 ||
       bps[1] + cycle - bps[length(bps)] <2 ||
       sigma <= 0){
        
        log.lik = -1e32
        # sse = 1e32
        
    } else {
        
        #### Simulate time series
        sim_ts = sapply(1:length(data), multi_stepfunc,
                        cycle = cycle, c(bps,vals))
        
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
fit_multi_stepfunc = function(data,N,cycle,
                              bps_lower_limit,
                              bps_upper_limit,
                              vals_lower_limit,
                              vals_upper_limit,
                              sigma_lower_limit, 
                              sigma_upper_limit, 
                              iters = 10){
    
    #### Optimize from random start parameters
    # browser()
    all_fits = list()
    
    all_fits = pbmclapply(1:iters, function(x){
        
        start_value = params_initialize(N = N, 
                                        cycle = cycle, 
                                        bps_lower_limit = bps_lower_limit,
                                        bps_upper_limit = bps_upper_limit,
                                        vals_lower_limit = vals_lower_limit,
                                        vals_upper_limit = vals_upper_limit, 
                                        sigma_lower = sigma_lower_limit, 
                                        sigma_upper = sigma_upper_limit)
        #### Use NM in optim 
        fit = optim(par = start_value,
                    multi_stepfunc_aic,
                    data = data,
                    cycle = cycle,
                    method = "Nelder-Mead",
                    control = list(maxit = 20000))
        
        return(list(start_value = start_value,fit = fit))
        
    },mc.cores = 5)
    
}

select_best_fit = function(all_fits, N, cycle, iters){
    
    all_converge = sapply(1:length(all_fits),function(x){all_fits[[x]]$fit$convergence})
    
    converged_fits = all_fits[which(all_converge == 0)]
    
    if(length(converged_fits) == 0){
        best_fit = list()
    } else {
        
        all_aic = sapply(1:length(converged_fits),function(x){converged_fits[[x]]$fit$value})
        all_params = lapply(1:length(converged_fits),function(x){converged_fits[[x]]$fit$par})
        all_params = do.call(rbind,all_params)
        
        all_results = cbind(all_params,all_aic)
        
        best_fit = list(all_results = all_params[which.min(all_results[,"all_aic"]),],
                        aic = min(all_aic),
                        N = N,
                        cycle = cycle,
                        converged_rate = length(converged_fits)/iters)
    }
    
    return(best_fit)
}



#
# Simple example using code with simulated data
#

# cyc.length = 3250
# dat = sample(seq(0.5, 2, by = 0.1), 2*cyc.length, replace = TRUE)
# 
# z = fit_multi_stepfunc(data = dat, N = 3, cycle = cyc.length, 
# 	bps_lower_limit = 5, bps_upper_limit = 240, vals_lower_limit = 0.1, vals_upper_limit = 3, 
# 	sigma_lower_limit = 1, sigma_upper_limit = 2, iters = 5)
# z
# 	
# plot(dat, type = 'l')
# best_step = multi_stepfunc(c(1:length(dat)), cycle = cyc.length, z$best_par)
# points(best_step, type = 'l', col = 'red', lwd = 2)
# 
# 
# 
# # Example using codd wih VN ILI data
# 
# load('~/Desktop/Flu_Vietnam/Descr/ILI_pct_zeta.RData')
# 
# cyc.length = 365
# dat = zeta_perc_avg_df$smoothed_zeta_score
# 
# z = fit_multi_stepfunc(data = dat, N = 2, cycle = cyc.length, 
# 	bps_lower_limit = 5, bps_upper_limit = 200, vals_lower_limit = 0.1, vals_upper_limit = 3, 
# 	sigma_lower_limit = 1, sigma_upper_limit = 2, iters = 10)
# z
# 	
# plot(dat, type = 'l')
# best_step = multi_stepfunc(c(1:length(dat)), cycle = cyc.length, z$best_par)
# points(best_step, type = 'l', col = 'red', lwd = 2)
# 




# coarse grid search for initial values

# bp.vals = round(seq(5, 360, len = 10))
# step.vals = seq(0.6, 1.3, by = 0.1)
# sigma.vals = seq(0.05, 0.3, by = 0.05)
# 
# cons.pars = expand.grid('bp1' = bp.vals, 'bp2' = bp.vals, 'val1' = step.vals, 
# 	'val2' = step.vals, 'sigma' = sigma.vals)
# dim(cons.pars)
# cons.pars = cons.pars[which(cons.pars$bp1 < cons.pars$bp2), ]
# dim(cons.pars)
# # cons.pars = cons.pars[sort(sample(c(1:nrow(cons.pars)), 1000)), ]
# 
# grid.aics = list()
#     
# grid.aics = pbmclapply(1:nrow(cons.pars), function(x){
#     	
# 	multi_stepfunc_aic(data = dat, cycle = 365, cons.pars[x, 1:5])
# 
# }, mc.cores = 4)
# 
# cons.pars$aic = unlist(grid.aics)
# 
# 
# 
# # rerunning analysis with parameter values given lowest 10% of aic values
# 
# cons.pars = cons.pars[which(cons.pars$aic < quantile(cons.pars$aic, 0.1)), ]
# 
# z2 = fit_multi_stepfunc(data = dat, N = 2, cycle = cyc.length, 
# 	bps_lower_limit = max(floor(0.9*min(cons.pars$bp1)), 3), 
# 	bps_upper_limit = min(ceiling(1.1*max(cons.pars$bp2)), 360), 
# 	vals_lower_limit = max(0.9*min(c(cons.pars$val1, cons.pars$val2)), 0), 
# 	vals_upper_limit = 1.1*max(c(cons.pars$val1, cons.pars$val2)), 
# 	sigma_lower_limit = max(0.9*min(cons.pars$sigma), 0), 
# 	sigma_upper_limit = 1.1*max(cons.pars$sigma), 
# 	iters = 10)
# z2
# 	
# plot(dat, type = 'l')
# best_step = multi_stepfunc(c(1:length(dat)), cycle = cyc.length, z2$best_par)
# points(best_step, type = 'l', col = 'red', lwd = 2)
# 
# 
# 
