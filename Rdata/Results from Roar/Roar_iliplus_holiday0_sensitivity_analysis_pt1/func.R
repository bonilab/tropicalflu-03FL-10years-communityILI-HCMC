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
                             sigma = NULL,
                             seed = NULL){
    
    if(!is.null(seed)) {
        
        set.seed(seed)
    }
    
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

multi_stepfunc_aic = function(data, cycle, ...){
    
    ellipsis_list = as.list(...)
    
    bps = list_to_vector(ellipsis_list, prefix = 'bp')
    vals = list_to_vector(ellipsis_list, prefix = 'val')
    sigma = list_to_vector(ellipsis_list, prefix = 'sigma')
    
    # cat('\n bps: ', bps, '\n vals: ', vals, '\n sigma: ', sigma, '\n\n')
    
    #### bps cannot be the start or the end of a cycle 
    ## bps must be sorted
    ## every bps is at least 15 days distant from its neighbor values
    
    if(bps[1] <= 0 || bps[length(bps)] >= cycle || prod(bps == sort(bps)) != 1 ||
       prod(sapply(1:(length(bps)-1),function(x){bps[x+1] - bps[x] >=15})) != 1 ||
       bps[1] + cycle - bps[length(bps)] <15){
        
        log.lik = -1e32
        # sse = 1e32
        
    } else {
        
        #### Simulate time series
        sim_ts = sapply(1:length(data), multi_stepfunc,
                        cycle = cycle, c(bps,vals))
        
        log.lik = sum(dnorm(data, mean = sim_ts, sd = sigma, log = TRUE),
                      na.rm = T)
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

#### Main function ####
fit_multi_stepfunc = function(data,N,cycle,
                              bps_lower_limit,
                              bps_upper_limit,
                              vals_lower_limit,
                              vals_upper_limit,
                              sigma_lower_limit, 
                              sigma_upper_limit, 
                              seed = seed){
    
    #### Optimize from random start parameters
    # browser()
    
    start_value = params_initialize(N = N, 
                                    cycle = cycle, 
                                    bps_lower_limit = bps_lower_limit,
                                    bps_upper_limit = bps_upper_limit,
                                    vals_lower_limit = vals_lower_limit,
                                    vals_upper_limit = vals_upper_limit, 
                                    sigma_lower = sigma_lower_limit, 
                                    sigma_upper = sigma_upper_limit,
                                    seed = seed)
    #### Use NM in optim 
    fit = optim(par = start_value,
                multi_stepfunc_aic,
                data = data,
                cycle = cycle,
                method = "Nelder-Mead",
                control = list(maxit = 25000))
    
    return(fit)
    
}

select_best_fit = function(all_fits, N, cycle, iters){
    
    all_converge = sapply(1:length(all_fits),function(x){all_fits[[x]]$convergence})
    
    converged_fits = all_fits[which(all_converge == 0)]
    
    if(length(converged_fits) == 0){
        best_fit = list()
    } else {
        
        all_aic = sapply(1:length(converged_fits),function(x){converged_fits[[x]]$value})
        all_params = lapply(1:length(converged_fits),function(x){converged_fits[[x]]$par})
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

