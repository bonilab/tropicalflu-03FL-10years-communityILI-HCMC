#### Profile likelihood is applied on 2-step function only #####################
library(pbmcapply)

##### 2-step function simulation ####
multi_stepfunc_2N = function(t, cycle, bp1, bp2, val1, val2){
    
    #### Capture varying length of bps and vals
    
    bps = c(bp1,bp2)
    vals = c(val1,val2)
    
    t = t %% cycle
    
    ### The first step starts with the first day of the cycle
    ### The last step ends with the last day in the cycle
    
    y = val1
    
    y = y + (val2 - val1) * (t > bp1 & t <= bp2)
    
    return(y)
    
}

##### profile AIC(Likelihood) function on each parameter #####
profile_bp1_2_stepfunc = function(data, cycle, params, bp1, sd = 0.1) {
    
    # cat('\n bps: ', bps, '\n vals: ', vals, '\n sigma: ', sigma, '\n\n')
    
    #### bps cannot be the start or the end of a cycle 
    ## bps must be sorted
    ## every bps is at least 15 days distant from its neighbor values
    
    bp2 = params['bp2']
    val1 = params['val1']
    val2 = params['val2']
    # sigma = params['sigma']
    
    bps = c(bp1,bp2)
    vals = c(val1,val2)
    
    if(bp1 <= 0 || bp2 >= cycle || bp1 > bp2 ||
       bp2 - bp1 < 15 ||
       bp1 + cycle - bp2 <15){
        
        log.lik = -1e32
        # sse = 1e32
        
    } else {
        
        #### Simulate time series
        sim_ts = sapply(1:length(data), multi_stepfunc_2N,
                        cycle = cycle, bp1 = bp1,
                        bp2 = bp2,val1 = val1,
                        val2 = val2)
        
        log.lik = sum(dnorm(data, mean = sim_ts, sd = sd, log = TRUE),na.rm = T)
        # sse = sum((sim_ts - data)^2)
        
    }
    
    # k = (length(sigma) + length(bps) + length(vals)) 
    
    # cat('\n k = ', k, ', log.lik = ', log.lik, '\n\n')
    
    # aic = (2 * k) - (2 * log.lik)
    
    return(log.lik)
    
}

profile_bp2_2_stepfunc = function(data, cycle, params, bp2, sd = 0.1) {
    
    # cat('\n bps: ', bps, '\n vals: ', vals, '\n sigma: ', sigma, '\n\n')
    
    #### bps cannot be the start or the end of a cycle 
    ## bps must be sorted
    ## every bps is at least 15 days distant from its neighbor values
    
    bp1 = params['bp1']
    val1 = params['val1']
    val2 = params['val2']
    # sigma = params['sigma']
    
    bps = c(bp1,bp2)
    vals = c(val1,val2)
    
    if(bp1 <= 0 || bp2 >= cycle || bp1 > bp2 ||
       bp2 - bp1 < 15 ||
       bp1 + cycle - bp2 <15){
        
        log.lik = -1e32
        # sse = 1e32
        
    } else {
        
        #### Simulate time series
        sim_ts = sapply(1:length(data), multi_stepfunc_2N,
                        cycle = cycle, bp1 = bp1,
                        bp2 = bp2,val1 = val1,
                        val2 = val2)
        
        log.lik = sum(dnorm(data, mean = sim_ts, sd = sd, log = TRUE),na.rm = T)
        # sse = sum((sim_ts - data)^2)
        
    }
    
    # k = (length(sigma) + length(bps) + length(vals)) 
    
    # cat('\n k = ', k, ', log.lik = ', log.lik, '\n\n')
    
    # aic = (2 * k) - (2 * log.lik)
    
    return(log.lik)
    
}

profile_val1_2_stepfunc = function(data, cycle, params, val1, sd = 0.1) {
    
    # cat('\n bps: ', bps, '\n vals: ', vals, '\n sigma: ', sigma, '\n\n')
    
    #### bps cannot be the start or the end of a cycle 
    ## bps must be sorted
    ## every bps is at least 15 days distant from its neighbor values
    
    bp1 = params['bp1']
    bp2 = params['bp2']
    val2 = params['val2']
    
    bps = c(bp1,bp2)
    vals = c(val1,val2)
    
    if(bp1 <= 0 || bp2 >= cycle || bp1 > bp2 ||
       bp2 - bp1 < 15 ||
       bp1 + cycle - bp2 <15){
        
        log.lik = -1e32
        # sse = 1e32
        
    } else {
        
        #### Simulate time series
        sim_ts = sapply(1:length(data), multi_stepfunc_2N,
                        cycle = cycle, bp1 = bp1,
                        bp2 = bp2,val1 = val1,
                        val2 = val2)
        
        log.lik = sum(dnorm(data, mean = sim_ts, sd = sd, log = TRUE),na.rm = T)
        # sse = sum((sim_ts - data)^2)
        
    }
    
    k = (length(sigma) + length(bps) + length(vals)) 
    
    # cat('\n k = ', k, ', log.lik = ', log.lik, '\n\n')
    
    # aic = (2 * k) - (2 * log.lik)
    
    return(log.lik)
    
}

profile_val2_2_stepfunc = function(data, cycle, params, val2, sd = 0.1) {
    
    # cat('\n bps: ', bps, '\n vals: ', vals, '\n sigma: ', sigma, '\n\n')
    
    #### bps cannot be the start or the end of a cycle 
    ## bps must be sorted
    ## every bps is at least 15 days distant from its neighbor values
    
    bp1 = params['bp1']
    bp2 = params['bp2']
    val1 = params['val1']
    # sigma = params['sigma']
    
    bps = c(bp1,bp2)
    vals = c(val1,val2)
    
    if(bp1 <= 0 || bp2 >= cycle || bp1 > bp2 ||
       bp2 - bp1 < 15 ||
       bp1 + cycle - bp2 <15){
        
        log.lik = -1e32
        # sse = 1e32
        
    } else {
        
        #### Simulate time series
        sim_ts = sapply(1:length(data), multi_stepfunc_2N,
                        cycle = cycle, bp1 = bp1,
                        bp2 = bp2,val1 = val1,
                        val2 = val2)
        
        log.lik = sum(dnorm(data, mean = sim_ts, sd = sd, log = TRUE),na.rm = T)
        # sse = sum((sim_ts - data)^2)
        
    }
    
    # k = (length(sigma) + length(bps) + length(vals)) 
    
    # cat('\n k = ', k, ', log.lik = ', log.lik, '\n\n')
    
    # aic = (2 * k) - (2 * log.lik)
    
    return(log.lik)
    
}

select_best_fit_llk = function(all_fits, repeats){
    
    all_converge = sapply(1:length(all_fits),function(x){all_fits[[x]]$convergence})
    
    converged_fits = all_fits[which(all_converge == 0)]
    
    if(length(converged_fits) == 0){
        best_fit = list()
    } else {
        
        all_llk = sapply(1:length(converged_fits),function(x){converged_fits[[x]]$value})
        all_params = lapply(1:length(converged_fits),function(x){converged_fits[[x]]$par})
        all_params = do.call(rbind,all_params)
        
        all_results = cbind(all_params,all_llk)
        
        best_fit = list(all_results = all_params[which.max(all_results[,"all_llk"]),],
                        llk = max(all_llk),
                        converged_rate = length(converged_fits)/repeats)
    }
    
    return(best_fit)
}

#### wrapping function of all the profile parameter functions ######
profile_llk_2_stepfunc = function(data,cycle,
                                  param_to_profile,
                                  profile_param_value,
                                  bp1_lower_limit = NULL,
                                  bp1_upper_limit = NULL,
                                  bp2_lower_limit = NULL,
                                  bp2_upper_limit = NULL,
                                  val1_lower_limit = NULL,
                                  val1_upper_limit = NULL,
                                  val2_lower_limit = NULL,
                                  val2_upper_limit = NULL,
                                  sd = 0.1,
                                  repeats = 10) {
    if(param_to_profile == 'bp1') {
        
        all_fits = pbmclapply(1:repeats,function(x) {
            
            bp2_start = sample(bp2_lower_limit:bp2_upper_limit,1)
            val1_start = sample(val1_lower_limit:val1_upper_limit,1)
            val2_start = sample(val2_lower_limit:val2_upper_limit,1)
            # sigma_start = sample(sigma_lower_limit:sigma_upper_limit,1)
            
            start_params = c(bp2 = bp2_start,val1 = val1_start,val2 = val2_start)
            
            fit = optim(start_params,
                        profile_bp1_2_stepfunc,
                        data = data,
                        cycle = cycle,
                        bp1 = profile_param_value,
                        sd = sd, 
                        method = 'Nelder-Mead',
                        control = list(maxit = 10000,
                                       fnscale = -1))
            
            return(fit)
        },mc.cores = 5)
        
    } else if(param_to_profile == 'bp2') {
        
        all_fits = pbmclapply(1:repeats,function(x) {
            
            bp1_start = sample(bp1_lower_limit:bp1_upper_limit,1)
            val1_start = sample(val1_lower_limit:val1_upper_limit,1)
            val2_start = sample(val2_lower_limit:val2_upper_limit,1)
            # sigma_start = sample(sigma_lower_limit:sigma_upper_limit,1)
            
            start_params = c(bp1 = bp1_start,val1 = val1_start,val2 = val2_start)
            
            fit = optim(start_params,
                        profile_bp2_2_stepfunc,
                        data = data,
                        cycle = cycle,
                        bp2 = profile_param_value,
                        sd = sd,
                        method = 'Nelder-Mead',
                        control = list(maxit = 10000,
                                       fnscale = -1))
            
            return(fit)
        },mc.cores = 5) 
        
    } else if(param_to_profile == 'val1') {
        
        all_fits = pbmclapply(1:repeats,function(x) {
            
            bp1_start = sample(bp1_lower_limit:bp1_upper_limit,1)
            bp2_start = sample(bp2_lower_limit:bp2_upper_limit,1)
            val2_start = sample(val2_lower_limit:val2_upper_limit,1)
            # sigma_start = sample(sigma_lower_limit:sigma_upper_limit,1)
            
            start_params = c(bp1 = bp1_start,bp2 = bp2_start,val2 = val2_start)
            
            fit = optim(start_params,
                        profile_val1_2_stepfunc,
                        data = data,
                        cycle = cycle,
                        val1 = profile_param_value,
                        sd = sd,
                        method = 'Nelder-Mead',
                        control = list(maxit = 10000,
                                       fnscale = -1))
            
            return(fit)
        },mc.cores = 5) 
        
    } else if (param_to_profile == 'val2') {
        
        all_fits = pbmclapply(1:repeats,function(x) {
            
            bp1_start = sample(bp1_lower_limit:bp1_upper_limit,1)
            bp2_start = sample(bp2_lower_limit:bp2_upper_limit,1)
            val1_start = sample(val1_lower_limit:val1_upper_limit,1)
            # sigma_start = sample(sigma_lower_limit:sigma_upper_limit,1)
            
            start_params = c(bp1 = bp1_start,bp2 = bp2_start,val1 = val1_start)
            
            fit = optim(start_params,
                        profile_val2_2_stepfunc,
                        data = data,
                        cycle = cycle,
                        val2 = profile_param_value,
                        sd = sd,
                        method = 'Nelder-Mead',
                        control = list(maxit = 10000,
                                       fnscale = -1))
            
            return(fit)
        },mc.cores = 5)
        
    } else {
        
        stop('profile parameters must be one of the step function parameters!')
    } 
    
    best_fit = select_best_fit_llk(all_fits = all_fits,repeats = repeats)
    return(best_fit)
    
}
