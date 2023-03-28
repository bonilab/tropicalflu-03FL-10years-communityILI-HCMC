
### Simulate new cases via SIRS model ----

case_simulate = function(params) {
    
    source("SIRS model.R")
    
    ### Input parameters 
    R_num = params[["R_num"]]
    beta = params[["beta"]]
    reporting_rate = params[["reporting_rate"]]
    daynum = params[["daynum"]]
    
    ### Set burn-in and total time fixed, same length as obs
    burn_in = 1460
    total = 1460 + daynum
    
    ### Solve ODE
    out = sirs_multiR_solve(num_Rs = R_num, beta = beta, burn_in = burn_in, total = total)
    out = data.frame(out)
    
    ### Get cumulative incidence
    J = out$J
    
    ### Get new cases
    new_cases = round(J[2:length(J)] -J[1:length(J)-1])
    
    ### Add reporting stochasticity, following binomial
    reported_cases = sapply(new_cases, function(x){rbinom(n = 1, size = x, prob = reporting_rate)})
    
    return(reported_cases)
    
}


get_nllk_pois_case = function(params,obs,daynum) {
    source("SIRS model.R")
    
    ### Input parameters 
    R_num = params[["R_num"]]
    
    if(R_num %%1 !=0) {
        
        nllk = Inf
    } else {
        
        ### Constrain R_num to be an integer
        beta = params[["beta"]]
        reporting_rate = params[["reporting_rate"]]
        daynum = daynum
        
        ### Set burn-in and total fixed, same length as obs
        burn_in = 1460
        total = 1460 + daynum
        
        ### Solve ODE
        out = sirs_multiR_solve(num_Rs = R_num, beta = beta, burn_in = burn_in, total = total)
        out = data.frame(out)
        
        ### Get cumulative incidence
        J = out$J
        
        ### Get new cases
        new_cases = round(J[2:length(J)] -J[1:length(J)-1])
        
        ### Add reporting stochasticity, following binomial
        reported_cases = sapply(new_cases, function(x){rbinom(n = 1, size = x, prob = reporting_rate)})
        
        ### Confront simulated reported cases with observed incidence converted from iliplus
        
        ### poisson distribution
        llk = mapply(function(x,lambda){
            
            llk = dpois(x,lambda = lambda,log = T)
            return(llk)
        },x = obs,lambda = reported_cases)
        
        nllk = -sum(llk)
    }
    
    return(nllk)
    
}

#### Add dispersion parameter k in negative binomial regression

get_nllk_negbinom = function(params,obs,daynum) {
    
    ### Input parameters 
    R_num = params[["R_num"]]
    beta = params[["beta"]]
    reporting_rate = params[["reporting_rate"]]
    k = params[["k"]]
    daynum = daynum
    
    if(beta >1 |reporting_rate >1 | R_num %% 1 !=0) {
        
        nllk = Inf
    } else {
        
        ### Solve ODE
        
        ### Set burn-in and total time fixed, same length as obs
        burn_in = 1460
        total = 1460 + daynum
        
        ### Solve ODE
        source("SIRS model.R")
        out = sirs_multiR_solve(num_Rs = R_num, beta = beta, burn_in = burn_in, total = total)
        out = data.frame(out)
        
        ### Get cumulative incidence
        J = out$J
        
        ### Get new cases
        new_cases = round(J[2:length(J)] -J[1:length(J)-1])
        
        ### Add reporting stochasticity, following binomial
        reported_cases = sapply(new_cases, function(x){rbinom(n = 1, size = x, prob = reporting_rate)})
        
        ### Confront simulated reported cases with observed incidence converted from iliplus
        
        ### negative binomial distribution, also need to estimate dispersion parameter?
        llk = mapply(function(x,size,mu){
            
            llk = dnbinom(x = x,size = size,mu = mu,log = T)
            return(llk)
        },x = obs,size = k,mu = reported_cases)
        
        nllk = -sum(llk)
    }
    
    return(nllk)
}

### Input parameters seperately:
nllk_negbinom = function(R_num,beta,reporting_rate,k,obs,daynum) {
    
    ### Input parameters 
    R_num = R_num
    beta = beta
    reporting_rate = reporting_rate
    k = k
    daynum = daynum
    
    if(beta >1 |reporting_rate >1 | R_num %% 1 !=0) {
        
        nllk = Inf
    } else {
        
        ### Solve ODE
        
        ### Set burn-in and total time fixed, same length as obs
        burn_in = 1460
        total = 1460 + daynum
        
        ### Solve ODE
        source("SIRS model.R")
        out = sirs_multiR_solve(num_Rs = R_num, beta = beta, burn_in = burn_in, total = total)
        out = data.frame(out)
        
        ### Get cumulative incidence
        J = out$J
        
        ### Get new cases
        new_cases = round(J[2:length(J)] -J[1:length(J)-1])
        
        ### Add reporting stochasticity, following binomial
        reported_cases = sapply(new_cases, function(x){rbinom(n = 1, size = x, prob = reporting_rate)})
        
        ### Confront simulated reported cases with observed incidence converted from iliplus
        
        ### negative binomial distribution, also need to estimate dispersion parameter?
        llk = mapply(function(x,size,mu){
            
            llk = dnbinom(x = x,size = size,mu = mu,log = T)
            return(llk)
        },x = obs,size = k,mu = reported_cases)
        
        nllk = -sum(llk)
    }
    
    return(nllk)
}



