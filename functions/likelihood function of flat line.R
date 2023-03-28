
flat_line_ll = function(data,params) {
    
    val = params['val']
    sigma = params['sigma']
    
    sim_ts = rep(val,length(data))
    #sim_ts = rnorm(length(data),mean = val, sd = sigma)
    
    log.lik = sum(dnorm(data,mean = sim_ts, sd = sigma, log = T),na.rm = T)
    
    # number of parameters 
    # k = length(val) + length(sigma)
    # 
    # # AIC 
    # aic = (2*k) - (2*log.lik)
    
    return(log.lik)
}
