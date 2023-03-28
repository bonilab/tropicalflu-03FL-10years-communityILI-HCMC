library(pbmcapply)
library(deSolve)

setwd("~/Dropbox/Influenza and Respiratory Disease in Vietnam/_DATA/2019 12 31 - All SMS ILI Reports/data")
load("clean ILINUM.RData")

total_ilinum = rowSums(ILINUM,na.rm = T)
total_ilinum_df = data.frame(Date = seq(as.Date("2010-01-01"),as.Date("2019-12-31"),"1 day"),
                             ilinum = total_ilinum)

setwd("~/Dropbox/Influenza and Respiratory Disease in Vietnam/_DATA/Influenza PCR Data/Fuhan New Data Assembly October 2020 and later/Rdata")
load("ILI plus based on zeta score 1 of ILI percent.RData")

ili_plus_num_smoothed = data.frame(date = ili_plus_zeta_perc_df$Date,
                                   ili_plus = rep(NA,nrow(ili_plus_zeta_perc_df)))

for(i in 1:nrow(ili_plus_zeta_perc_df)){
    
    ind = which(total_ilinum_df$Date == ili_plus_zeta_perc_df$Date[i])
    ili_plus_num_smoothed$ili_plus[i] = ceiling(total_ilinum_df$ilinum[ind]*ili_plus_zeta_perc_df$smoothed_ili_plus[i])
}

# plot(ili_plus_num_smoothed$date,ili_plus_num_smoothed$ili_plus,type = "l")
# 
# summary(ili_plus_num_smoothed$ili_plus)
# hist(ili_plus_num_smoothed$ili_plus)



sin_plus3 = function(t,sigma,mean = 0, sd = 1){
    
    
    times = length(t) %/% 365
    
    if(times == 0){
        
        if(length(t) < (365*sigma/2)){
            
            phi = round(rlnorm(1,mean,sd))
            y = sin(2*pi*(t - phi)/(365*sigma))
            
        } else {
            
            phi = round(rlnorm(1,mean,sd))
            cycle_tmp = t - phi
            y = sin(2*pi*cycle_tmp/(365*sigma))
            high_period = seq(1+phi,365*sigma/2 + phi,1)
            y[-high_period] = 0
        }
    } else {
        cycle = c()
        total_high_period = c()
        
        for(i in 1:times){
            
            phi = round(rlnorm(1,mean,sd))
            cycle_tmp = c(1:365) - phi
            cycle = append(cycle,cycle_tmp)
            
            high_period = seq(1 +phi + 365*(i -1),phi +365*sigma/2 +365*(i -1),1)
            total_high_period = append(total_high_period,high_period)
            
        }
        
        y = sin(2*pi*cycle/(365*sigma))
        y[-total_high_period] = 0
    }
    
    return(y)
}


seasonal_sirs = function(t,y,params){
    
    S = y[1]
    I = y[2]
    R = y[3]
    
    beta = params["beta"]
    a = params["a"]
    sigma = params["sigma"]
    gamma = params["gamma"]
    vu = params["vu"]
    
    dS = -beta*(1 + a*sin_plus(t,sigma,0))*S*I + gamma*R
    dI = beta*(1 + a*sin_plus(t,sigma,0))*S*I -vu*I
    dR = vu*I - gamma*R
    
    result = c(dS,dI,dR)
    return(list(result))
}


neverfail_ode = function(out = T,params = params){
    
    times = seq(1,365*10,1)
    start = c(S = 0.999,I = 0.001,R = 0)
    
    while (out) {
        
        out = tryCatch(expr = ode(y = start, times = times, func = seasonal_sirs,parms = params),
                       warning = function(x) T)
        
        if (inherits(out, 'matrix')) {
            break
        } 
    }
    
    return(out)
}

neverfail_ode_wrapper = function(placeholder,out = T,params = params){
    out = neverfail_ode(out = T,params = params)
    out = data.frame(out)
    
    return(out$I)
}

report_sampling_wrapper = function(placeholder,I,n = 1,rr = rr){
    
    reported = sapply(I,rbinom,n = 1,prob = rr)
    
    return(reported)
}

params = c(beta = 2/5, a = 0.01,sigma = 1,vu = 1/5,gamma = 1/365)
out = neverfail_ode(out = T,params = params)
out = neverfail_ode_wrapper(1,out = T,params = params)

estimate = function(params,data = ili_plus_num_smoothed$ili_plus){
    
    a = params["a"]
    sigma = params["sigma"]
    rr = params["rr"]
    sd = params["sd"]
    
    times = seq(1,365*10,1)
    params = c(beta = 2/5, a,sigma,vu = 1/5,gamma = 1/365)
    
    start = c(S = 0.999,I = 0.001,R = 0)
    
    #### Integrate ODE
    # out = pbmclapply(1:100,neverfail_ode_wrapper,out = T,params = params,mc.cores = 8)
    out = neverfail_ode(out = T,params = params)
    
    out = data.frame(out)
    
    I = out$I[1000:length(out$I)]
    # m_I = rowMeans(out$I)
    
    # I = m_I[1000:length(m_I)]
    
    ###### Convert I from fraction to number 
    
    ### Choose the digits manually
    I = ceiling(I*100000)
    
    reported = pbmclapply(1:100,report_sampling_wrapper,I,n = 1,rr = rr,mc.cores = 8)
    reported = data.frame(reported)
    
    m_reported = ceiling(rowMeans(reported))
    
    nll = c()
    
    for(i in 1:length(data)){
        
        nll[i] = dnorm(data[i],m_reported[i],sd = sd,log = T)
    }
    
    
    return(sum(-nll))
    
}

start = c(a = 0.01,sigma = 0.5,rr = 0.1,sd = 1)
estimate(start,ili_plus_num_smoothed$ili_plus)

Sys.time()
fit = optim(par = start,fn = estimate,data = ili_plus_num_smoothed$ili_plus)
Sys.time()
