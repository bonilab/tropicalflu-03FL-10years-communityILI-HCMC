
# ###### Use method lsode because default method lsoda cannot solve ODE when beta is large
# ### lsode is used to solve stiff ODE, lsoda switches between stiff and nonstiff
# ### More need to understand
# 
library(deSolve)
sirs_1R = function(t,y,params){

    S = y[1]
    I = y[2]
    R = y[3]
    J = y[4]

    beta = params["beta"]
    gamma = params["gamma"]
    vu = params["vu"]

    dS = -beta*S*I + vu*R
    dI = beta*S*I - gamma*I
    dR = gamma*I - vu*R
    dJ = beta*S*I

    result = c(dS,dI,dR,dJ)
    return(list(result))
}

sirs_1R_solve = function(beta_cand,burn_in,total){

    parms = c(beta = beta_cand/10e7, vu = 1/365, gamma = 1/5)
    start = c(S = 10e7,I = 3000,R = 0,J = 0)
    times = seq(1,total,1)
    out = ode(y = start, times = times, parms = parms,func = sirs_1R,method = "lsode")
    out = data.frame(out)

    #### Exclude burn-in period
    out = out[burn_in:total,]

    return(out)

}


out = sirs_1R_solve(beta_cand = 0.3,burn_in = 5000,total = 10000)
plot(out$I)

# 
sirs_2R = function(t,y,params){

    S = y[1]
    I = y[2]
    R1 = y[3]
    R2 = y[4]
    # N = 10e7 + 3000
    J = y[5]


    beta = params["beta"]
    gamma = params["gamma"]
    vu = params["vu"]
    mu = params["mu"]

    dS = -beta*S*I + vu*R2
    dI = beta*S*I - gamma*I
    dR1 = gamma*I - vu*R1
    dR2 = vu*R1 - vu*R2
    dJ = beta*S*I

    result = c(dS,dI,dR1,dR2,dJ)

    return(list(result))
}

sirs_2R_solve = function(beta_cand,burn_in,total){

    params = c(beta = beta_cand/10e7,vu = 2/365,gamma = 1/5)
    start = c(S = 10e7,I = 3000, R1 = 0,R2 = 0,J = 0)
    time = seq(1,total,1)

    out = ode(y = start, times = time, func = sirs_2R, parms = params,method = "lsode")
    out = data.frame(out)


    #### Exclude burn-in period
    out = out[burn_in:total,]

    return(out)

}

out_2R = sirs_2R_solve(beta_cand = 0.3,burn_in = 1000,total = 5000)
plot(out_2R$I)

# 
# 
sirs_3R = function(t,y,params){

    S = y[1]
    I = y[2]
    R1 = y[3]
    R2 = y[4]
    R3 = y[5]
    # N = 10e7 + 3000
    J = y[6]


    beta = params[["beta"]]
    gamma = params[["gamma"]]
    vu = params[["vu"]]
    # mu = params["mu"]

    dS = -beta*S*I + vu*R3
    dI = beta*S*I - gamma*I
    dR1 = gamma*I - vu*R1
    dR2 = vu*R1 - vu*R2
    dR3 = vu*R2 - vu*R3
    dJ = beta*S*I

    result = c(dS, dI,dR1,dR2,dR3,dJ)

    return(list(result))
}

sirs_3R_solve = function(beta_cand,burn_in,total){

    params = c(beta = beta_cand/10e7,vu = 3/365,gamma = 1/5)
    start = c(S = 10e7,I = 3000, R1 = 0,R2 = 0,R3 = 0,J = 0)
    start_dy = c(S = 10e7,I = 3000, J = 0,R1 = 0,R2 = 0,R3 = 0)
    time = seq(1,total,1)

    out1 = ode(y = start, times = time, func = sirs_3R, parms = params,method = "lsode")
    out2 = ode(y = start_dy, times = time, func = sirs_dynamic_R, parms = params,method = "lsode")
    out = data.frame(out)
    #### Exclude burn-in period
    out = out[burn_in:total,]

    return(out)
}


# sirs_4R = function(t,y,params){
#     
#     S = y[1]
#     I = y[2]
#     R1 = y[3]
#     R2 = y[4]
#     R3 = y[5]
#     R4 = y[6]
#     # N = 10e7 + 3000
#     J = y[7]
#    
#     
#     
#     beta = params["beta"]
#     gamma = params["gamma"]
#     vu = params["vu"]
#     # mu = params["mu"]
#     
#     dS = -beta*S*I + vu*R4 
#     dI = beta*S*I - gamma*I
#     dR1 = gamma*I - vu*R1
#     dR2 = vu*R1 - vu*R2
#     dR3 = vu*R2 - vu*R3
#     dR4 = vu*R3 - vu*R4
#     dJ = beta*S*I
#     
#     result = c(dS,dI,dR1,dR2,dR3,dR4,dJ)
#     
#     return(list(result))
# }
# 
# sirs_4R_solve = function(beta_cand,burn_in,total){
#     params = c(beta = beta_cand/10e7,vu = 4/365,gamma = 1/5)
#     start = c(S = 10e7,I = 3000, R1 = 0,R2 = 0,R3 = 0,R4 = 0,J = 0)
#     time = seq(1,total,1)
#     
#     out = ode(y = start, times = time, func = sirs_4R, parms = params,method = "lsode")
#     out = data.frame(out)
#     #### Exclude burn-in period
#     out = out[burn_in:total,]
#     
#     return(out)
# }
# 
# sirs_5R = function(t,y,params){
# 
#     S = y[1]
#     I = y[2]
#     R1 = y[3]
#     R2 = y[4]
#     R3 = y[5]
#     R4 = y[6]
#     R5 = y[7]
#     # N = 10e7 + 3000
#     J = y[8]
# 
# 
#     beta = params["beta"]
#     gamma = params["gamma"]
#     vu = params["vu"]
#     # mu = params["mu"]
# 
#     dS = -beta*S*I + vu*R5
#     dI = beta*S*I - gamma*I
#     dR1 = gamma*I - vu*R1
#     dR2 = vu*R1 - vu*R2
#     dR3 = vu*R2 - vu*R3
#     dR4 = vu*R3 - vu*R4
#     dR5 = vu*R4 - vu*R5
#     dJ = beta*S*I
# 
#     result = c(dS,dI,dR1,dR2,dR3,dR4,dR5,dJ)
# 
#     return(list(result))
# }
# 
# sirs_5R_solve = function(beta_cand,burn_in,total){
# 
#     params = c(beta = beta_cand/10e7,vu = 5/365,gamma = 1/5)
#     start = c(S = 10e7,I = 3000, R1 = 0,R2 = 0,R3 = 0,R4 = 0,R5 = 0,J = 0)
#     time = seq(1,total,1)
# 
#     out = ode(y = start, times = time, func = sirs_5R, parms = params,method = "lsode")
#     out = data.frame(out)
# 
#     #### Exclude burn-in period
#     out = out[burn_in:total,]
# 
#     return(out)
# 
# }
# 

sirs_multiR = function(t, y, params) {
    
    if(names(y)[1] != 'S' | names(y)[2] != 'I' | names(y)[3] != 'J') {
        
        stop(paste('start index does not match'))
    }
    
    S = y[1]
    I = y[2]
    J = y[3]
    
    # Assume that from the 4th position in y, they are treated as R
    y_len = length(y)
    num_Rs = y_len - 3
    
    beta = params["beta"]
    gamma = params["gamma"]
    vu = params["vu"]
    
    dS = -beta*S*I + vu*y[y_len]
    dI = beta*S*I - gamma*I
    dJ = beta*S*I
    
    # Calculate the following terms
    #
    # dR1 = gamma * I     - vu  * R1
    # dR2 = vu    * R1    - vu  * R2
    # dR3 = vu    * R2    - vu  * R3
    # ...
    # dR  = parms1   * cp1   - parms2 * cp2
    
    if(num_Rs == 1) {
        
        dR1 = gamma*I - vu*y[y_len]
        
        result = c(dS,dI,dJ,dR1)
    } else {
        parms1 = c(gamma, rep(vu, length=num_Rs-1))
        cp1 = c(I, y[4:(y_len-1)])
        parms2 = rep(vu, length=num_Rs)
        cp2 = y[4:y_len]
        
        dR = parms1 * cp1 - parms2 * cp2
        result = c(dS,dI,dJ,dR)
        
    }
    
    return(list(result))
}

library(deSolve)
sirs_multiR_solve = function(num_Rs,beta,burn_in,total) {
    
    R = rep(0,num_Rs)
    names(R) = paste0("R",1:num_Rs)
    params = c(beta = beta/10e7,vu = num_Rs/365,gamma = 1/5)
    start = c(S = 10e7,I = 3000, J = 0, R)
    time = seq(1,total,1)
    
    out = ode(y = start, times = time, func = sirs_multiR, parms = params,method = "lsode")
    out = data.frame(out)
    
    #### Exclude burn-in period
    out = out[burn_in:total,]
    
    return(out)
    
}
