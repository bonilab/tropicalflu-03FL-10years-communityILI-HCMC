#### parametric boostrap for 95% PI of gamma hurdle model ####
rm(list = ls())
library(pbmcapply)
setwd("~/Dropbox/Influenza and Respiratory Disease in Vietnam/code")
load('Rdata/Overall.ILIplus.21-day.aggregated.pcr.RData')
source('functions/get gamma hurdle model.R')

#### cyc330 has the lowest AIC, load N-stepfunc from cyc330 cyc385 ###
load('Rdata/ILI+_best_stepfunc_for_sensitivity_analysis_cyc150_450.Rdata')
load('Rdata/Data of overall iliplus_UV_WS_absent_scaled_weather.AR_21.human for gamma hurdle model.Rdata')
ili_plus_21d_df$sim_ts_330 = all_iliplus_sim_ts[which(rownames(all_iliplus_sim_ts) == '330'),]
ili_plus_21d_df$sim_ts_385 = all_iliplus_sim_ts[which(rownames(all_iliplus_sim_ts) == '385'),]

data_for_fit$sim_ts_330 = ili_plus_21d_df$sim_ts_330[which(ili_plus_21d_df$date %in% data_for_fit$date)]
data_for_fit$sim_ts_385 = ili_plus_21d_df$sim_ts_385[which(ili_plus_21d_df$date %in% data_for_fit$date)]

### Get the params MLE from real data ###
full_gamma_hurdle = fit_gamma_hurdle_model(data = data_for_fit[,-1],
                                           response_name = 'smth_7d_iliplus',
                                           interaction = F)

### need to rebuild logistic model based on full gamma hurdle model ###

## get logistic model from gamma hurdle model ##
logi_model = get_logi(data = data_for_fit[,-1],
                         response_name = 'smth_7d_iliplus',
                         logi_call = as.character(full_gamma_hurdle$call[5]))

logi_call = as.character(full_gamma_hurdle$call[5])

### get gamma model from gamma hurdle model ###
gamma_model = get_gamma(data = data_for_fit[,-1],
                           response_name = 'smth_7d_iliplus',
                           gamma_call = as.character(full_gamma_hurdle$call[2]))

gamma_call = as.character(full_gamma_hurdle$call[2])

gamma_disp = sd(data_for_fit$iliplus_lag_21d)

### Get M parametric boostrap samples of data ##
predictor_data = data_for_fit[,-which(colnames(data_for_fit) == 'smth_7d_iliplus')]

M = 200

all_info = pbmclapply(1:M,function(m) {
    
    ### simulate M time series ###
    sim_resp = c()
    
    for(i in 1:nrow(data_for_fit)) {
        
        
        ### binomial draw ###
        pred_prob = predict(logi_model, 
                            type = 'response',
                            newdata = data_for_fit[i,])
        
        z = rbinom(1,size = 1, prob = pred_prob)
        
        if(z == 1) {
            
            pred_value = predict(gamma_model,
                                 type = 'response',
                                 newdata = data_for_fit[i,])
            
            ## re-parameterize shape and rate, shape = mean^2/sd^2, rate = mean/sd^2
            sim_resp[i] = rgamma(1,shape = pred_value^2/gamma_disp^2,
                                 rate = pred_value/gamma_disp^2)
        } else {
            sim_resp[i] = 0
        }
        
    }
   
    ### Get param estimate from sim data ###
    sim_data = cbind(smth_7d_iliplus = sim_resp,
                     predictor_data)
    
    ### get logistic model ###
    logi_model = get_logi(data = sim_data,
                          response_name = 'smth_7d_iliplus',
                          logi_call = logi_call)
    
    ### get gamma model ###
    gamma_model = get_gamma(data = sim_data,
                            response_name = 'smth_7d_iliplus',
                            gamma_call = gamma_call)
    
    return(list(sim_resp = sim_resp,
                logi_mod = logi_model,
                gamma_mod = gamma_model))

},mc.cores = 5)

save(all_info,file = 'Rdata/parametric_boostrap_model_TimeSeries.Rdata')

### Get prediction interval for each point based on 500 samples of params MLE ##
load('Rdata/parametric_boostrap_model_TimeSeries.Rdata')

all_sim_mle_data = pbmclapply(1:M, function(m) {
   
    logi_mod = all_info[[m]]$logi_mod
    gamma_mod = all_info[[m]]$gamma_mod
    
    ### simulate M time series ###
    sim_resp_mle = c()
    gamma_disp_mle = sd(all_info[[m]]$sim_resp,na.rm = T)
    
    for(i in 1:nrow(data_for_fit)) {
        
        ### binomial draw ###
        pred_prob_mle = predict(logi_mod,
                                type = 'response',
                                newdata = data_for_fit[i,])
        
        z = rbinom(1,size = 1, prob = pred_prob_mle)
        
        if(z == 1) {
            
            pred_value_mle = predict(gamma_mod,
                                     type = 'response',
                                     newdata = data_for_fit[i,])
            
            ## re-parameterize shape and rate, shape = mean^2/sd^2, rate = mean/sd^2
            sim_resp_mle[i] = rgamma(1,shape = pred_value_mle^2/gamma_disp_mle^2,
                                     rate = pred_value_mle/gamma_disp_mle^2)
            
        } else {
            sim_resp_mle[i] = 0
        }
        
    }
   
   return(sim_resp_mle) 
},mc.cores = 5)

all_sim_mle_data_df = do.call(rbind,all_sim_mle_data)
save(all_sim_mle_data_df,file = 'Rdata/simulated_mle_data.Rdata')

### 95% PI ###

fitted_data = data.frame(date = data_for_fit$date,
                         fitted = fitted(full_gamma_hurdle),
                         lower = NA,
                         upper = NA)

for(k in 1:nrow(data_for_fit)) {
    
    fitted_data$lower[k] = quantile(all_sim_mle_data_df[,k],0.025)
    fitted_data$upper[k] = quantile(all_sim_mle_data_df[,k],0.975)
    
}

plot(fitted_data$fitted,type = 'l',ylim = c(0,max(fitted_data$upper)))
lines(fitted_data$lower,col = 'red')
lines(fitted_data$upper,col = 'blue')

save(fitted_data,file = 'Rdata/ILI+_param_boostrap.Rdata')




