# load('../Rdata/Data and Fits for Regression/Data of overall iliplus.UV_absent_weather.AR_21.human for gamma hurdle model.Rdata')
# source('functions/plotting.R')
# 
# data_logistic = data_for_fit
# data_logistic$no_flu = ifelse(data_logistic$smth_7d_iliplus == 0,1,0)
# data_logistic$smth_7d_iliplus = NULL
# 
# non_zero_data = data_for_fit[data_for_fit$smth_7d_iliplus != 0,]

source('functions/predictor transformation function.R')
library(stringr)

check_predictors = function(all_non_linear_terms,selected_fits,response_name,family,data) {
    
    original_predictors = gsub('.*\\((.*)\\).*','\\1',all_non_linear_terms)
    
    ### remove comma if there is any ###
    if(length( grep('.*,.*',original_predictors)) >0) {
        
        original_predictors = gsub('(.*),.*','\\1',original_predictors)
    }
    
    if(any(duplicated(original_predictors))) {
        
        unique_original_predictors = unique(original_predictors)
        all_selected_transforms = c()
        for(i in 1:length(unique_original_predictors)) {
            pattern = paste0('.*',unique_original_predictors[i],'.*')
            this_non_linear_term = grep(pattern,all_non_linear_terms,value = T)
            
            aic_each_transform = data.frame()
            for(j in 1:length(this_non_linear_term)) {
                
                non_linear_pattern = paste0('.*',this_non_linear_term[j],'.*')
                non_linear_pattern = str_replace(non_linear_pattern,'\\(','\\\\(')
                non_linear_pattern = str_replace(non_linear_pattern,'\\)','\\\\)')
                
                this_formula = grep(non_linear_pattern,selected_fits,value = T)
                this_fit = glm(paste0(response_name,'~',this_formula),
                               family = family,
                               data = data)
                
                aic_each_transform = rbind(aic_each_transform,
                                           data.frame(term = this_non_linear_term[j],
                                                      aic = AIC(this_fit)))
                
                
            }
            selected_transform = aic_each_transform$term[which.min(aic_each_transform$aic)]
            all_selected_transforms = c(all_selected_transforms,selected_transform)
            
        }
        
        return(all_selected_transforms)
        
    } else {
        
        return(all_non_linear_terms)
    }
    
}

##### logistic model on the probability of absence of flu ###
logistic_model_selection = function(data,response_name) {
    
    ### null model ###
    null_model = do.call(glm,list(formula = paste(response_name,'~1'),
                                  family = 'binomial',
                                  data = data))
    
    fit = do.call(glm,list(formula = paste(response_name,'~.'),
                           family = 'binomial',
                           data = data))
    
    ### AIC-based forward model selection ###
    fit2 = step(null_model,scope = list(lower = null_model,upper = fit),
                trace = 0)
    
    ### Add non-linear transformation if AIC decreased ####
    transforms = c("square root","log","inverse","cube root",
                   "fourth root","fifth root")
    
    ### Only transform predictors, remove intercept ###
    predictors = names(fit2$coefficients)
    removed_predictors = c("(Intercept)",'iliplus_lag_21d','sim_ts_330','sim_ts_385','school')
    
    if(length(which(predictors %in% removed_predictors))>0){
        predictors = predictors[-which(predictors %in% removed_predictors)]
    }
    
    
    better_fits = list()
    better_fits_df = data.frame()
    for(transform in transforms){
        
        i = which(transforms == transform)
        better_fits[[i]] = transform_predictors(predictors = predictors,
                                                response_name = response_name,
                                                family = 'binomial',
                                                transform = transform,
                                                data = data,
                                                fit_aic = AIC(fit2))
        if(length(better_fits[[i]])>1) {
            better_fits_df = rbind(better_fits_df,
                                   data.frame(formula = better_fits[[i]]$formula,
                                              aic = better_fits[[i]]$aic))
        }
    }
    
    
    #### Only the non-linear terms that decreases AIC more than 9 points will be chosen ###
    if(any(better_fits_df$aic < AIC(fit2)-9)) {
        
        selected_fits = better_fits_df$formula[which(better_fits_df$aic < AIC(fit2) - 9)]
        
        if(length(selected_fits) == 1) {
            
            fit3 = do.call(glm,list(formula = paste0(response_name,'~',selected_fits),
                                    family = 'binomial',
                                    data = data))
        } else {
            
            ### Extract the predictors from formula ### 
            all_non_linear_terms = c()
            all_predictors_list = c()
            
            for(j in 1:length(selected_fits)) {
                
                predictors_list = strsplit(selected_fits[j],split = '+',fixed = T)[[1]]
                
                ## non-linear terms ##
                non_linear_term = grep('.*\\(.*\\).*',predictors_list,value = T)
                all_non_linear_terms = c(all_non_linear_terms,non_linear_term)
                all_predictors_list = c(all_predictors_list,predictors_list)
            }
            
            #### if there are multiple transformation on one predictor, choose the 
            #### the transformation that has the lowest AIC ###
            all_non_linear_terms = check_predictors(all_non_linear_terms = all_non_linear_terms,
                                                    selected_fits = selected_fits,
                                                    response_name = response_name,
                                                    family = 'binomial',
                                                    data = data)
            
            other_predictors = all_predictors_list[duplicated(all_predictors_list)]
            
            fit3_predictor_formula = paste(c(all_non_linear_terms,other_predictors),collapse = "+")
            
            fit3 = do.call(glm,list(formula = paste0(response_name,'~',fit3_predictor_formula),
                                    family = 'binomial',
                                   data = data))
        }
        
    } else {
        
        fit3 = fit2
    }
    
    
    ### Step 4: AIC-based model selection ###
    fit4 = step(null_model, scope = list(lower = null_model, upper = fit3))
    
    return(fit4)
}


# logistic_model = logistic_model_selection(data = data_logistic[,-1],
#                                           response = 'no_flu')

gamma_model_selection = function(data,response_name) {
    
    ## null model ##
    null_model = do.call(glm,list(formula = paste(response_name,'~1'),
                                  family = Gamma(link = 'log'),
                                  data = data))
    
    ## exhaustive model ##
    fit = do.call(glm,list(formula = paste(response_name,'~.'),
                           family = Gamma(link = 'log'),
                           data = data))
    
    ## AIC-based forward model selection ##
    fit2 = step(null_model,scope = list(lower = null_model,upper = fit),
                trace = 0)
    
    ### Add non-linear transformation if AIC decreased ####
    transforms = c("square root","log","inverse","cube root",
                   "fourth root","fifth root")
    
    ### Only transform predictors, remove intercept ###
    predictors = names(fit2$coefficients)
    removed_predictors = c("(Intercept)",'iliplus_lag_21d','sim_ts_330','sim_ts_385','school')
    
    if(length(which(predictors %in% removed_predictors))>0){
        predictors = predictors[-which(predictors %in% removed_predictors)]
    }
    
    
    better_fits = list()
    better_fits_df = data.frame()
    for(transform in transforms){
        
        i = which(transforms == transform)
        better_fits[[i]] = transform_predictors(predictors = predictors,
                                                response_name = response_name,
                                                family = Gamma(link = 'log'),
                                                transform = transform,
                                                data = data,
                                                fit_aic = AIC(fit2))
        if(length(better_fits[[i]])>1) {
            better_fits_df = rbind(better_fits_df,
                                   data.frame(formula = better_fits[[i]]$formula,
                                              aic = better_fits[[i]]$aic))
        }
    }
    
    
    #### Only the non-linear terms that decreases AIC more than 9 points will be chosen ###
    if(any(better_fits_df$aic < AIC(fit2)-9)) {
        
        selected_fits = better_fits_df$formula[which(better_fits_df$aic < AIC(fit2) - 9)]
        
        if(length(selected_fits) == 1) {
            
            fit3 = do.call(glm,list(formula = paste0(response_name,'~',selected_fits),
                                    family = Gamma(link = 'log'),
                                    data = data))
        } else {
            
            ### Extract the predictors from formula ### 
            all_non_linear_terms = c()
            all_predictors_list = c()
            
            for(j in 1:length(selected_fits)) {
                
                predictors_list = strsplit(selected_fits[j],split = '+',fixed = T)[[1]]
                
                ## non-linear terms ##
                non_linear_term = grep('.*\\(.*\\).*',predictors_list,value = T)
                all_non_linear_terms = c(all_non_linear_terms,non_linear_term)
                all_predictors_list = c(all_predictors_list,predictors_list)
            }
            
            #### if there are multiple transformation on one predictor, choose the 
            #### the transformation that has the lowest AIC ###
            all_non_linear_terms = check_predictors(all_non_linear_terms = all_non_linear_terms,
                                                    selected_fits = selected_fits,
                                                    response_name = response_name,
                                                    family = Gamma(link = 'log'),
                                                    data = data)
            
            other_predictors = all_predictors_list[duplicated(all_predictors_list)]
            
            fit3_predictor_formula = paste(c(all_non_linear_terms,other_predictors),collapse = "+")
            
            fit3 = do.call(glm,list(formula = paste0(response_name,'~',fit3_predictor_formula),
                                    family = Gamma(link = 'log'),
                                   data = data))
        }
        
    } else {
        
        fit3 = fit2
    }
    
    
    ### Step 4: AIC-based model selection ###
    fit4 = step(null_model, scope = list(lower = null_model, upper = fit3))
    
    ## Add interaction terms ##
    fit4_formula = as.list(fit4$call)$formula
    fit4_predictors = as.character(fit4_formula)[3]
    
    all_interactions_formula = paste0('(',fit4_predictors,')^2')
    full_model = do.call(glm,list(formula = paste(response_name,'~',all_interactions_formula),
                                  family = Gamma(link = 'log'),
                                  data = data))
    
    fit5 = step(null_model,scope = list(lower = null_model, upper = full_model),
                trace = 0)
    
    return(fit5)
    
}

gamma_model_selection_no_interaction = function(data,response_name) {
    
    ## null model ##
    null_model = do.call(glm,list(formula = paste(response_name,'~1'),
                                  family = Gamma(link = 'log'),
                                  data = data))
    
    ## exhaustive model ##
    fit = do.call(glm,list(formula = paste(response_name,'~.'),
                           family = Gamma(link = 'log'),
                           data = data))
    
    ## AIC-based forward model selection ##
    fit2 = step(null_model,scope = list(lower = null_model,upper = fit),
                trace = 0)
    
    
    ### Add non-linear transformation if AIC decreased ####
    transforms = c("square root","log","inverse","cube root",
                   "fourth root","fifth root")
    
    ### Only transform predictors, remove intercept ###
    predictors = names(fit2$coefficients)
    removed_predictors = c("(Intercept)",'iliplus_lag_21d','sim_ts_330','sim_ts_385','school')
    
    if(length(which(predictors %in% removed_predictors))>0){
        predictors = predictors[-which(predictors %in% removed_predictors)]
    }
    
    
    better_fits = list()
    better_fits_df = data.frame()
    for(transform in transforms){
        
        i = which(transforms == transform)
        better_fits[[i]] = transform_predictors(predictors = predictors,
                                                response_name = response_name,
                                                family = Gamma(link = 'log'),
                                                transform = transform,
                                                data = data,
                                                fit_aic = AIC(fit2))
        if(length(better_fits[[i]])>1) {
            better_fits_df = rbind(better_fits_df,
                                   data.frame(formula = better_fits[[i]]$formula,
                                              aic = better_fits[[i]]$aic))
        }
    }
    
    
    #### Only the non-linear terms that decreases AIC more than 9 points will be chosen ###
    if(any(better_fits_df$aic < AIC(fit2)-9)) {
        
        selected_fits = better_fits_df$formula[which(better_fits_df$aic < AIC(fit2) - 9)]
        
        if(length(selected_fits) == 1) {
            
            fit3 = do.call(glm,list(formula = paste0(response_name,'~',selected_fits),
                                    family = Gamma(link = 'log'),
                                    data = data))
        } else {
            
            ### Extract the predictors from formula ### 
            all_non_linear_terms = c()
            all_predictors_list = c()
            
            for(j in 1:length(selected_fits)) {
                
                predictors_list = strsplit(selected_fits[j],split = '+',fixed = T)[[1]]
                
                ## non-linear terms ##
                non_linear_term = grep('.*\\(.*\\).*',predictors_list,value = T)
                all_non_linear_terms = c(all_non_linear_terms,non_linear_term)
                all_predictors_list = c(all_predictors_list,predictors_list)
            }
            
            #### if there are multiple transformation on one predictor, choose the 
            #### the transformation that has the lowest AIC ###
            all_non_linear_terms = check_predictors(all_non_linear_terms = all_non_linear_terms,
                                                    selected_fits = selected_fits,
                                                    response_name = response_name,
                                                    family = Gamma(link = 'log'),
                                                    data = data)
            
            other_predictors = all_predictors_list[duplicated(all_predictors_list)]
            
            fit3_predictor_formula = paste(c(all_non_linear_terms,other_predictors),collapse = "+")
            
            fit3 = do.call(glm,list(formula = paste0(response_name,'~',fit3_predictor_formula),
                                    family = Gamma(link = 'log'),
                                   data = data))
        }
        
    } else {
        
        fit3 = fit2
    }
    
    
    ### Step 4: AIC-based model selection ###
    fit4 = step(null_model, scope = list(lower = null_model, upper = fit3))
    
    
    return(fit4)
    
}

# gamma_model = gamma_model_selection(data = non_zero_data[,-1],
#                                     response = 'smth_7d_iliplus')

#### Combine logistic model and gamma model together ###
library(glmmTMB)
get_gamma_hurdle_model = function(logistic_model,gamma_model,data) {
    
    logistic_call = as.list(logistic_model$call)
    logistic_formula = as.character(logistic_call$formula)[3]
    
    gamma_call = as.list(gamma_model$call)
    gamma_formula = gamma_call$formula
    
    gamma_hurdle = glmmTMB(formula = gamma_formula,
                           ziformula = formula(paste('~',logistic_formula)),
                           family = ziGamma(link = 'log'),
                           data = data,
                           control = glmmTMBControl(optCtrl = list(iter.max = 1e4,eval.max = 1e4)))
    
    ### Check model convergence : TRUE indicates convergence ###
    if(gamma_hurdle$sdr$pdHess) {
        print("Model Converges Successfully!")
    } else {
        print('Model Fails to Converge!')
    }
    
    return(gamma_hurdle)
}

fit_gamma_hurdle_model = function(data,response_name,interaction = T) {
    
    ### Get data for logistic model ###
    data_logi = data
    data_logi$zero = ifelse(data[,response_name] == 0,1,0)
    data_logi[,response_name] = NULL
    
    ### Get data for gamma model ###
    non_zero_data = data[data[,response_name] != 0,]
    
    ### Logistic model ###
    m_logi = logistic_model_selection(data = data_logi,response_name = 'zero')
    
    ### gamma model ###
    if(interaction == T) {
        
        m_gamma = gamma_model_selection(data = non_zero_data,response_name = response_name)
        
    } else {
        
        m_gamma = gamma_model_selection_no_interaction(data = non_zero_data,
                                                       response_name = response_name)
    }
    
    ### gamma hurdle model ###
    m_gamma_hurdle = get_gamma_hurdle_model(logistic_model = m_logi,
                                            gamma_model = m_gamma,
                                            data = data)
    return(m_gamma_hurdle)
}

### Use pseudo-r2 to evaluate goodness of fit ###
get_pseudo_r2_ols = function(model,response,n_parms) {
    
    ### r2 = 1 - RSS/SST ###
    RSS = sum(resid(model)^2)
    SST = sum((response - mean(response))^2)
    
    r2 = 1 - RSS/SST
    ### Adj-r2 = 1 - [(1 - r2)*(n - 1)/(n - k - 1)] ###
    n = length(response)
    k = n_parms
    
    adj_r2 = 1 - ((1 - r2)*(n - 1)/(n - k - 1))
    
    return(round(adj_r2,2))    
    
}

### get logistic model given call ###
get_logi = function(data,response_name,logi_call) {
    
    data$present = ifelse(data[,response_name] !=0,1,0)
        
    data[,response_name] = NULL
    
    logi_model = glm(paste0('present',logi_call),
                          family = 'binomial',
                          data = data)
    
    return(logi_model)
    
}

### get gamma model given call ###
get_gamma = function(data, response_name, gamma_call) {
    
    data = data[-which(data[,response_name] == 0),]
    
    gamma_model = glm(gamma_call,
                      family = Gamma('log'),
                      data = data)
    
    return(gamma_model)
    
}



