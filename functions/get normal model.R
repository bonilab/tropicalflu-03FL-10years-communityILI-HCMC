#################### lognormal regression pipeline #############################
# 0. log transfer the response 
# 1. fit the exhaustive model including all the predictors 
# 2. AIC-based model selection
# 3. Add non-linear terms based on AIC
# 4. AIC-based model selection

source('functions/predictor transformation function.R')
library(stringr)

get_normal_model = function(response_name, data) {
    
    # remove any NA/Inf #
    
    if(length(which(is.na(data[,response_name])))>0) {
        data = data[-which(is.na(data[,response_name])),]
    }
    if(length(which(is.infinite(data[,response_name])))>0) {
        data = data[-which(is.infinite(data[,response_name]))]
    }
    
    
    ### null model ###
    null_model = do.call(lm,list(formula = paste(response_name,'~1'),
                                 data = data))
    ## Step 1: fit the exhaustive model including all the predictors ##
    fit = do.call(lm,list(formula = paste(response_name,'~.'),
                          data = data))
    
    ## Step 2: AIC-based model selection ##
    fit2 = step(null_model, scope = list(lower = null_model, upper = fit))
    
    ## Step 3: Add non-linear terms ##
    transforms = c("square root","log","inverse","cube root",
                   "fourth root","fifth root")
    
    ### Only transform predictors, remove intercept ###
    predictors = names(fit2$coefficients)
    removed_predictors = c("(Intercept)")
    
    if(length(which(predictors %in% removed_predictors))>0){
        predictors = predictors[-which(predictors %in% removed_predictors)]
    }
    
    
    better_fits = list()
    better_fits_df = data.frame()
    for(transform in transforms){
        
        i = which(transforms == transform)
        better_fits[[i]] = transform_predictors(predictors = predictors,
                                                response = response_name,
                                                family = 'gaussian',
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
            
            fit3 = do.call(lm,list(formula = paste0(response_name,'~',selected_fits),
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
                                                    family = 'gaussian',
                                                    data = data)
            
            other_predictors = all_predictors_list[duplicated(all_predictors_list)]
            
            fit3_predictor_formula = paste(c(all_non_linear_terms,other_predictors),collapse = "+")
            
            fit3 = do.call(lm,list(formula = paste0(response_name,'~',fit3_predictor_formula),
                                   data = data))
        }
        
    } else {
        
        fit3 = fit2
    }
    
    
    ### Step 4: AIC-based model selection ###
    fit4 = step(null_model, scope = list(lower = null_model, upper = fit3))
    
    return(fit4)
}


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