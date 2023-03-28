diff_rss = function(predictor_to_assess,model,data) {
    
    rss_model = sum(resid(model)^2)
    
    model_adding_pred_call = paste(response_name,'~',paste0(names(coef(model))[-1],collapse = '+'),
                                   '+',predictor_to_assess)
    
    model_adding_pred = lm(model_adding_pred_call,data = data)
    
    rss_model_adding_pred = sum(resid(model_adding_pred)^2)
    
    rss_model - rss_model_adding_pred
}

diff_aic_forward = function(predictor_to_assess,response_name,model,data) {
    
    aic_model = AIC(model)
    
    model_adding_pred_call = paste(response_name,'~',paste0(names(coef(model))[-1],collapse = '+'),
                                    '+',predictor_to_assess)
    
    model_adding_pred = lm(model_adding_pred_call,data = data)
    
    aic_model_adding_pred = AIC(model_adding_pred)
    
    aic_model - aic_model_adding_pred
}

diff_aic_backward = function(predictor_to_assess,response_name,full_model,data) {
    
    aic_full_model = AIC(full_model)
    
    preds_full_model = names(coef(full_model))[-1]
    
    if(predictor_to_assess %in% preds_full_model) {
        
        preds_current_model = preds_full_model[-which(preds_full_model == predictor_to_assess)]
    } else {
        
        stop(paste('predictor to assess is not in the model'))
    }
    
    model_subtracting_pred_call = paste(response_name,'~',paste0(preds_current_model,collapse = '+'))
    
    model_subtracting_pred = lm(model_subtracting_pred_call,data = data)
    
    aic_model_subtracting_pred = AIC(model_subtracting_pred)
    
    aic_full_model - aic_model_subtracting_pred
}

diff_aic_gamma_hurdle_forward = function(predictor_to_assess,gamma_hurdle_model,data) {
    
    aic_model = AIC(gamma_hurdle_model)
    
    logi_call_adding_predictor = paste(as.character(gamma_hurdle_model$call)[5],'+',predictor_to_assess)
    
    gamma_call_adding_predictor = paste(as.character(gamma_hurdle_model$call)[2],'+',predictor_to_assess)
    
    gamma_hurdle_adding_predictor = glmmTMB(formula = formula(gamma_call_adding_predictor),
                                            ziformula = formula(logi_call_adding_predictor),
                                            family = ziGamma('log'),
                                            data = data,
                                            control = glmmTMBControl(optCtrl = list(iter.max = 1e4,eval.max = 1e4)))
    
    aic_gamma_hurdle_adding_predictor = AIC(gamma_hurdle_adding_predictor)
    
    aic_model - aic_gamma_hurdle_adding_predictor
}

diff_aic_backward_glm = function(predictor_to_assess,response_name,full_model,family, data) {
    
    aic_full_model = AIC(full_model)
    
    preds_full_model = names(coef(full_model))[-1]
    
    if(predictor_to_assess %in% preds_full_model) {
        
        preds_current_model = preds_full_model[-which(preds_full_model == predictor_to_assess)]
        
    } else {
        
        stop(paste('predictor to assess is not in the model'))
    }
    
    model_subtracting_pred_call = paste(response_name,'~',paste0(preds_current_model,collapse = '+'))
    
    model_subtracting_pred = glm(model_subtracting_pred_call,family = family, data = data)
    
    aic_model_subtracting_pred = AIC(model_subtracting_pred)
    
    aic_full_model - aic_model_subtracting_pred
}
