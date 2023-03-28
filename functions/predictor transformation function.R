library(pracma)
transform_predictors = function(predictors,response_name,family,transform,data,fit_aic){
    
    fit_list = list()
    allfits_aic = c()
    
    for(i in 1:length(predictors)){
        
        #### warning still pops up for NA/Inf after transformation ###
        if(transform == "square root"){
            ## Remove the predictors that will produce NA or Inf
            if(any(is.infinite(sqrt(data[,predictors[i]]))) || any(is.na(sqrt(data[,predictors[i]])))){
                break()
            } else {
                transformed_pred = paste0("sqrt(",predictors[i],")")
            }
        }
        
        if(transform == "2nd-order"){transformed_pred = paste0("I(",predictors[i],")^2")}
        
        if(transform == "log"){
            
            ## Remove the predictors that will produce NA or Inf
            if(any(is.infinite(log(data[,predictors[i]]))) || any(is.na(log(data[,predictors[i]])))){
                break()
            } else {
                transformed_pred = paste0("log(",predictors[i],")")
            }
        }
        
        if(transform == "inverse"){
            
            ## Remove the predictors that will produce NA or Inf
            if(any(is.infinite(1/(data[,predictors[i]]))) || any(is.na(1/(data[,predictors[i]])))){
                break()
            } else {
                transformed_pred = paste0("1/(",predictors[i],")")
            }
        }
        
        if(transform == "cube root"){transformed_pred = paste0("nthroot(",predictors[i],",3)")}
        
        if(transform == "fourth root") {
            ## Remove the predictors that will produce NA or Inf
            if(any(is.infinite(sqrt(data[,predictors[i]]))) || any(is.na(sqrt(data[,predictors[i]])))){
                break()
            } else {
                transformed_pred = paste0("nthroot(",predictors[i],",4)")
            }
            
        }
        
        if(transform == "fifth root"){transformed_pred = paste0("nthroot(",predictors[i],",5)")}
        
        form = ""
        
        for (j in 1:length(predictors)){
            
            if(j == i){
                form = paste(form,"+",transformed_pred)
            } else {
                
                predictor = paste(predictors[j])
                form = paste(form,"+",predictor)
            }
        }
        
        formula = gsub("^ \\+ (.*)",paste0(response_name,'~\\1'),form)
        fit_list[[i]] = glm(formula = formula,family = family,data)
        # fit_list[[i]] = do.call(glm,list(formula = formula,family = family, data = data))
        allfits_aic[i] = AIC(fit_list[[i]])
    }
    ind = which(allfits_aic< fit_aic)
    
    if(length(ind) >0){
        selected_fit = lapply(ind,function(x){fit_list[[x]]})
        
        return(list(formula = sapply(1:length(selected_fit),function(x){
            
            paste(names(selected_fit[[x]]$coefficients)[-1],collapse = '+')
        }),
        aic = sapply(1:length(selected_fit),function(x){AIC(selected_fit[[x]])})))
    } else {
        print(paste(transform,"does not give a lower AIC"))
    }
}


# transform_predictors_gamma_hurdle = function(ziformula,predictors,response_name,transform,data,fit_aic){
#     
#     fit_list = list()
#     allfits_aic = c()
#     
#     for(i in 1:length(predictors)){
#         
#         if(transform == "square root"){
#             ## Remove the predictors that will produce NA or Inf
#             if(any(is.infinite(sqrt(data[,predictors[i]]))) || any(is.na(sqrt(data[,predictors[i]])))){
#                 break()
#             } else {
#                 transformed_pred = paste0("sqrt(",predictors[i],")")
#             }
#         }
#         
#         if(transform == "2nd-order"){transformed_pred = paste0("I(",predictors[i],")^2")}
#         
#         if(transform == "log"){
#             
#             ## Remove the predictors that will produce NA or Inf
#             if(any(is.infinite(log(data[,predictors[i]]))) || any(is.na(log(data[,predictors[i]])))){
#                 break()
#             } else {
#                 transformed_pred = paste0("log(",predictors[i],")")
#             }
#         }
#         
#         if(transform == "inverse"){
#             
#             ## Remove the predictors that will produce NA or Inf
#             if(any(is.infinite(1/(data[,predictors[i]]))) || any(is.na(1/(data[,predictors[i]])))){
#                 break()
#             } else {
#                 transformed_pred = paste0("1/(",predictors[i],")")
#             }
#         }
#         
#         if(transform == "cube root"){transformed_pred = paste0("nthroot(",predictors[i],",3)")}
#         
#         if(transform == "fourth root") {
#             ## Remove the predictors that will produce NA or Inf
#             if(any(is.infinite(sqrt(data[,predictors[i]]))) || any(is.na(sqrt(data[,predictors[i]])))){
#                 break()
#             } else {
#                 transformed_pred = paste0("nthroot(",predictors[i],",4)")
#             }
#             
#         }
#         
#         if(transform == "fifth root"){transformed_pred = paste0("nthroot(",predictors[i],",5)")}
#         
#         form = ""
#         
#         for (j in 1:length(predictors)){
#             
#             if(j == i){
#                 form = paste(form,"+",transformed_pred)
#             } else {
#                 
#                 predictor = paste(predictors[j])
#                 form = paste(form,"+",predictor)
#             }
#         }
#         
#         formula = gsub("^ \\+ (.*)",paste0(response_name,'~\\1'),form)
#         formula = formula(formula)
#         fit_list[[i]] = glmmTMB(formula = formula, 
#                                 family = ziGamma(link = 'log'),
#                                 ziformula = ziformula,
#                                 data)
#         
#         allfits_aic[i] = AIC(fit_list[[i]])
#     }
#     ind = which(allfits_aic< fit_aic)
#     
#     if(length(ind) >0){
#         selected_fit = lapply(ind,function(x){fit_list[[x]]})
#         
#         return(list(formula = sapply(1:length(selected_fit),function(x){names(selected_fit[[x]]$coefficients)[-1]}),
#                     aic = sapply(1:length(selected_fit),function(x){AIC(selected_fit[[x]])})))
#     } else {
#         print(paste(transform,"does not give a lower AIC"))
#     }
# }


# transform_predictors_vglm = function(predictors,transform,data,fit_aic){
#     
#     fit_list = list()
#     allfits_aic = c()
#     
#     for(i in 1:length(predictors)){
#         
#         if(transform == "square root"){
#             ## Remove the predictors that will produce NA or Inf
#             if(any(is.infinite(sqrt(data[,predictors[i]]))) || any(is.na(sqrt(data[,predictors[i]])))){
#                 break()
#             } else {
#                 transformed_pred = paste0("sqrt(",predictors[i],")")
#             }
#         }
#         
#         if(transform == "2nd-order"){transformed_pred = paste0("I(",predictors[i],")^2")}
#         
#         if(transform == "log"){
#             
#             ## Remove the predictors that will produce NA or Inf
#             if(any(is.infinite(log(data[,predictors[i]]))) || any(is.na(log(data[,predictors[i]])))){
#                 break()
#             } else {
#                 transformed_pred = paste0("log(",predictors[i],")")
#             }
#         }
#         
#         if(transform == "inverse"){
#             
#             ## Remove the predictors that will produce NA or Inf
#             if(any(is.infinite(1/(data[,predictors[i]]))) || any(is.na(1/(data[,predictors[i]])))){
#                 break()
#             } else {
#                 transformed_pred = paste0("1/(",predictors[i],")")
#             }
#         }
#         
#         if(transform == "cube root"){transformed_pred = paste0("nthroot(",predictors[i],",3)")}
#         
#         if(transform == "fourth root") {
#             ## Remove the predictors that will produce NA or Inf
#             if(any(is.infinite(sqrt(data[,predictors[i]]))) || any(is.na(sqrt(data[,predictors[i]])))){
#                 break()
#             } else {
#                 transformed_pred = paste0("nthroot(",predictors[i],",4)")
#             }
#             
#         }
#         
#         if(transform == "fifth root"){transformed_pred = paste0("nthroot(",predictors[i],",5)")}
#         
#         form = ""
#         
#         for (j in 1:length(predictors)){
#             
#             if(j == i){
#                 form = paste(form,"+",transformed_pred)
#             } else {
#                 
#                 predictor = paste(predictors[j])
#                 form = paste(form,"+",predictor)
#             }
#         }
#         
#         formula = gsub("^ \\+ (.*)","iliplus~\\1",form)
#         fit_list[[i]] = vglm(formula = formula,
#                              family = tobit(Lower = 0,Upper = Inf, lmu = 'identitylink',type.fitted = 'censored'),
#                              data)
#         allfits_aic[i] = AIC(fit_list[[i]])
#     }
#     ind = which(allfits_aic< fit_aic)
#     
#     if(length(ind) >0){
#         selected_fit = lapply(ind,function(x){fit_list[[x]]})
#         
#         return(list(formula = sapply(1:length(selected_fit),function(x){names(coef(selected_fit[[x]]))[-c(1,2)]}),
#                     aic = sapply(1:length(selected_fit),function(x){AIC(selected_fit[[x]])})))
#     } else {
#         print(paste(transform,"does not give a lower AIC"))
#     }
# }
