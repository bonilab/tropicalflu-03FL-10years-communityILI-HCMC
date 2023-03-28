###### predictor importance on temperate ILI zeta-score pipeline ##############
library(ggplot2)
source('functions/Multiple_stepfunc_weekly_data.R')
source('functions/get normal model.R')
### 1. Get weather data ####

get_weather_data = function(weather_csv_location,
                            rows_to_skip,
                            ncol, 
                            plot_raw = T,
                            plot_weekly = T,
                            save = T,
                            filename) {
    
    weather = read.csv(weather_csv_location,skip = rows_to_skip)
    
    if(ncol(weather) != ncol) {
        
        stop('wrong csv reading!')
    }
    
    weather$date = as.Date(paste(weather$YEAR,weather$DOY,sep = "-"),format = "%Y-%j")
    
    ## check the weather data before running functions: plotting, any NA, any left out dates
    if(plot_raw) {
        
        par(mfrow = c(3,1),mar = c(2,2,2,2))
        plot(weather$date,weather$T2M,pch = 20, cex = 0.5,main = 'Temp')
        plot(weather$date,weather$RH2M,pch = 20, cex = 0.5,main = 'RH')
        plot(weather$date,weather$PRECTOTCORR,pch = 20, cex = 0.5,main = "RF")
        
    }
    
    complete_ts = seq(min(weather$date),max(weather$date),"1 day")
    
    if(nrow(weather) != length(complete_ts)) {
        
        diff_dates = complete_ts[which(!complete_ts %in% weather$date)]
        warning(paste('Missing dates in weather data:',diff_dates))
    }
    
    weather_temp = cbind(weather$T2M,weather$RH2M,weather$PRECTOTCORR)
    
    if(any(is.na(weather_temp)) | any(weather_temp == -99)){
        
        stop(paste('Missing data in the weather'))
    }
    
    ## 1a. Get daily AH from daily Temp and RH
    weather$AH = 0.611*exp(17.502*weather$T2M/(240.97+weather$T2M))*18.02*weather$RH2M/(8.31*(273.15+weather$T2M))
    
    ## 1b. Get weekly average weather data, save weekly weather data 
    weekly_dates = seq(min(weather$date),max(weather$date),'7 days')
    weather$ind = findInterval(weather$date,weekly_dates)
    index = unique(weather$ind)
    
    weekly_weather = data.frame()
    for(i in index) {
        
        weather_one_week = weather[weather$ind == i,]   
        
        #### take arithemic mean of the weather in one week ###
        weekly_weather = rbind(weekly_weather,
                               data.frame(date = weather_one_week$date[1],
                                          temp = mean(weather_one_week$T2M),
                                          AH = mean(weather_one_week$AH),
                                          RF = mean(weather_one_week$PRECTOTCORR)))
    }
    
    #### check the data by plotting ###
    if(plot_weekly) {
        
        par(mfrow = c(3,1),mar = c(2,2,2,2))
        plot(weekly_weather$date,weekly_weather$temp,pch = 20, cex = 0.5,main = 'weekly temp')
        plot(weekly_weather$date,weekly_weather$AH,pch = 20, cex = 0.5,main = 'weekly AH')
        plot(weekly_weather$date,weekly_weather$RF,pch = 20, cex = 0.5,main = 'weekly RF')
        
    }
    
    if(save) {
        
        save(weekly_weather,file = filename)
        print('Saving completed!')
    } else {
        
        return(weekly_weather)
    }
    
}

# get_weather_data(weather_csv_location = '../../_DATA/Weather Data NASA/POWER_Point_Daily_19971005_20191231_point_US_HHS_r1_boston.csv',
#                  save = T,
#                  filename = '../Rdata/Weather Data/weekly weather of Boston_19997_2019.Rdata')


### 2. Get step function ####
###### Read US HHS zeta step function given cycle 150 - 450 #################

get_stepfunc = function(results_location,
                        N_fits,
                        plot_stepfunc_grid = T,
                        data_to_fit,
                        save = T,
                        filename) {
    
    all_fits = list()
    
    for(i in 1:N_fits) {
        
        load(paste0(results_location,i,'.RData'))
        all_fits[[i]] = best_fit
        
    }
    
    all_N = sapply(1:length(all_fits),function(x) {all_fits[[x]]$N})
    all_cycle = sapply(1:length(all_fits),function(x) {all_fits[[x]]$cycle})
    all_aic = sapply(1:length(all_fits),function(x) {all_fits[[x]]$aic})
    
    all_results = data.frame(cycle = all_cycle,
                             N = all_N,
                             aic = all_aic)

    if(plot_stepfunc_grid) {
    
    selected_results = all_results[all_results$aic < 1e+32,]
    ggplot(selected_results) +
        geom_raster(aes(x = cycle, y = N, fill = aic)) +
        scale_fill_viridis_c(option = 'H')
    
    }
    
    #### save all the step functions ####
    cycle_aic = data.frame()
    all_sim_ts = c()
    
    for(cycle in unique(all_cycle)) {
        
        indices_aic = data.frame(index = which(all_results$cycle == cycle),
                                 aic = all_results$aic[all_results$cycle == cycle])
        index = indices_aic$index[which.min(indices_aic$aic)]
        
        cycle_aic = rbind(cycle_aic, data.frame(cycle = cycle, 
                                                aic = min(indices_aic$aic)))
        
        sim_ts = sapply(1:length(data_to_fit),multi_stepfunc,cycle = cycle,
                        all_fits[[index]]$all_results)
        
        all_sim_ts = rbind(all_sim_ts,sim_ts)
    }
    
    cycle_aic = cycle_aic[order(cycle_aic$aic),]
    rownames(all_sim_ts) = unique(all_cycle)
    
    stepfunc_info = list(cycle_aic = cycle_aic, all_sim_ts = all_sim_ts)
    
    if(save) {
        
        save(stepfunc_info, file = filename)
        
    } else {
        
        return(stepfunc_info)
    }
    
}
### 3. Generate regression dataset ####

# load('../Rdata/Weather Data/weekly weather of Boston_19997_2019.Rdata')
# load('../Rdata/Global ILI/US HHS Regions ILI percent and Zeta score.Rdata')
# load('../Rdata/sensitivity analysis_temperate regions/hhs_r1_stepfunc_info.Rdata')

# all_sim_ts = hhs_r1_info$all_sim_ts
# 
# annual_cycle = all_sim_ts[which(rownames(all_sim_ts) == '52'),]
# non_annual_cycle = all_sim_ts[which(rownames(all_sim_ts) == '27'),]
# 
# hhs1 = hhs_regions[[1]]
# 
# plot(hhs1$zeta_score,type = 'l')
# lines(annual_cycle,col = 'red')
# lines(non_annual_cycle,col = 'green')

get_regression_data = function(weather_data,ili_zeta_df,
                               annual_cycle,non_annual_cycle, 
                               time_name, response_name, 
                               save = T, filename) {
    
    ## 2a. Normalize the weather data
    scaled_weather =  data.frame(date = weather_data$date,
                                 apply(weather_data[,-1],2,scale,center = T,scale = T))
    
    ## 2b. create regression dataset, add zeta, weather, stepfunc
    ili_zeta_df$annual = annual_cycle
    ili_zeta_df$non_annual = non_annual_cycle
    
    #### SHOLD ADD NA FOR ANY MISSING DATA ####
    data_for_fit = data.frame(date = scaled_weather$date,
                              ili_zeta = NA,
                              annual = NA, 
                              non_annual = NA,
                              scaled_weather[,-1])
    
    for(i in 1:nrow(data_for_fit)) {
        
        ind = which(ili_zeta_df[,time_name] == data_for_fit$date[i])
        
        if(length(ind) == 0) {
            next
        } else {
            
            data_for_fit$ili_zeta[i] = ili_zeta_df[,response_name][ind]
            data_for_fit$annual[i] = ili_zeta_df$annual[ind]
            data_for_fit$non_annual[i] = ili_zeta_df$non_annual[ind]
        }
    }
    
    ## 2c. Add school term: May 1st - August 20th is holiday every year
    data_for_fit$school = 1
    data_for_fit$year = format(data_for_fit$date,"%Y")
    years = unique(data_for_fit$year)
    
    #### School is open from Aug 20 to April 30 every year
    for(year in years){
        
        holiday = seq(as.Date(paste0(year,"-05-01")),as.Date(paste0(year,"-08-19")),"1 day")
        data_for_fit$school[data_for_fit$date %in% holiday] = 0
        
    }
    
    data_for_fit$year = NULL
    
    ## 2d. Generate 1w-3w lagged weather data 
    max_lag = 3
    week_lags = seq(1,max_lag,1)
    
    weather_to_lag = data_for_fit[,which(colnames(data_for_fit) %in% c('temp','AH','RF'))]
    
    new_dataset = data_for_fit[(max_lag + 1):nrow(data_for_fit),]
    
    for(week_lag in week_lags){
        
        lagged_start = 1 + max_lag - week_lag
        lagged_end = lagged_start + nrow(new_dataset) -1
        
        lag_df = weather_to_lag[lagged_start:lagged_end,]
        
        colnames(lag_df) = paste0(colnames(weather_to_lag),'_lag_',week_lag,'wk')
        
        new_dataset = cbind(new_dataset,lag_df)
    }
    
    ## 2e. Add AR(1wk) of ILI zeta-score
    zeta_lagged_period = new_dataset$date - 7

    new_dataset$zeta_lag1w = data_for_fit$ili_zeta[which(data_for_fit$date %in% zeta_lagged_period)]
    
    ## 2f. log transform the response variable, save dataset ##
    #### NOT CONSIDERED BECAUSE OF ZERO IN THE DATA #######
    
    if(save) {
        
        save(new_dataset,file = filename)
        
    } else {
        
        return(new_dataset)
    }
}

# get_regression_data(weather_data = weekly_weather,
#                     ili_zeta_df = hhs1,
#                     annual_cycle = annual_cycle,
#                     non_annual_cycle = non_annual_cycle,
#                     filename = '../Rdata/sensitivity analysis_temperate regions/regression_dataset_hhs_r1.Rdata')


### 4. Multiple regression ####
## using get_normal_model, force non-annual step function to be in the model
## in case it gets removed by stepwise model selection



### No NA in the response and predictor ###
# new_dataset_complete = na.omit(new_dataset)
# full_model = get_normal_model(response_name = 'ili_zeta',
#                               data = new_dataset_complete[,-which(colnames(new_dataset) %in% c('date'))])
# 
# summary(full_model)
# 
# step(full_model)
# 
# fit = lm(ili_zeta~.,
#          data = new_dataset_complete[,-1])
# summary(fit)
# 
# null = lm(ili_zeta~1, data = new_dataset_complete[,-1])
# fit2 = step(null, scope = list(lower = null, upper = fit))
# 
# summary(fit2)
# step(fit2)
# 
# 
# plot(new_dataset_complete$ili_zeta,type = 'l')
# lines(new_dataset_complete$annual,col = 'red')
# 
# 
# 
# 
# 
# 
# ### 5. Predictor importance in the full model ####
# ## use functions in predictor importance.R ##
source('functions/predictor importance.R')

get_pred_impo_full_model = function(model,response_name,family,data,barplot = T) {
    
    predictors = names(coef(model)[-1])
    
    preds_aic = data.frame()
    for(predictor in predictors) {
        
        preds_aic = rbind(preds_aic,
                          data.frame(preds = predictor,
                                     diff_aic = diff_aic_backward_glm(predictor_to_assess = predictor,
                                                                      response_name = response_name,
                                                                      full_model = model,
                                                                      family = family,
                                                                      data = data)))
                                     
    }
    
    if(barplot) {
        barplot(-preds_aic$diff_aic,
                names.arg = preds_aic$preds)
        
    }
    
    return(preds_aic)
    
    
}
# 
# ##### 5. Relative importance ####
# library(relaimpo)
# 
# relaimpo_r1 = calc.relimp(full_model)
# 
# barplot(t@lmg)
# 
# 

