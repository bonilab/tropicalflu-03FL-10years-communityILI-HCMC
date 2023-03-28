rm(list = ls())
library(relaimpo)
source('functions/pred_impo_pipeline functions.R')
source('functions/predictor importance.R')
load('Rdata/Global ILI/US HHS Regions ILI percent and Zeta score.Rdata')
load('Rdata/Weekly zeta from averaging daily zeta in a week.Rdata')

names(hhs_regions) = paste0('HHS R',1:10)

regions_index = 1:10

#### Get the step function of HHS ###
hhs_regions_all_stepfunc = lapply(regions_index,function(x) {
    
    if(x == 1) {
        
        results_location = 'Rdata/Results from Roar/Roar_us_HHS_ILI_r1/best_fit_'
        
    } else {
        
        results_location = paste0('Rdata/Results from Roar/Roar_us_HHS_ILI_r',x,'/best_fit_r',x,'_')
    }
    
    hhs_stepfunc = get_stepfunc(results_location = results_location,
                                N_fits = 88, 
                                data_to_fit = hhs_regions[[x]]$zeta_score,
                                save = F)
    
    return(hhs_stepfunc)
})
names(hhs_regions_all_stepfunc) = paste0('r',regions_index)

#### check and choose the step function for each region ###
lapply(regions_index, function(x) {
    
    hhs_regions_all_stepfunc[[x]]$cycle_aic
})


cycle_indexes = lapply(regions_index,function(x) {
  
  cycle_aic_per_region = hhs_regions_all_stepfunc[[x]]$cycle_aic
  
  annual_cycle = cycle_aic_per_region$cycle[which.min(cycle_aic_per_region$aic)]
  
  if(annual_cycle != 52) {
    
    warning(paste('best cycle is not annual in region',x))
  }
  
  non_annual_cycle_range = cycle_aic_per_region[which(cycle_aic_per_region$cycle %in% c(27:31)),]
  non_annual_cycle = non_annual_cycle_range$cycle[which.min(non_annual_cycle_range$aic)]
  
  return(as.character(c(annual_cycle, non_annual_cycle)))
})

cycle_indexes = do.call(rbind, cycle_indexes)
colnames(cycle_indexes) = c('annual','non_annual')
cycle_indexes

par(mfcol = c(5,2),mar = c(2,2,2,2))

hhs_regions_cycles = lapply(regions_index, function(x) {
    
    plot(hhs_regions[[x]]$firstday,hhs_regions[[x]]$zeta_score,type = 'l',main = names(hhs_regions)[x])
    
    this_r_stepfunc = hhs_regions_all_stepfunc[[x]]$all_sim_ts
    
    lines(hhs_regions[[x]]$firstday,this_r_stepfunc[which(rownames(this_r_stepfunc) == cycle_indexes[x,'annual']),],col = 'red')
    lines(hhs_regions[[x]]$firstday,this_r_stepfunc[which(rownames(this_r_stepfunc) == cycle_indexes[x,'non_annual']),],col = 'blue')
    
    return(list(annual = this_r_stepfunc[which(rownames(this_r_stepfunc) == cycle_indexes[x,'annual']),],
                non_annual = this_r_stepfunc[which(rownames(this_r_stepfunc) == cycle_indexes[x,'non_annual']),]))
})

#### Create regression dataset ###
hhs_regression_data = lapply(regions_index, function(x) {
    
    hhs_weather = get_weather_data(weather_csv_location = paste0('../_DATA/Weather Data NASA/HHS regions NASA POWER Point/POWER_Point_Daily_US_HHS_r',x,'.csv'),
                                   rows_to_skip = 11,
                                   ncol = 5,
                                   plot_raw = F,
                                   plot_weekly = F,
                                   save = F)
    
    annual_cycle = hhs_regions_cycles[[x]]$annual
    non_annual_cycle = hhs_regions_cycles[[x]]$non_annual
    
    get_regression_data(weather_data = hhs_weather,
                        ili_zeta_df = hhs_regions[[x]],
                        annual_cycle = annual_cycle,
                        time_name = 'firstday',
                        response_name = 'zeta_score',
                        non_annual_cycle = non_annual_cycle,
                        save = F)
    
})

names(hhs_regression_data) = paste0('r',regions_index)
save(hhs_regression_data,file = 'Rdata/HHS_r1_r10_regression_dataset.Rdata')
load('Rdata/HHS_r1_r10_regression_dataset.Rdata')

### didn't use lognormal regression, because proportion of zeroes cannot be ignored(~10%) ###
# hhs_regression_data = lapply(1:length(hhs_regression_data),function(x) {
#     
#     hhs_regression_data[[x]]$ili_zeta = log(hhs_regression_data[[x]]$ili_zeta)
#     
#     return(hhs_regression_data[[x]])
# })

### remove NA ###
hhs_complete_reg_data = lapply(regions_index, function(x) {
    
    na.omit(hhs_regression_data[[x]])
    # sum(is.infinite(hhs_regression_data[[x]]$ili_zeta))
})

lapply(regions_index,function(x) {sum(is.infinite(hhs_regression_data[[x]]$ili_zeta))})
lapply(regions_index,function(x) {sum(is.na(hhs_regression_data[[x]]$ili_zeta))})


hhs_selected_model = lapply(regions_index, function(x) {
    
    get_normal_model(response_name = 'ili_zeta',
                     data = hhs_complete_reg_data[[x]][,-1])
})

## Normal regression leads to few negative values ##
lapply(regions_index,function(x) {sum(fitted(hhs_selected_model[[x]]) <0)})

i = 9
summary(hhs_selected_model[[i]])

### plot data with model ##
par(mfcol = c(5,2),mar = c(2,2,2,2)) 
lapply(regions_index, function(x){
    
    plot(hhs_complete_reg_data[[x]]$date,hhs_complete_reg_data[[x]]$ili_zeta,type = 'l',main = names(hhs_regions)[x])
    lines(hhs_complete_reg_data[[x]]$date,hhs_selected_model[[x]]$fitted.values,col = 'red')
})


lapply(regions_index, function(x) {
    
    plot(hhs_regions[[x]]$weighted_ILIpercent,type = 'l')
})

### predictor importance in the full model ###
hhs_pred_impo_full = lapply(regions_index, function(x) {
    
    get_pred_impo_full_model(model = hhs_selected_model[[x]],
                             response_name = 'ili_zeta',
                             family = 'gaussian',
                             data = hhs_complete_reg_data[[x]][-1],
                             barplot = F)
})



lapply(regions_index,function(x) {
    
    barplot(-hhs_pred_impo_full[[x]]$diff_aic,
            names.arg = hhs_pred_impo_full[[x]]$preds,
            main = names(hhs_regions)[x])
})

### relative importance ###
hhs_rela_impo = lapply(regions_index, function(x) {
    
  calc.relimp(hhs_selected_model[[x]])
})
save(hhs_rela_impo,file = 'Rdata/rela_impo_zeta_reg_HHS_regions.Rdata')

### barplot of rela impo ###

hhs_lmg = lapply(regions_index,function(x) {
    
    each_lmg = hhs_rela_impo[[x]]@lmg
    each_lmg = each_lmg[order(each_lmg, decreasing = T)]
    
    
    each_lmg
})

par(mfrow = c(5,2),mar = c(2,2,2,2))
lapply(regions_index,function(x) {
    
    barplot(hhs_lmg[[x]],main = names(hhs_regions)[x])
})


