####### Generate regression dataset for zeta ####################
rm(list = ls())
source('functions/Detrend and Smooth.R')

load('Rdata/Daily weather in HCMC from NASA.Rdata')
load('Rdata/7d.smth.zeta.ilipercent.Rdata')

##### 1. Scale daily weather ########
scaled_weather = data.frame(date = weather$date,
                            apply(weather[,-1],2,scale,center = T,scale = T))

##### 2. Smooth weather and zeta #####
data_for_fit = data.frame(date = weather$date,
                          zeta = zeta_perc_avg_df$smoothed_zeta_score,
                          apply(scaled_weather[,-1],2,moving_average,
                                window_width = 7))

##### 3. Add school term ######
# lunar = c(as.Date("2010-02-14"),as.Date("2011-02-03"),as.Date("2012-01-23"),
#           as.Date("2013-02-10"),as.Date("2014-01-31"),as.Date("2015-02-19"),
#           as.Date("2016-02-08"),as.Date("2017-01-28"),as.Date("2018-02-16"),
#           as.Date("2019-02-05"))

data_for_fit$school = 1
data_for_fit$year = format(data_for_fit$date,"%Y")
years = unique(data_for_fit$year)

#### School is open from Aug 15 to June 1 every year
### Lunar New Year starts from 7 days before New Year and lasts 6 days
for(year in years){
    
    holiday = seq(as.Date(paste0(year,"-06-02")),as.Date(paste0(year,"-08-15")),"1 day")
    # newyear_start = lunar[which(format(lunar,"%Y") == year)]
    # newyear = seq(newyear_start - 7,newyear_start + 6,"1 day")
    data_for_fit$school[data_for_fit$date %in% holiday] = 0
    # data_for_fit$holiday[data_for_fit$date %in% newyear] = 1
}

data_for_fit$year = NULL

##### 4. Add lagged version of weather ####

### lagged weather ###
max_lag = 21
day_lags = seq(7,max_lag,7)

weather_to_lag = data_for_fit[,which(colnames(data_for_fit) %in% c('temperature_2m','AH','RF'))]

new_dataset = data_for_fit[(max_lag + 1):nrow(data_for_fit),]

for(day_lag in day_lags){
    
    lagged_start = 1 + max_lag - day_lag
    lagged_end = lagged_start + nrow(new_dataset) -1
    
    lag_df = weather_to_lag[lagged_start:lagged_end,]
    
    colnames(lag_df) = paste0(colnames(weather_to_lag),'_lag_',day_lag,'d')
    
    new_dataset = cbind(new_dataset,lag_df)
}

###### 5. Log transform the zeta ########
zeta_lagged_period = new_dataset$date - 7
data_for_fit$logzeta = log(data_for_fit$zeta)

new_dataset$zeta_lag7d = data_for_fit$logzeta[data_for_fit$date %in% zeta_lagged_period]
new_dataset$logzeta = log(new_dataset$zeta)

## remove any NA/Inf ##
if(any(is.na(new_dataset$logzeta)) | any(is.infinite(new_dataset$logzeta))) {
    
    new_dataset = new_dataset[-which(is.na(new_dataset$logzeta) | is.infinite(new_dataset$logzeta))]
}

new_dataset$zeta = NULL
save(new_dataset,file = 'Rdata/Dataset_zeta_minAvailDay_150_scaled_weather_human_AR7.RData')
