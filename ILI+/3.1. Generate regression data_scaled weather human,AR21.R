########## Regression on overall ILI+ ############
rm(list = ls())
library(car)
setwd("~/Dropbox/Influenza and Respiratory Disease in Vietnam/code")
source('functions/Detrend and Smooth.R')
load('Rdata/Daily weather in HCMC from NASA.Rdata')
load('Rdata/Overall.ILIplus.21-day.aggregated.pcr.RData')

#### Processing regression data set ####
# 1. Add school term as categorical predictors
# 2. Scale weather predictors
# 3. Smooth the numerical predictors
# 4. Add lagged version of smoothed numerical predictors(weather, iliplus)
# 5. Log transformation on response variable if needed
# 6. Check and remove any NA Inf


##### 1.  Add school term ######
data_for_fit = data.frame(date = ili_plus_21d_df$date,
                          smth_7d_iliplus = ili_plus_21d_df$ili_plus_smth_7d)

data_for_fit$school = 1
data_for_fit$year = format(data_for_fit$date,"%Y")
years = unique(data_for_fit$year)

#### School is open from Aug 15 to June 1 every year
for(year in years){
    
    holiday = seq(as.Date(paste0(year,"-06-02")),as.Date(paste0(year,"-08-15")),"1 day")
    data_for_fit$school[data_for_fit$date %in% holiday] = 0
   
}

data_for_fit$year = NULL

###### 2. Scale weather predictors #####
scaled_weather = data.frame(date = weather$date,
                            apply(weather[,-1],2,scale,center = T,scale = T))

###### 3. Smooth weather predictors ###
data_for_fit = cbind(data_for_fit,
                     apply(scaled_weather[,-1],2,moving_average,
                           window_width = 7))
                         

####### 2. Smooth the numerical variables: predictors and response variable ####

# Smooth the entire weather time series at first #
smth_weather = data.frame(date = scaled_weather$date)

for(j in 2:ncol(scaled_weather)) {
    
    smth_weather = cbind(smth_weather, moving_average(scaled_weather[,j],7))
    
}
colnames(smth_weather)[2:ncol(smth_weather)] = colnames(scaled_weather)[2:ncol(scaled_weather)]

### Add smoothed weather in the data ###
weather_same_period = smth_weather[smth_weather$date %in% data_for_fit$date,]

data_for_fit = data.frame(data_for_fit,
                          weather_same_period[,-1])


###### 3. Add lagged predictors #########
#### 3.1 Generate lagged weather #####

# Lagged weather is from 7d to 21d

max_lag = 21
day_lags = seq(7,max_lag,7)

weather_to_lag = smth_weather

#### Since the weather data prior to the ili data is available, 
# we just index the earlier weather data #

for(day_lag in day_lags){
    
    weather_lagged_period = data_for_fit$date - day_lag
    
    lag_df = weather_to_lag[weather_to_lag$date %in% weather_lagged_period,-1]
    
    colnames(lag_df) = paste0(colnames(weather_to_lag)[2:ncol(weather_to_lag)],'_lag_',day_lag,'d')
    
    data_for_fit = cbind(data_for_fit,lag_df)
}



##### 3.2 Generate 21d-lagged iliplus #####

day_lag = 21

iliplus_to_lag = ili_plus_21d_df$ili_plus_smth_7d

data_for_fit = data_for_fit[(day_lag + 1):nrow(data_for_fit),]

data_for_fit$iliplus_lag_21d = iliplus_to_lag[1:nrow(data_for_fit)]


#### 4. Check and remove NA and Inf #######

if(any(is.na(data_for_fit))) {
    
    missing_ind = c()
    
    for(i in 1:nrow(data_for_fit)) {
        
        if(any(is.na(data_for_fit[i,]))) {
            
            missing_ind = c(missing_ind, i)
        }
        
    }
    
}


inf_ind = c()

for(j in 1:ncol(data_for_fit)) {
    
    if(any(is.infinite(data_for_fit[,j]))) {
        
        inf_ind = c(inf_ind,j)
    }
}

data_for_fit = na.omit(data_for_fit)

save(data_for_fit,file = 'Rdata/Data of overall iliplus_UV_WS_absent_scaled_weather.AR_21.human for gamma hurdle model.Rdata')

rm(list = setdiff(ls(),'data_for_fit'))







