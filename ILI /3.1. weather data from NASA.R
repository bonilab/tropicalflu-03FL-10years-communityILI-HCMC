rm(list = ls())
setwd("~/Dropbox/Influenza and Respiratory Disease in Vietnam/code")
weather = read.csv("datasets/POWER_Point_Daily_20100101_20191231_010d7714N_106d6895E_LST_VN.csv",
                   header = T,skip = 22)

weather$date = as.Date(paste(weather$YEAR,weather$DOY,sep = "-"),format = "%Y-%j")

############ Sanity Check ##########################
complete_ts = seq(as.Date("2010-01-01"),as.Date("2019-12-31"),"1 day")
nrow(weather) == length(complete_ts)

#### RH is measured given the temperature, so it is better to use AH to represent humidity

weather$AH = 0.611*exp(17.502*weather$T2M/(240.97+weather$T2M))*18.02*weather$RH2M/(8.31*(273.15+weather$T2M))

### Dew point is calculated from RH and temperature, remove it because of redundancy
### Specific humidity is removed, only use AH to represent humidity

### Plot each weather variable ###
weather_vars = c("T2M","AH","PRECTOTCORR","ALLSKY_SFC_UV_INDEX",
                 "WS2M")

par(mfrow = c(3,2),mar = c(2,2,2,2))

for(weather_var in weather_vars){
    
    plot(weather$date,weather[,weather_var],type = "l",main = weather_var)
}


weather = data.frame(date = weather$date,
                     temperature_2m = weather$T2M,
                     AH = weather$AH,
                     RF = weather$PRECTOTCORR)

### Check missing values in climate data ####

if(any(weather == -999)) {
    
    missing_ind = c()
    for(i in 1:nrow(weather)) {
        
        if(any(weather[i,] == -999)) {

            missing_ind = c(missing_ind, i)
        }
    }
    
}

#### Change -999 to NA in the missing data ###
# weather$UV[missing_ind] = NA
# weather[missing_ind,]

save(weather, file = "Rdata/Daily weather in HCMC from NASA.Rdata")

