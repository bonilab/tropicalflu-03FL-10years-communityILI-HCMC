rm(list = ls())
setwd("~/Dropbox/Influenza and Respiratory Disease in Vietnam/code")
load("Rdata/clean ILIPERCENT.RData")
source("functions/Detrend and Smooth.R")

#### Calculate zeta score for each clinic ###
ZETA_MAT = matrix(NA, nrow = nrow(PERCENT),ncol = ncol(PERCENT))

for(j in 1:ncol(PERCENT)) {
    
    ZETA_MAT[,j] = zeta_score(PERCENT[,j],window_width = 365,min_available_days = 150)
}

rownames(ZETA_MAT) = rownames(PERCENT)
colnames(ZETA_MAT) = colnames(PERCENT)

#### Get weekly averaged zeta score from each clinic ###
AVG_ZETA = c()
week_index = seq(1,nrow(ZETA_MAT),7)

for(j in 1:ncol(ZETA_MAT)) {
    
    one_clinic_daily_zeta = ZETA_MAT[,j]
    weekly_zeta = c()
    
    for(i in 1:length(week_index)) {
        
        if(i == length(week_index)){
            
            one_week_zeta = one_clinic_daily_zeta[week_index[i]:length(one_clinic_daily_zeta)]
            
        } else {
        
        one_week_zeta = one_clinic_daily_zeta[week_index[i]:(week_index[i+1]-1)]
        
        }
        
        weekly_zeta = c(weekly_zeta,mean(one_week_zeta,na.rm = T))
    }
    
    AVG_ZETA = cbind(AVG_ZETA,weekly_zeta)
}

rownames(AVG_ZETA) = week_index
colnames(AVG_ZETA) = colnames(ZETA_MAT)

### Plot weekly zeta with daily zeta ###
par(mfrow = c(7,5),mar = c(2,2,2,2))

for(j in 1:ncol(ZETA_MAT)){
    
    plot(ZETA_MAT[,j],type = "l",main = paste("clinic",colnames(ZETA_MAT)[j]))
    lines(x = week_index,y = AVG_ZETA[,j],col = "red")
    
}

### Generate average zeta across the weekly zeta from all the clinics ###
week_days = seq(as.Date("2010-01-01"),as.Date("2019-12-31"),"7 days")

zeta_weekly_avg_df = data.frame(date = week_days,
                                weekly_zeta = rowMeans(AVG_ZETA,na.rm = T))

save(zeta_weekly_avg_df,file = "Rdata/Weekly zeta from averaging daily zeta in a week.Rdata")





