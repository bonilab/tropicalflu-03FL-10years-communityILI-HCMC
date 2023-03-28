######### Compare with Sentinel data ############
rm(list = ls())
source('functions/Detrend and Smooth.R')
source('functions/Data aggregation.R')

######## Calculate ili+ from sentinel data ##########
sentinel_data = read.csv('datasets/Influenza PCR Data/Southern_Vn_Flu_Sentinel.csv')

# Sentinel data includes ILI patients, total patients, flu cases from 4 hospitals in southern VN
# location 14 15 are two hospitals in HCMC
# Where there are no visits, there is no valid data
# Where there are no tests, there is no valid data 

# Replace -99 to NA
sentinel_data[sentinel_data == -99] = NA

# Get weekly ILI percent 
sentinel_data$ILIperc = round(sentinel_data$Total_ILI/sentinel_data$Num_Patients,2)
all_locations = unique(sentinel_data$Location)
sentinel_list = list()
i = 1

for(location in all_locations) {
    
    sentinel_list[[i]] = sentinel_data[which(sentinel_data$Location == location),]
    i = i + 1
}

names(sentinel_list) = all_locations

# Add the firstday of the week in the data, the firstday is Sunday 
for(i in 1:length(sentinel_list)) {
    
    sentinel_list[[i]]$Firstday = NA
    for(j in 1:nrow(sentinel_list[[i]])) {
        
        if(j == 1) {
            
            sentinel_list[[i]]$Firstday[j] = as.Date('2006-01-01')
        } else {
            
            sentinel_list[[i]]$Firstday[j] = sentinel_list[[i]]$Firstday[j - 1] + 7
        }
    }
    
    sentinel_list[[i]]$Firstday = as.Date(sentinel_list[[i]]$Firstday,origin = "1970-01-01")
}

### All locations have the same period, select period of location 13
whole_period = sentinel_list[[1]]$Firstday

####### Get weekly zeta ########### 
sentinel_zeta_list = list()

for(i in 1:length(sentinel_list)) {
    
    Zeta = zeta_score(sentinel_list[[i]]$ILIperc,53,21)
    Zeta_num = zeta_score(sentinel_list[[i]]$Total_ILI,53,21)
    
    zeta_df = data.frame(firstday = whole_period,
                         weekly_zeta_perc = Zeta,
                         weekly_zeta_num = Zeta_num)
    
    sentinel_zeta_list[[i]] = zeta_df
}

names(sentinel_zeta_list) = all_locations

##### Plotting Zeta Score ######
par(mfcol = c(2,2),mar = c(2,2,2,2))
for(i in 1:length(sentinel_list)) {
    
    plot(x = sentinel_zeta_list[[i]]$firstday,
         y = sentinel_zeta_list[[i]]$weekly_zeta_perc,
         type = 'l',
         main = all_locations[i],
         ylim = c(min(sentinel_zeta_list[[i]][,c(2:3)],na.rm = T),
                  max(sentinel_zeta_list[[i]][,c(2:3)],na.rm = T)))
    
    lines(x = sentinel_zeta_list[[i]]$firstday,
          y = sentinel_zeta_list[[i]]$weekly_zeta_num,
          col = 'red')
    
}

# zeta of ILI num is more stable than zeta of ILI percent
save(sentinel_zeta_list,file = 'Rdata/weekly.zeta.sentinel.southern.VN.Rdata')
rm('Zeta','Zeta_num','location','i','j','sentinel_data','zeta_df')


### Get 21-day aggregated positive rate ###

sentinel_21d_aggregated_pcr_list = lapply(1:length(sentinel_list),function(x) {
    
    index = indexing_start_to_end(start = min(sentinel_list[[x]]$Firstday),
                                  end = max(sentinel_list[[x]]$Firstday),
                                  interval = '21 days',
                                  date_in_data = sentinel_list[[x]]$Firstday)
    
    aggregated_flu = data.frame()
    
    for(i in 1:length(index$date_labels)) {
        
        df_tmp = sentinel_list[[x]][which(findInterval(sentinel_list[[x]]$Firstday,index$date_labels) == i),]
        
        total_tests = ifelse(all(is.na(df_tmp$Total_tests)),NA,sum(df_tmp$Total_tests,na.rm = T))
        total_H1 = ifelse(all(is.na(df_tmp$H1_pos)),NA,sum(df_tmp$H1_pos,na.rm = T))
        total_H3 = ifelse(all(is.na(df_tmp$H3_pos)),NA,sum(df_tmp$H3_pos,na.rm = T))
        total_B = ifelse(all(is.na(df_tmp$B_pos)),NA,sum(df_tmp$B_pos,na.rm = T))
        
        if(all(is.na(c(total_H1,total_H3,total_B)))) {
            
            positive_counts = NA
            
        } else {
            
            positive_counts = sum(total_H1,total_H3,total_B,na.rm = T)
        }
        
        positive_perc = round(positive_counts/total_tests,2)
        H1_perc = round(total_H1/total_tests,2)
        H3_perc = round(total_H3/total_tests,2)
        B_perc = round(total_B/total_tests,2)
        
        aggregated_flu = rbind(aggregated_flu,
                               data.frame(firstday = df_tmp$Firstday[1],
                                          TestsNum = total_tests,
                                          fluNum = positive_counts,
                                          fluPerc = positive_perc,
                                          H1Num = total_H1,
                                          H1Perc = H1_perc,
                                          H3Num = total_H3,
                                          H3Perc = H3_perc,
                                          BNum = total_B,
                                          BPerc = B_perc))    
        
    }
    
    return(aggregated_flu)
    
})

# Plotting the variables in the raw data ##
vars = colnames(sentinel_list[[1]])[4:10]
years = unique(sentinel_list[[1]]$Year)

par(mfcol = c(7,4),mar = c(2,2,2,2))

for(i in 1:length(sentinel_list)) {
    
    for(var in vars) {
        
        plot(x = sentinel_list[[i]]$Firstday,y = sentinel_list[[i]][,var],type = "l",
             main = paste(var,'in location',all_locations[i]))
        abline(v = seq(as.Date(paste0(years[1],'-01-01')),as.Date(paste0(years[length(years)],'-01-01')),'1 year'),
               lty = 'dashed',col = 'grey')
    }
}

## Compare aggregated flu data with raw data 
all_vars = data.frame(raw_vars = c('Total_tests','H1_pos','H3_pos','B_pos'),
                      aggregated_vars = c('TestsNum','H1Num','H3Num','BNum'))

par(mfcol = c(4,4),mar = c(2,2,2,2))

for(i in 1:length(sentinel_list)) {
    
    for(j in 1:nrow(all_vars)) {
        
        plot(x = sentinel_list[[i]]$Firstday,y = sentinel_list[[i]][,all_vars[j,]$raw_vars],
             type = 'l',
             ylim = c(min(sentinel_list[[i]][,all_vars[j,]$raw_vars],na.rm = T),
                      max(sentinel_21d_aggregated_pcr_list[[i]][,all_vars[j,]$aggregated_vars],na.rm = T)),
             main = paste(all_vars[j,]$raw_vars,'in Location',all_locations[i]))
        lines(x = sentinel_21d_aggregated_pcr_list[[i]]$firstday,
              y = sentinel_21d_aggregated_pcr_list[[i]][,all_vars[j,]$aggregated_vars],
              col = 'red')
        abline(v = seq(as.Date(paste0(years[1],'-01-01')),as.Date(paste0(years[length(years)],'-01-01')),'1 year'),
               lty = 'dashed',col = 'grey')
        
    }
}

save(sentinel_21d_aggregated_pcr_list,file = 'Rdata/21-day.aggregated.pcr.from.sentinel.southern.VN.Rdata')
rm('all_vars','i','j','var','vars')


##### Calculate ILI+ from Zeta_ILIperc and 21-day aggregated pcr rate ######

sentinel_iliplus_list = list()

for(i in 1:length(sentinel_list)) {
    
    sentinel_zeta_list[[i]]$ind = findInterval(sentinel_zeta_list[[i]]$firstday,
                                          sentinel_21d_aggregated_pcr_list[[i]]$firstday)
    iliplus_df = data.frame()
    for(j in 1:nrow(sentinel_21d_aggregated_pcr_list[[i]])) {
        
        df_tmp = sentinel_zeta_list[[i]][sentinel_zeta_list[[i]]$ind == j,]
        iliplus = round(df_tmp$weekly_zeta_perc*sentinel_21d_aggregated_pcr_list[[i]]$fluPerc[j],2)
        H1plus = round(df_tmp$weekly_zeta_perc*sentinel_21d_aggregated_pcr_list[[i]]$H1Perc[j],2)
        H3plus = round(df_tmp$weekly_zeta_perc*sentinel_21d_aggregated_pcr_list[[i]]$H3Perc[j],2)
        Bplus = round(df_tmp$weekly_zeta_perc*sentinel_21d_aggregated_pcr_list[[i]]$BPerc[j],2)
        
        
        iliplus_df = rbind(iliplus_df,
                           data.frame(firstday = df_tmp$firstday,
                                      iliplus = iliplus,
                                      H1plus = H1plus,
                                      H3plus = H3plus,
                                      Bplus = Bplus))
    }
    
    sentinel_iliplus_list[[i]] = iliplus_df
}

names(sentinel_iliplus_list) = all_locations

vars = colnames(sentinel_iliplus_list[[1]])[2:5]

par(mfcol = c(4,4),mar = c(2,2,2,2))
for(i in 1:length(sentinel_iliplus_list)) {
    
    for(var in vars) {
        
        plot(x = sentinel_iliplus_list[[i]]$firstday,
             y = sentinel_iliplus_list[[i]][,var],
             type = 'l',
             main = paste(var,'from Location',all_locations[i]))
        abline(v = seq(as.Date(paste0(years[1],'-01-01')),as.Date(paste0(years[length(years)],'-01-01')),'1 year'),
               lty = 'dashed',col = 'grey')
    }
}

save(sentinel_iliplus_list,file = 'Rdata/All.iliplus.weekly.zeta_iliperc.21-day.pcr.sentinel.southern.VN.Rdata')
rm('df_tmp','iliplus_df','sentinel_iliplus_list','Bplus','H1plus','H3plus','iliplus')


####### Calculate ILI+ from Zeta_ilinum and 21-day aggregated pcr rate
sentinel_iliplus_list = list()

for(i in 1:length(sentinel_list)) {
    
    sentinel_zeta_list[[i]]$ind = findInterval(sentinel_zeta_list[[i]]$firstday,
                                               sentinel_21d_aggregated_pcr_list[[i]]$firstday)
    iliplus_df = data.frame()
    for(j in 1:nrow(sentinel_21d_aggregated_pcr_list[[i]])) {
        
        df_tmp = sentinel_zeta_list[[i]][sentinel_zeta_list[[i]]$ind == j,]
        iliplus = round(df_tmp$weekly_zeta_num*sentinel_21d_aggregated_pcr_list[[i]]$fluPerc[j],2)
        H1plus = round(df_tmp$weekly_zeta_num*sentinel_21d_aggregated_pcr_list[[i]]$H1Perc[j],2)
        H3plus = round(df_tmp$weekly_zeta_num*sentinel_21d_aggregated_pcr_list[[i]]$H3Perc[j],2)
        Bplus = round(df_tmp$weekly_zeta_num*sentinel_21d_aggregated_pcr_list[[i]]$BPerc[j],2)
        
        
        iliplus_df = rbind(iliplus_df,
                           data.frame(firstday = df_tmp$firstday,
                                      iliplus = iliplus,
                                      H1plus = H1plus,
                                      H3plus = H3plus,
                                      Bplus = Bplus))
    }
    
    sentinel_iliplus_list[[i]] = iliplus_df
}

vars = colnames(sentinel_iliplus_list[[1]])[2:5]

par(mfcol = c(4,4),mar = c(2,2,2,2))
for(i in 1:length(sentinel_iliplus_list)) {
    
    for(var in vars) {
        
        plot(x = sentinel_iliplus_list[[i]]$firstday,
             y = sentinel_iliplus_list[[i]][,var],
             type = 'l',
             main = paste(var,'from Location',all_locations[i]))
        abline(v = seq(as.Date(paste0(years[1],'-01-01')),as.Date(paste0(years[length(years)],'-01-01')),'1 year'),
               lty = 'dashed',col = 'grey')
    }
}

save(sentinel_iliplus_list,file = 'Rdata/All.iliplus.weekly.zeta_ilinum.21-day.pcr.sentinel.southern.VN.Rdata')
rm(list = ls())
