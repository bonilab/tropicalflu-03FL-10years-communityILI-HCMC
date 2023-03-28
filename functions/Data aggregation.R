# Work flow: 
# 1. Create time interval 
# 2. Get the date label for each interval 
# 3. Get the index for the date in the data based on the date label
# 4. Divide the data into intervals according to the index number in the data
# 5. Aggregate the data by interval


### Index starts from the beginning of the time series, and end at the end of the time series ###

indexing_start_to_end = function(start,end,interval,date_in_data) {
    
    all_date_labels = seq(start,end,by = interval)
    
    if(max(date_in_data,na.rm = T) > max(all_date_labels)) {
        
        num_days = as.numeric(gsub('(*\\d) days','\\1',interval))
        all_date_labels = c(all_date_labels,max(all_date_labels) + num_days)
    }
    
    date_index_data = findInterval(date_in_data,all_date_labels)
    
    return(list(date_labels = all_date_labels,
                date_index_data = date_index_data))
}

### Data aggregation
aggregate_pcr = function(data, 
                         date_labels,
                         data_date_index) {
    
    ili_aggregated = data.frame()
    data$int_index = data_date_index
    
    for (i in 1:length(date_labels)) {
        
        df_tmp = data[which(data$int_index == i),]
        
        if(nrow(df_tmp) == 0) {
            
            ili_aggregated= rbind(ili_aggregated,
                                  data.frame(firstday = date_labels[i],
                                             num_np_swab = NA,
                                             pos_num_np_swab = NA,
                                             pos_num_A_np_swab = NA,
                                             pos_num_B_np_swab = NA,
                                             pos_num_H1_np_swab = NA,
                                             pos_num_H3_np_swab = NA,
                                             pos_perc_np_swab = NA,
                                             pos_perc_A_np_swab = NA,
                                             pos_perc_B_np_swab = NA,
                                             pos_perc_H1_np_swab = NA,
                                             pos_perc_H3_np_swab = NA))
            
        } else {
            
            total_counts = nrow(df_tmp)
            total_A = ifelse(all(is.na(df_tmp$FluA)),NA,sum(df_tmp$FluA == 1,na.rm = T))
            total_B = ifelse(all(is.na(df_tmp$FluB)),NA,sum(df_tmp$FluB == 1,na.rm = T))
            
            # NA in subtype data means not done, it is because the patient is positive for fluB or ORV, 
            # thus the number of H1/H3 positive patients should be 0.  
            total_H1 = ifelse(all(is.na(df_tmp$H1)),0,sum(df_tmp$H1 == 1,na.rm = T))
            total_H3 = ifelse(all(is.na(df_tmp$H3)),0,sum(df_tmp$H3 == 1,na.rm = T))
            

            ### This requires at least one flu data is available ###
            # All fluA info is available in HCMC clinics data
            positive_counts = sum(total_A,total_B,na.rm = T)
            
            ili = round(positive_counts/total_counts,digits = 4)*100
            
            ili_A = round(total_A/total_counts,digits = 4)*100
            ili_B = round(total_B/total_counts,digits = 4)*100
            
            ili_H1 = round(total_H1/total_counts,digits = 4)*100
            ili_H3 = round(total_H3/total_counts,digits = 4)*100
            
            ili_aggregated= rbind(ili_aggregated,
                                  data.frame(firstday = date_labels[i],
                                             num_np_swab = total_counts,
                                             pos_num_np_swab = positive_counts,
                                             pos_num_A_np_swab = total_A,
                                             pos_num_B_np_swab = total_B,
                                             pos_num_H1_np_swab = total_H1,
                                             pos_num_H3_np_swab = total_H3,
                                             pos_perc_np_swab = ili,
                                             pos_perc_A_np_swab = ili_A,
                                             pos_perc_B_np_swab = ili_B,
                                             pos_perc_H1_np_swab = ili_H1,
                                             pos_perc_H3_np_swab = ili_H3))
        }
        
    }
    
    data$int_index = NULL
    return(ili_aggregated)
}

aggregate_14d_annual_conserved = function(start_year,end_year,data) {
    
    whole_ts = c(start_year:end_year)
    aggregated_pcr = data.frame()
    
    for(year in whole_ts) {
        ## Specify the entire year
        whole_year = seq(as.Date(paste0(year,"-01-01")),as.Date(paste0(year,"-12-31")),"1 day")
        
        # Divide the year into 14-day intervals
        
        # Date labels: the firstday of every 14d interval
        date_labels = seq(whole_year[1],whole_year[length(whole_year)],"14 days")
        
        if(length(whole_year) == 366) {
            
            # To make sure the intervals end at 12-31, add one day in the interval 
            # that includes 02-29 in the leap year
            date_labels[6:length(date_labels)] = date_labels[6:length(date_labels)] + 1
        }
        
        # 365th day is in the last interval 27
        # Since we want to aggregate the exactly whole year, we include 365th day into 26th interval
        # Every year has 26 14-day intervals 
        date_labels = date_labels[1:26]
        
        
        # Find the interval of each day based on the date label
        whole_day_of_year = format(whole_year,"%j")
        ind_14d_in_one_year = findInterval(whole_year,date_labels)
        
        whole_year_df = data.frame(date = whole_year,
                                   index = ind_14d_in_one_year)
        
        # Use the whole year date index as a metric, find the interval of the sample date in the data
        pcr_this_year = data[data$SampleYear == year,]
        
        pcr_this_year$ind = sapply(1:nrow(pcr_this_year),function(x) {
            
            whole_year_df$index[which(whole_year_df$date == pcr_this_year$SampleDate[x])]
        })
        
        aggregated_ili_this_year = aggregate_pcr(data = pcr_this_year,
                                                 date_labels = date_labels,
                                                 data_date_index = pcr_this_year$ind)
        
        aggregated_pcr = rbind(aggregated_pcr,aggregated_ili_this_year)
    }
    
    return(aggregated_pcr)
    
}


aggregate_matrix_daily_to_weekly = function(MATRIX,start_date,end_date){
    
    ### Get week index for each day in the matrix
    week_indexes = seq(1,nrow(MATRIX),7)
    indexes = sapply(1:nrow(MATRIX),function(x){findInterval(x,week_indexes)})
    
    #### The remainder of 3652/7 is 5, thus the last weekly data is the sum of 5 days
    
    WEEKLY_DATA = c()
    
    for(i in 1:ncol(MATRIX)){
        
        one_clinic_daily = cbind(MATRIX[,i],indexes)
        one_clinic_weekly = c()
        
        for(index in 1:length(week_indexes)){
            
            one_clinic_one_week = sum(one_clinic_daily[which(one_clinic_daily[,"indexes"] == index),1],na.rm = T)
            one_clinic_weekly = rbind(one_clinic_weekly,one_clinic_one_week)
        }
        
        WEEKLY_DATA = cbind(WEEKLY_DATA,one_clinic_weekly)
    }
    
    rownames(WEEKLY_DATA) = paste("firstday",seq(start_date,end_date,"7 days"))
    colnames(WEEKLY_DATA) = colnames(MATRIX)
    
    return(WEEKLY_DATA)
}

aggregate_matrix_daily_to_any_length = function(MATRIX,start_date,end_date,window_length){
    
    ### Get day index in the window for each day in the matrix
    indexes = seq(1,nrow(MATRIX),window_length)
    indexes_in_mat = sapply(1:nrow(MATRIX),function(x){findInterval(x,indexes)})
    
    #### The last window might include fewer number of days #####
    
    AGG_DATA = c()
    
    for(i in 1:ncol(MATRIX)){
        
        one_clinic_daily = cbind(MATRIX[,i],indexes_in_mat)
        one_clinic_all_windows = c()
        
        for(index in 1:length(indexes)){
            
            one_clinic_one_window = sum(one_clinic_daily[which(one_clinic_daily[,"indexes_in_mat"] == index),1],na.rm = T)
            one_clinic_all_windows = rbind(one_clinic_one_window,one_clinic_all_windows)
        }
        
        AGG_DATA = cbind(AGG_DATA,one_clinic_all_windows)
    }
    
    window_interval = paste(window_length,'days')
    rownames(AGG_DATA) = paste("firstday",seq(start_date,end_date,window_interval))
    colnames(AGG_DATA) = colnames(MATRIX)
    
    return(AGG_DATA)
}
# aggregate_num = function(df,start,end,interval){
#     
#     t = seq(start,end,by = interval)
#     complete_int = findInterval(t,t)
#     df$int_index = findInterval(df$Date,t)
#     
#     ili_aggregated = data.frame()
#     
#     for (i in 1:length(complete_int)){
#         
#         df_tmp = df[which(df$int_index == complete_int[i]),]
#         
#         if(nrow(df_tmp) == 0){
#             
#             if(i + 1 <= length(complete_int)){
#                 ili_aggregated= rbind(ili_aggregated,
#                                       data.frame(firstday = t[i],
#                                                  interval = paste(t[i],"-",t[i + 1] - 1),
#                                                  ilinum = NA))
#             } else {
#                 
#                 ili_aggregated= rbind(ili_aggregated,
#                                       data.frame(firstday = t[i],
#                                                  interval = paste(t[i],"-",end),
#                                                  ilinum = NA))
#             }
#             
#         } else {
#             
#             ilinum = sum(df_tmp$ilinum,na.rm = T)
#             
#             if(i + 1 <= length(complete_int)){
#                 
#                 ili_aggregated= rbind(ili_aggregated,
#                                       data.frame(firstday = t[i],
#                                                  interval = paste(t[i],"-",t[i + 1] - 1),
#                                                  ilinum = ilinum))
#             } else {
#                 ili_aggregated= rbind(ili_aggregated,
#                                       data.frame(firstday = t[i],
#                                                  interval = paste(t[i],"-",end),
#                                                  ilinum = ilinum))
#                 
#             }
#         }
#     }
#     
#     df$int_index = NULL
#     return(ili_aggregated)
# }