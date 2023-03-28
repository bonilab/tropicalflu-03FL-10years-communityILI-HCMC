
zeta_score = function(data,window_width = 365,min_available_days = 182){
    
    half_window_width = (window_width - 1)/2
    
    zeta_series = rep(NA, length(data))
    
    for (i in 1:length(data) ) {
        
        # Define window left and window right
        wl = i - half_window_width
        wr = i + half_window_width
        
        # at the beginning and the end of the ilipercent series,window is truncated
        # and i is not at the center
        if ( wl < 0 ) {
            
            wl = 1
        }
        if ( wr > length(data)) {
            
            wr = length(data)
        }
        
        window_temp = data[wl:wr]
        
        # are window-temp and window_temp empty?
        
        n_available_days = sum(!is.na(window_temp))
        if(n_available_days >= min_available_days && n_available_days >0) {
            
            zeta_series[i] = data[i]/mean(window_temp,na.rm = T)
        }
    }
    return(zeta_series)
}

z_score = function(data,window_width = 365,minAvailDays = 182){
    
    half_window_width = (window_width - 1)/2
    
    z_series = rep(NA, length(data))
    
    for (i in 1:length(data) ) {
        
        # Define window left and window right
        wl = i - half_window_width
        wr = i + half_window_width
        
        # at the beginning and the end of the ilipercent series,window is truncated
        # and i is not at the center
        if ( wl < 0 ) {
            
            wl = 1
        }
        if ( wr > length(data)) {
            
            wr = length(data)
        }
        
        window_temp = data[wl:wr]
        
        nAvailDays = sum(!is.na(window_temp))
        if(nAvailDays >= minAvailDays && nAvailDays >0 ) {
            
            z_series[i] = (data[i] - mean(window_temp,na.rm = T))/sd(window_temp,na.rm = T)
        }
    }
    return(z_series)
    }


moving_average = function(data,window_width = 15){
    
    if(window_width%%2 == 0) {
        
        stop("window width must be odd number")
    }
    
    half_window_width = (window_width - 1)/2
    
    average_z = c()
    
    for (center_in_window in 1:length(data)){
        
        wl = center_in_window - half_window_width
        wr = center_in_window + half_window_width
        
        if ( wl < 1 ) {
            
            wl = 1
        }
        
        if ( wr > length(data)) {
            
            wr = length(data)
        }
        
        # extract the data from ili number within the window
        window_temp = data[wl : wr] 
        
        
        # exclude -99 from the window
        if (all(is.na(window_temp))){
            
            average_z[center_in_window] = NA
            
            
        } else if (any(is.na(window_temp))) {
            
            # moving average
            average_z[center_in_window] = mean(window_temp[-which(is.na(window_temp))])
            
        } else {
            
            # moving average
            average_z[center_in_window] = mean(window_temp)
            
        }
        
        
    }
    return(average_z)
}


