########################################################################################
#                                                                                      #
# We calculate the number of days between two dates in three steps                     #   
#                                                                                      #
# 1. Define if startdate and the date are in the same year                             # 
# 1.1 Yes: daynum = days in month gap + days in date gap -> sameyear                   #
# 1.2 No: daynum = days in the year of stardate + days in the complete year sequence   #
#                  (1/1 of next year to 1/1 of the end year) + days in the end year    #
#                                                                                      #
# begin.  How many days left in the start year until 01/01 in next year                #
# middle. How many years from 01/01 in next year of the start year and 01/01 of        #
#         the end year                                                                 #
# end.    How many days from 01/01 in the end year until the end date                  #
#                                                                                      #
########################################################################################      

# Calculate the days from start date and the date in the same year
get_numdays_if_in_same_year <-  function(date,startdate){
        
        # Extract info from start date
        start_y = as.numeric(format(startdate, format = "%Y"))
        start_m = as.numeric(format(startdate, format = "%m"))
        start_d = as.numeric(format(startdate, format = "%d"))
        
        # Extract info from date
        date_y = as.numeric(format(date, format = "%Y"))
        date_m = as.numeric(format(date, format = "%m"))
        date_d = as.numeric(format(date, format = "%d"))
        
        # Define the number of days between 01/01 and the 1st day in each month
        d_thru_m_normal = c(0,31,59,90,120,151,181,212,243,273,304,334) 
        d_thru_m_leap = c(0,31,60,91,121,152,182,213,244,274,305,335)
        
        
        # Test if start date is in leap year
        if (start_y %% 4 == 0){
                d_thru_m = d_thru_m_leap    
        } else {
                d_thru_m = d_thru_m_normal    
        }
        
        day_diff = d_thru_m[date_m] - d_thru_m[start_m] + date_d - start_d 
        return(day_diff)
        
}
# Calculate the days from the startdate to 01/01 in next year
get_numdays_to_end_of_year <- function(date,startdate){
        
        # Extract info from start date
        start_y = as.numeric(format(startdate, format = "%Y"))
        start_m = as.numeric(format(startdate, format = "%m"))
        start_d = as.numeric(format(startdate, format = "%d"))
        
        # Extract info from date
        date_y = as.numeric(format(date, format = "%Y"))
        date_m = as.numeric(format(date, format = "%m"))
        date_d = as.numeric(format(date, format = "%d"))
        
        # Define the number of days for each month in leap and common years
        d_in_m_leap = c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
        d_in_m_normal = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
        
        # How many days from this month to 01/01 next year
        day_month_gap = 0
        
        # Test if start date is in leap year
        if (start_y %% 4 == 0){
                day_month_gap = sum(d_in_m_leap[start_m:12])       
        } else {
                day_month_gap = sum(d_in_m_normal[start_m:12]) 
        }
        
        # 
        day_begin_year = day_month_gap - start_d + 1
        return(day_begin_year)
}

# Calculate the days from 01/01 next year to 01/01 of the end year
get_num_days_in_middle_section <- function(start_y_next, date, startdate){
        
        # Define the year of start date and date
        start_y = as.numeric(format(startdate, format = "%Y"))
        date_y = as.numeric(format(date, format = "%Y"))
        
        # Define complete year sequence
        year_sequence = c(start_y_next : date_y)
        
        day_middle = 0
        day_temp = 0
        d_in_year = 0
        
        # Calculate the days in the year sequence
        for ( i in year_sequence){
                
                # Test if the sequence contains leap year
                if (i %% 4 == 0) {
                        d_in_year = 366
                } else {
                        d_in_year = 365
                }
                day_temp = sum(day_temp) + d_in_year 
        }
        day_middle = day_temp - d_in_year
        return(day_middle)
}


# Calculate the days from the 01/01 until the end date in the last year
get_numdays_from_beg_of_year <- function(date){
        
        # Extract date information
        date_y = as.numeric(format(date, format = "%Y"))
        date_m = as.numeric(format(date, format = "%m"))
        date_d = as.numeric(format(date, format = "%d"))
        
        # Define the number of days between 01/01 and the 1st day in each month
        d_thru_m_normal = c(0,31,59,90,120,151,181,212,243,273,304,334) 
        d_thru_m_leap = c(0,31,60,91,121,152,182,213,244,274,305,335)
        
        d_thru_m = 0
        
        # Test if the date in in leap year
        if (date_y %% 4 == 0){
                d_thru_m = d_thru_m_leap
        } else {
                d_thru_m = d_thru_m_normal
        }
        day_end = d_thru_m[date_m] + date_d - 1
        return(day_end)
}

daynum = function(date, startdate){
        
        start_y = as.numeric(format(startdate, format = "%Y"))
        date_y = as.numeric(format(date, format = "%Y"))
        
        if (isTRUE(date < startdate)){
                stop("date should be larger than startdate")
        }
        # This is how many days that we have counted
        total_days = 0
        
        # Test if startdate and date are in the same year
        if (isTRUE(start_y == date_y)) {
                total_days = get_numdays_if_in_same_year(date,startdate)
                return(total_days)
        } else {
        
        # How many days from the startdate to 01/01 next year?
        day_begin_year = get_numdays_to_end_of_year(date,startdate)
        start_y_next = start_y + 1
        
        # How many days in the complete-year sequence?
        day_middle = get_num_days_in_middle_section(start_y_next, date,startdate)
        
        # How many days from the 1/1 to the date in the last year?
        day_end = get_numdays_from_beg_of_year(date)
        
        total_days = day_begin_year + day_middle + day_end
        return(total_days)
        }
}






















