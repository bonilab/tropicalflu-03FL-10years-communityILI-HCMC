rm(list = ls())
source("functions/Data aggregation.R")
library(ggplot2)
library(dplyr)
library(reshape2)

ili_pcr = read.csv("datasets/ILI/AllPCRandSubtypeResults.03FLand17FL.20220210.csv",
                   header = T,stringsAsFactors = F,
                   sep = ",")

ili_pcr$SampleDate = paste0(ili_pcr$SampleYear,"-",ili_pcr$SampleMonth,"-",ili_pcr$SampleDay)
ili_pcr$SampleDate = as.Date(ili_pcr$SampleDate,tz = "UTC")

#### Only choose the period before 2020-01-01 ###
ili_pcr = ili_pcr[ili_pcr$SampleDate < as.Date('2020-01-01'),]


start = min(ili_pcr$SampleDate,na.rm = T)
end = max(ili_pcr$SampleDate,na.rm = T)

### Replace -99 with NA
ili_pcr[ili_pcr == -99] = NA

### Check if there are patients with both flu info is missing
any(is.na(ili_pcr$FluA) & is.na(ili_pcr$FluB))

### any missing fluA
any(is.na(ili_pcr$FluA))

### any missing fluB
# FluB is missing after 2019-07-01
if(any(is.na(ili_pcr$FluB))) {
    
    ili_pcr$SampleDate[which(is.na(ili_pcr$FluB))]
}

## Check if there are missing subtype info for fluA patients 
any(is.na(ili_pcr$H1[ili_pcr$FluA == 1]) & is.na(ili_pcr$H3[ili_pcr$FluA == 1]))

#### Only choose the period before 2020-01-01 ###
# ili_pcr = ili_pcr[ili_pcr$SampleDate < as.Date('2020-01-01'),]

#### Descriptive statistics ####

### Total flu ###
sum(ili_pcr$FluA == 1 | ili_pcr$FluB == 1)
sum(ili_pcr$FluA == 1 | ili_pcr$FluB == 1)/2604

### A ###
sum(ili_pcr$FluA == 1,na.rm = T)

### H1N1 ### 
sum(ili_pcr$H1 == 1,na.rm = T)

### H3N2 ###
sum(ili_pcr$H3 == 1,na.rm = T)

### B ###
sum(ili_pcr$FluB == 1,na.rm = T)


#### Aggregate ili_pcr data given different time interval ####

## interval = 21 days
date_index_list = indexing_start_to_end(start = start,
                                        end = end,
                                        interval = "21 days",
                                        date_in_data = ili_pcr$SampleDate)

aggregated_pcr_21d = aggregate_pcr(data = ili_pcr,
                  date_labels = date_index_list$date_labels,
                  data_date_index = date_index_list$date_index_data)



### sanity check ###
# for(i in 2:nrow(aggregated_pcr_21d)) {
#     
#     date_diff = aggregated_pcr_21d$firstday[i] - aggregated_pcr_21d$firstday[i-1]
#     if(date_diff !=21) {
#         
#         paste(i)
#     }
#     
# }

### Exclude the under-reported period during Lunar New Year ###
### Few tests and 0 positive test on these days because of no reporting ###
### Put zeroes in these days ###
removed_dates = c(as.Date('2014-01-22'),as.Date('2015-02-04'),as.Date('2016-01-27'),
                  as.Date('2017-01-18'),as.Date('2018-01-31'),as.Date('2019-01-23'))

aggregated_pcr_21d[aggregated_pcr_21d$firstday %in% removed_dates,]
# aggregated_pcr_21d[aggregated_pcr_21d$firstday %in% removed_dates,c(3:12)] = NA

plotted_vars = c('pos_num_np_swab','pos_num_B_np_swab','pos_num_H1_np_swab','pos_num_H3_np_swab',
                 'pos_perc_np_swab','pos_perc_B_np_swab','pos_perc_H1_np_swab','pos_perc_H3_np_swab')

par(mfcol = c(4,2),mar = c(2,2,2,2))
for(var in plotted_vars) {
    
    plot(aggregated_pcr_21d$firstday,aggregated_pcr_21d[,var],type = 'l',
         lwd = 1.5,main = var,xlab = '',ylab = '',
         xaxt = 'n')
    axis(1,at = seq(as.Date('2012-01-01'),as.Date('2021-01-01'),'1 year'),
         labels = format(seq(as.Date('2012-01-01'),as.Date('2021-01-01'),'1 year'),'%Y'))
    abline(v = seq(as.Date('2012-01-01'),as.Date('2021-01-01'),'1 year'),
           lty = 'dashed',col = 'grey')
}


save(aggregated_pcr_21d,file = "Rdata/PCR_data_aggregated_by_21d.RData")

#### Descriptive Statistics ####
### Only consider the period before 2020-01-01
aggregated_pcr_21d = aggregated_pcr_21d[aggregated_pcr_21d$firstday < as.Date('2020-01-01'),]

### total swabs ###
sum(aggregated_pcr_21d$num_np_swab,na.rm = T)

### total flu positive ###
sum(aggregated_pcr_21d$pos_num_np_swab,na.rm = T)
sum(aggregated_pcr_21d$pos_num_np_swab,na.rm = T)/sum(aggregated_pcr_21d$num_np_swab,na.rm = T)

### total H1N1 ###
sum(aggregated_pcr_21d$pos_num_H1_np_swab,na.rm = T)
sum(aggregated_pcr_21d$pos_num_H1_np_swab,na.rm = T)/sum(aggregated_pcr_21d$num_np_swab,na.rm = T)

### total H3N2 ###
sum(aggregated_pcr_21d$pos_num_H3_np_swab,na.rm = T)
sum(aggregated_pcr_21d$pos_num_H3_np_swab,na.rm = T)/sum(aggregated_pcr_21d$num_np_swab,na.rm = T)

### total B ###
sum(aggregated_pcr_21d$pos_num_B_np_swab,na.rm = T)
sum(aggregated_pcr_21d$pos_num_B_np_swab,na.rm = T)/sum(aggregated_pcr_21d$num_np_swab,na.rm = T)


### total swabs by year ###
library(dplyr)

aggregated_pcr_21d$year = format(aggregated_pcr_21d$firstday,'%Y')
aggregated_pcr_21d

summary_pcr = aggregated_pcr_21d %>% 
    group_by(year) %>%
    summarise(pos_num = sum(num_np_swab,na.rm = T),
              pos_perc = round(sum(pos_num_np_swab,na.rm = T)/sum(num_np_swab,na.rm = T),3),
              pos_perc_H1 = round(sum(pos_num_H1_np_swab,na.rm = T)/sum(num_np_swab,na.rm = T),3),
              pos_perc_H3 = round(sum(pos_num_H3_np_swab,na.rm = T)/sum(num_np_swab,na.rm = T),3),
              pos_perc_B = round(sum(pos_num_B_np_swab,na.rm = T)/sum(num_np_swab,na.rm = T),3))
### remove 2021 ###
summary_pcr = summary_pcr[-which(summary_pcr$year == 2021),]

#### flu statistics per year ####
summary(summary_pcr$pos_num)
summary(summary_pcr$pos_perc)
summary(summary_pcr$pos_perc_H1)
summary(summary_pcr$pos_perc_H3)
summary(summary_pcr$pos_perc_B)

