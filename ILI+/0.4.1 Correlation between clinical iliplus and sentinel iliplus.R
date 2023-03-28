####### correlation between clinical iliplus and sentinel iliplus ######
rm(list = ls())
setwd("~/Dropbox/Influenza and Respiratory Disease in Vietnam/code")
source('functions/Data aggregation.R')
source('functions/Detrend and Smooth.R')
source('functions/iliplus data analysis.R')
load('Rdata/All.iliplus.weekly.zeta_iliperc.21-day.pcr.sentinel.southern.VN.Rdata')
load('Rdata/clean ILINUM.RData')
load('Rdata/clean PATIENTS.RData')


ili_pcr = read.csv('../_DATA/Influenza PCR Data/Fuhan New Data Assembly October 2020 and later/AllPCRandSubtypeResults.03FLand17FL.20220210.csv')
ili_pcr$SampleDate = as.Date(paste0(ili_pcr$SampleYear,'-',
                                    ili_pcr$SampleMonth,'-',
                                    ili_pcr$SampleDay))

### Correlate clinical iliplus and sentinel iliplus in the same period
sentinel_period = sentinel_iliplus_list[[1]]$firstday

clinical_period = seq(min(ili_pcr$SampleDate,na.rm = T),
                      max(ili_pcr$SampleDate,na.rm = T),
                      '1 day')

common_period = clinical_period[clinical_period %in% sentinel_period]

### common period should be daily
common_period = seq(common_period[1],common_period[length(common_period)] + 6 ,'1 day')

#### aggregate daily ilinum/patients to weekly ilinum/patients #######
date_in_matrix = seq(as.Date('2010-01-01'),as.Date('2019-12-31'),'1 day')
date_index = findInterval(common_period,date_in_matrix)

selected_ILINUM = ILINUM[date_index,]

selected_WEEKLY_ILINUM = aggregate_matrix_daily_to_weekly(MATRIX = selected_ILINUM,
                                                          start_date = common_period[1],
                                                          end_date = common_period[length(common_period)])
selected_PATIENT = PATIENT[date_index,]

selected_WEEKLY_PATIENT = aggregate_matrix_daily_to_weekly(MATRIX = selected_PATIENT,
                                                          start_date = common_period[1],
                                                          end_date = common_period[length(common_period)])

selected_WEEKLY_ILIPERC = selected_WEEKLY_ILINUM/selected_WEEKLY_PATIENT

#### average ILIPERC across all clinics ####
weekly_iliperc_df = data.frame(firstday = seq(common_period[1],common_period[length(common_period)],'7 days'),
                               weekly_iliperc = rowMeans(selected_WEEKLY_ILIPERC,na.rm = T))

### weekly zeta 

weekly_zeta_perc_df = data.frame(firstday = seq(common_period[1],common_period[length(common_period)],'7 days'),
                                 weekly_zeta = zeta_score(weekly_iliperc_df$weekly_iliperc,
                                                          window_width = 53,
                                                          min_available_days = 21))

plot(weekly_zeta_perc_df$weekly_zeta,type = 'l')


rm('ILINUM','PATIENT','selected_ILINUM','selected_PATIENT','selected_WEEKLY_ILINUM','selected_WEEKLY_PATIENT',
   'selected_WEEKLY_ILIPERC','sentinel_period','clinical_period','weekly_iliperc_df',
   'date_in_matrix','date_index')

####### Aggregate 21-day pcr in the common period ##########
selected_ili_pcr = ili_pcr[ili_pcr$SampleDate %in% common_period,]

date_index_list = indexing_start_to_end(start = common_period[1],
                                        end = common_period[length(common_period)],
                                        interval = "21 days",
                                        date_in_data = selected_ili_pcr$SampleDate)

selected_21d_pcr = aggregate_pcr(data = selected_ili_pcr,
                                 date_labels = date_index_list$date_labels,
                                 data_date_index = date_index_list$date_index_data)

####### Calculate iliplus from weekly zeta and 21-day pcr data ######
clinical_weekly_iliplus = data.frame()
weekly_zeta_perc_df$ind_21d = findInterval(weekly_zeta_perc_df$firstday,selected_21d_pcr$firstday)

for(i in 1:nrow(selected_21d_pcr)) {
    
    df_tmp = weekly_zeta_perc_df[weekly_zeta_perc_df$ind_21d == i,]
    
    iliplus = round(df_tmp$weekly_zeta*selected_21d_pcr$pos_perc_np_swab[i]/100,2)
    H1plus = round(df_tmp$weekly_zeta*selected_21d_pcr$pos_perc_H1_np_swab[i]/100,2)
    H3plus = round(df_tmp$weekly_zeta*selected_21d_pcr$pos_perc_H3_np_swab[i]/100,2)
    Bplus = round(df_tmp$weekly_zeta*selected_21d_pcr$pos_perc_B_np_swab[i]/100,2)
    
    clinical_weekly_iliplus = rbind(clinical_weekly_iliplus,
                                    data.frame(firstday = df_tmp$firstday,
                                               iliplus = iliplus,
                                               H1plus = H1plus,
                                               H3plus = H3plus,
                                               Bplus = Bplus))
    
    }

rm('date_index_list','df_tmp','ili_pcr','selected_ili_pcr','Bplus',"H1plus","H3plus",'iliplus',
   'selected_21d_pcr','weekly_zeta_perc_df')


###### Plot and compare clinical iliplus and sentinel iliplus #####

for(i in 1:length(sentinel_iliplus_list)) {
    
    plot_sentinel_iliplus_with_clinical_iliplus(sentinel_data = sentinel_iliplus_list[[i]])
}

### plot weekly data causes gaps between the bars
### convert weekly iliplus to daily iliplus for plotting

plot_daily_sentinel_iliplus_list = list()
common_sentinel_iliplus_list = list()

for(i in 1:length(sentinel_iliplus_list)) {
    
    common_sentinel_iliplus_list[[i]] = sentinel_iliplus_list[[i]][sentinel_iliplus_list[[i]]$firstday %in% common_period,]
    plot_daily_sentinel_iliplus_list[[i]] = weekly_to_daily_iliplus(common_sentinel_iliplus_list[[i]],
                                                                    daily_period_in_data = common_period)
}

plot_daily_clinical_iliplus = weekly_to_daily_iliplus(iliplus_data = clinical_weekly_iliplus,
                                                      daily_period_in_data = common_period)

### plot by year 
## clinical iliplus ##

p1 = plot_stacked_iliplus_by_year(plot_daily_clinical_iliplus,
                             Bplus_label = "Clinical ILI+ of influenza B",
                             H1plus_label = "Clinical ILI+ of influenza H1N1",
                             H3plus_label = "Clinical ILI+ of influenza H3N2",
                             tag = '')

p1
#### Hospital for Tropical Diseases in HCMC ####
p2 = plot_stacked_iliplus_by_year(plot_daily_sentinel_iliplus_list[[2]],
                             Bplus_label = "Sentinel ILI+ of influenza B",
                             H1plus_label = "Sentinel ILI+ of influenza H1N1",
                             H3plus_label = "Sentinel ILI+ of influenza H3N2")

p2

library(cowplot)
ps = plot_grid(p1,p2,align = 'hv',axis = 'tblr')
ggsave(filename = 'plots/SuppFig7.jpg',plot = ps,
       width = 10,height = 6,units = 'in',dpi = 300)

### Correlation ###
correlation_all_iliplus(iliplus_data_1 = clinical_weekly_iliplus,
                        iliplus_data_2 = common_sentinel_iliplus_list[[2]])


