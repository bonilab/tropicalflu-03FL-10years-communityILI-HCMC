######## Correlation between zeta of ILIPERC in 21-day and pcr rate ############
rm(list = ls())
#### aggregate ILI patients and total patients in every 21 days ####
load('Rdata/clean ILINUM.RData')
load('Rdata/clean PATIENTS.RData')
source('functions/Data aggregation.R')

#### 21-day aggregated total patient matrix ####

#### select the clinical data in the same time period as the 21-day pcr ###
load('Rdata/PCR data aggregated by 21d_holiday_zeroed_20220210.RData')

clinical_data_ts = seq(as.Date('2010-01-01'),as.Date('2019-12-31'),'1 day')
pcr_data_ts = seq(min(aggregated_pcr_21d$firstday),max(aggregated_pcr_21d$firstday), '1 day')

common_date = clinical_data_ts[clinical_data_ts %in% pcr_data_ts]
common_date_ind_in_mat = which(clinical_data_ts %in% common_date)

SELECTED_PATIENT = PATIENT[common_date_ind_in_mat,]
SELECTED_ILINUM = ILINUM[common_date_ind_in_mat,]

PATIENTS_21d = aggregate_matrix_daily_to_any_length(MATRIX = SELECTED_PATIENT,
                                                    start_date = min(common_date),
                                                    end_date = max(common_date),
                                                    window_length = 21)
ILINUM_21d = aggregate_matrix_daily_to_any_length(MATRIX = SELECTED_ILINUM,
                                                  start_date = min(common_date),
                                                  end_date = max(common_date),
                                                  window_length = 21)

PERCENT_21d = ILINUM_21d/PATIENTS_21d

par(mfcol = c(7,5),mar = c(2,2,2,2))
for(i in 1:ncol(PERCENT_21d)) {
    
    if(all(is.na(PERCENT_21d[,i]))) {
        
        next()
    } else {
    plot(PERCENT_21d[,i],type = 'l')
    }
}

### Zeta matrix for each clinic ###
source('functions/Detrend and Smooth.R')

ZETA_21d = apply(PERCENT_21d,2,zeta_score,
                 window_width = 17,
                 min_available_days = 7)

par(mfcol = c(7,5),mar = c(2,2,2,2))

for(i in 1:ncol(ZETA_21d)) {
    
    if(all(is.na(ZETA_21d[,i]))) {
        
        next()
    } else {
        plot(ZETA_21d[,i],type = 'l')
    }
}

zeta_21d_avg_clinic = rowMeans(ZETA_21d,na.rm = T)
zeta_21d_avg_clinic_df = data.frame(firstday = seq(min(common_date),max(common_date),'21 days'),
                                    zeta = zeta_21d_avg_clinic)

#### compare zeta-21 day with 21-day pcr ###
selected_pcr = aggregated_pcr_21d[aggregated_pcr_21d$firstday %in% zeta_21d_avg_clinic_df$firstday,]

scaled_selected_pcr = scale(selected_pcr$pos_perc_np_swab)

cor.test(zeta_21d_avg_clinic_df$zeta,scaled_selected_pcr,
         method = 'pearson')

#### compare 21-day ILIpercent with 21-day pcr ###
ilipercent_21d_avg_clinic_df = data.frame(firstday = seq(min(common_date),max(common_date),'21 days'),
                                          ilipercent = rowMeans(PERCENT_21d,na.rm = T))
cor.test(ilipercent_21d_avg_clinic_df$ilipercent,
         selected_pcr$pos_perc_np_swab,
         method = 'pearson')

### compare 21-day ILIpercent with fluA perc ###
cor.test(ilipercent_21d_avg_clinic_df$ilipercent,
         selected_pcr$pos_perc_A_np_swab,
         mmethod = 'pearson')


### compare 21-day ILIpercent with H1N1 perc ###
cor.test(ilipercent_21d_avg_clinic_df$ilipercent,
         selected_pcr$pos_perc_H1_np_swab,
         mmethod = 'pearson')


### compare 21-day ILIpercent with H3N2 perc ###
cor.test(ilipercent_21d_avg_clinic_df$ilipercent,
         selected_pcr$pos_perc_H3_np_swab,
         mmethod = 'pearson')


### compare 21-day ILIpercent with B perc ###
cor.test(ilipercent_21d_avg_clinic_df$ilipercent,
         selected_pcr$pos_perc_B_np_swab,
         mmethod = 'pearson')
