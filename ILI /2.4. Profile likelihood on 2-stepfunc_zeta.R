######### Profile likelihood on stepfunction on zeta ###########################
rm(list = ls())
load('Rdata/7d.smth.zeta.ilipercent.Rdata')
source('functions/profile likelihood function on 2-stepfunction.R')

#### Profile on startpoint, bp1 #######
bp1_candidates = seq(35,75,10)

profile_llk_bp1 = list()
i = 1

for(bp1 in bp1_candidates) {
    
    profile_llk_bp1[[i]] = profile_llk_2_stepfunc(data = zeta_perc_avg_df$smoothed_zeta_score,
                                                  cycle = 365,
                                                  param_to_profile = 'bp1',
                                                  profile_param_value = bp1,
                                                  bp1_lower_limit = 50,
                                                  bp1_upper_limit = 60,
                                                  bp2_lower_limit = 130,
                                                  bp2_upper_limit = 140,
                                                  val1_lower_limit = 0.9,
                                                  val1_upper_limit = 1,
                                                  val2_lower_limit = 0.9,
                                                  val2_upper_limit = 1,
                                                  repeats = 100)
    i = i + 1
}
save(profile_llk_bp1,file = 'Rdata/Profile on bp1 2-stepfunc_zeta_sigma_0.1_cyc365.Rdata')

bp1_all_llk = sapply(1:length(profile_llk_bp1),function(x) {profile_llk_bp1[[x]]$llk})

#### Looks good! ###
plot(bp1_candidates,bp1_all_llk,type = 'l')

profile_llk_bp1
#### optimal likelihood - 1.92 ####
abline(h = max(bp1_all_llk) - 1.92,col = 'red')
abline(v = 53.63,col = 'grey')
abline(v = 55.4,col = 'grey')


#### Profile on endpoint, bp2 #######
bp2_candidates = seq(116,156,10)

profile_llk_bp2 = list()
i = 1

for(bp2 in bp2_candidates) {
    
    profile_llk_bp2[[i]] = profile_llk_2_stepfunc(data = zeta_perc_avg_df$smoothed_zeta_score,
                                                  cycle = 365,
                                                  param_to_profile = 'bp2',
                                                  profile_param_value = bp2,
                                                  bp1_lower_limit = 50,
                                                  bp1_upper_limit = 60,
                                                  bp2_lower_limit = 130,
                                                  bp2_upper_limit = 140,
                                                  val1_lower_limit = 0.9,
                                                  val1_upper_limit = 1,
                                                  val2_lower_limit = 0.9,
                                                  val2_upper_limit = 1,
                                                  repeats = 100)
    i = i + 1
}
save(profile_llk_bp2,file = 'Rdata/Profile on bp2 2-stepfunc_zeta_sigma_0.1_cyc365.Rdata')
profile_llk_bp2
bp2_all_llk = sapply(1:length(profile_llk_bp2),function(x) {profile_llk_bp2[[x]]$llk})

#### Looks good! ###
plot(bp2_candidates,bp2_all_llk,type = 'l')

#### optimal likelihood - 1.92 ####
abline(h = max(bp2_all_llk) - 1.92,col = 'red')
abline(v = 133.2,col = 'grey')
abline(v = 138.8,col = 'grey')


####### Profile on val1 ########
val1_cand = seq(0.81,1.21,0.1)

profile_llk_val1 = list()
i = 1

for(val1 in val1_cand) {
    
    profile_llk_val1[[i]] = profile_llk_2_stepfunc(data = zeta_perc_avg_df$smoothed_zeta_score,
                                                  cycle = 365,
                                                  param_to_profile = 'val1',
                                                  profile_param_value = val1,
                                                  bp1_lower_limit = 50,
                                                  bp1_upper_limit = 60,
                                                  bp2_lower_limit = 130,
                                                  bp2_upper_limit = 140,
                                                  val1_lower_limit = 0.9,
                                                  val1_upper_limit = 1,
                                                  val2_lower_limit = 0.9,
                                                  val2_upper_limit = 1,
                                                  repeats = 100)
    i = i + 1
}
save(profile_llk_val1,file = 'Rdata/Profile on val1 2-stepfunc_zeta_sigma_0.1_cyc365.Rdata')
profile_llk_val1

val1_all_llk = sapply(1:length(profile_llk_val1),function(x) {profile_llk_val1[[x]]$llk})

plot(val1_cand,val1_all_llk,type = 'l')

abline(h = max(val1_all_llk) - 1.92,col = 'red')
abline(v = 1.0095,col = 'grey')
abline(v = 1.0103,col = 'grey')


####### Profile on val2 ########
val2_cand = seq(0.72,1.12,0.1)

profile_llk_val2 = list()
i = 1

for(val2 in val2_cand) {
    
    profile_llk_val2[[i]] = profile_llk_2_stepfunc(data = zeta_perc_avg_df$smoothed_zeta_score,
                                                   cycle = 365,
                                                   param_to_profile = 'val2',
                                                   profile_param_value = val2,
                                                   bp1_lower_limit = 50,
                                                   bp1_upper_limit = 60,
                                                   bp2_lower_limit = 130,
                                                   bp2_upper_limit = 140,
                                                   val1_lower_limit = 0.9,
                                                   val1_upper_limit = 1,
                                                   val2_lower_limit = 0.9,
                                                   val2_upper_limit = 1,
                                                   repeats = 100)
    i = i + 1
}
save(profile_llk_val2,file = 'Rdata/Profile on val2 2-stepfunc_zeta_sigma_0.1_cyc365.Rdata')
profile_llk_val2

val2_all_llk = sapply(1:length(profile_llk_val2),function(x) {profile_llk_val2[[x]]$llk})

plot(val2_cand,val2_all_llk,type = 'l')

abline(h = max(val2_all_llk) - 1.92,col = 'red')
abline(v = 0.9193,col = 'grey')
abline(v = 0.9206,col = 'grey')



############## Cycle = 210 ###############
######### Profile likelihood on stepfunction on zeta ###########################
rm(list = ls())
load('Rdata/7d.smth.zeta.ilipercent.Rdata')
source('functions/profile likelihood function on 2-stepfunction.R')

#### Profile on startpoint, bp1 #######
bp1_candidates = seq(52,92,10)

profile_llk_cyc210_bp1 = list()
i = 1

for(bp1 in bp1_candidates) {
    
    profile_llk_cyc210_bp1[[i]] = profile_llk_2_stepfunc(data = zeta_perc_avg_df$smoothed_zeta_score,
                                                         cycle = 210,
                                                         param_to_profile = 'bp1',
                                                         profile_param_value = bp1,
                                                         bp1_lower_limit = 70,
                                                         bp1_upper_limit = 80,
                                                         bp2_lower_limit = 170,
                                                         bp2_upper_limit = 180,
                                                         val1_lower_limit = 0.9,
                                                         val1_upper_limit = 1,
                                                         val2_lower_limit = 0.9,
                                                         val2_upper_limit = 1,
                                                         repeats = 100)
    i = i + 1
}
save(profile_llk_cyc210_bp1,file = 'Rdata/Profile on bp1 2-stepfunc_zeta_sigma_0.1_cyc210.Rdata')

bp1_all_llk = sapply(1:length(profile_llk_cyc210_bp1),function(x) {profile_llk_cyc210_bp1[[x]]$llk})

profile_llk_cyc210_bp1
#### Looks good! ###
plot(bp1_candidates,bp1_all_llk,type = 'l',
     xlim = c(70,75),ylim = c(2330,2335))

#### optimal likelihood - 1.92 ####
abline(h = max(bp1_all_llk) - 1.92,col = 'red')
abline(v = 71.03,col = 'grey')
abline(v = 72.38,col = 'grey')


#### Profile on endpoint, bp2 #######
bp2_candidates = seq(156,196,10)

profile_llk_cyc210_bp2 = list()
i = 1

for(bp2 in bp2_candidates) {
    
    profile_llk_cyc210_bp2[[i]] = profile_llk_2_stepfunc(data = zeta_perc_avg_df$smoothed_zeta_score,
                                                  cycle = 210,
                                                  param_to_profile = 'bp2',
                                                  profile_param_value = bp2,
                                                  bp1_lower_limit = 70,
                                                  bp1_upper_limit = 80,
                                                  bp2_lower_limit = 170,
                                                  bp2_upper_limit = 180,
                                                  val1_lower_limit = 0.9,
                                                  val1_upper_limit = 1,
                                                  val2_lower_limit = 0.9,
                                                  val2_upper_limit = 1,
                                                  repeats = 100)
    i = i + 1
}
save(profile_llk_cyc210_bp2,file = 'Rdata/Profile on bp2 2-stepfunc_zeta_sigma_0.1_cyc210.Rdata')
profile_llk_cyc210_bp2
bp2_all_llk = sapply(1:length(profile_llk_cyc210_bp2),function(x) {profile_llk_cyc210_bp2[[x]]$llk})

#### Looks good! ###
plot(bp2_candidates,bp2_all_llk,type = 'l',
     ylim = c(2330,2335),xlim = c(170,180))

#### optimal likelihood - 1.92 ####
abline(h = max(bp2_all_llk) - 1.92,col = 'red')
abline(v = 175.2,col = 'grey')
abline(v = 176.8,col = 'grey')


####### Profile on val1 ########
val1_cand = seq(0.83,1.23,0.1)

profile_llk_cyc210_val1 = list()
i = 1

for(val1 in val1_cand) {
    
    profile_llk_cyc210_val1[[i]] = profile_llk_2_stepfunc(data = zeta_perc_avg_df$smoothed_zeta_score,
                                                   cycle = 210,
                                                   param_to_profile = 'val1',
                                                   profile_param_value = val1,
                                                   bp1_lower_limit = 70,
                                                   bp1_upper_limit = 80,
                                                   bp2_lower_limit = 170,
                                                   bp2_upper_limit = 180,
                                                   val1_lower_limit = 0.9,
                                                   val1_upper_limit = 1,
                                                   val2_lower_limit = 0.9,
                                                   val2_upper_limit = 1,
                                                   repeats = 100)
    i = i + 1
}
save(profile_llk_cyc210_val1,file = 'Rdata/Profile on val1 2-stepfunc_zeta_sigma_0.1_cyc210.Rdata')
profile_llk_cyc210_val1

val1_all_llk = sapply(1:length(profile_llk_cyc210_val1),function(x) {profile_llk_cyc210_val1[[x]]$llk})

plot(val1_cand,val1_all_llk,type = 'l',
     xlim = c(1.02,1.04),
     ylim = c(2330,2335))

abline(h = max(val1_all_llk) - 1.92,col = 'red')
abline(v = 1.0294,col = 'grey')
abline(v = 1.0307,col = 'grey')


####### Profile on val2 ########
val2_cand = seq(0.76,1.16,0.1)

profile_llk_cyc210_val2 = list()
i = 1

for(val2 in val2_cand) {
    
    profile_llk_cyc210_val2[[i]] = profile_llk_2_stepfunc(data = zeta_perc_avg_df$smoothed_zeta_score,
                                                   cycle = 210,
                                                   param_to_profile = 'val2',
                                                   profile_param_value = val2,
                                                   bp1_lower_limit = 70,
                                                   bp1_upper_limit = 80,
                                                   bp2_lower_limit = 170,
                                                   bp2_upper_limit = 180,
                                                   val1_lower_limit = 0.9,
                                                   val1_upper_limit = 1,
                                                   val2_lower_limit = 0.9,
                                                   val2_upper_limit = 1,
                                                   repeats = 100)
    i = i + 1
}
save(profile_llk_cyc210_val2,file = 'Rdata/Profile on val2 2-stepfunc_zeta_sigma_0.1_cyc210.Rdata')
profile_llk_cyc210_val2

val2_all_llk = sapply(1:length(profile_llk_cyc210_val2),function(x) {profile_llk_cyc210_val2[[x]]$llk})

plot(val2_cand,val2_all_llk,type = 'l',
     xlim = c(0.95,0.965),ylim = c(2330,2335))

abline(h = max(val2_all_llk) - 1.92,col = 'red')
abline(v = 0.959,col = 'grey')
abline(v = 0.961,col = 'grey')





