######### Profile likelihood on stepfunction on zeta ###########################
rm(list = ls())
load('Rdata/Overall.ILIplus.21-day.aggregated.pcr_holidayNA.RData')
source('functions/profile likelihood function on 2-stepfunction.R')

#### Profile on startpoint, bp1 #######
bp1_candidates = seq(35,75,10)

profile_llk_bp1 = list()
i = 1

for(bp1 in bp1_candidates) {
    
    profile_llk_bp1[[i]] = profile_llk_2_stepfunc(data = ili_plus_21d_df$ili_plus_smth_7d,
                                                  cycle = 365,
                                                  param_to_profile = 'bp1',
                                                  profile_param_value = bp1,
                                                  bp1_lower_limit = 50,
                                                  bp1_upper_limit = 60,
                                                  bp2_lower_limit = 290,
                                                  bp2_upper_limit = 330,
                                                  val1_lower_limit = 0.25,
                                                  val1_upper_limit = 0.3,
                                                  val2_lower_limit = 0.15,
                                                  val2_upper_limit = 0.2,
                                                  sd = 0.166,
                                                  repeats = 100)
    i = i + 1
}
save(profile_llk_bp1,file = 'Rdata/profile_lk_bp1_2-stepfunc_cyc365_iliplus_holidayNA.Rdata')

bp1_all_llk = sapply(1:length(profile_llk_bp1),function(x) {profile_llk_bp1[[x]]$llk})

#### Looks good! ###
plot(bp1_candidates,bp1_all_llk,type = 'l',xlim = c(50,60),ylim = c(1014,1018))

profile_llk_bp1
#### optimal likelihood - 1.92 ####
abline(h = max(bp1_all_llk) - 1.92,col = 'red')
abline(v = 53.3,col = 'grey')
abline(v = 56.1,col = 'grey')


#### Profile on endpoint, bp2 #######
bp2_candidates = seq(304,344,10)

profile_llk_bp2 = list()
i = 1

for(bp2 in bp2_candidates) {
    
    profile_llk_bp2[[i]] = profile_llk_2_stepfunc(data = ili_plus_21d_df$ili_plus_smth_7d,
                                                  cycle = 365,
                                                  param_to_profile = 'bp2',
                                                  profile_param_value = bp2,
                                                  bp1_lower_limit = 50,
                                                  bp1_upper_limit = 60,
                                                  bp2_lower_limit = 300,
                                                  bp2_upper_limit = 330,
                                                  val1_lower_limit = 0.25,
                                                  val1_upper_limit = 0.3,
                                                  val2_lower_limit = 0.15,
                                                  val2_upper_limit = 0.2,
                                                  sd = 0.166,
                                                  repeats = 100)
    i = i + 1
}
save(profile_llk_bp2,file = 'Rdata/profile_lk_bp2_2-stepfunc_cyc365_iliplus_holidayNA.Rdata')
profile_llk_bp2
bp2_all_llk = sapply(1:length(profile_llk_bp2),function(x) {profile_llk_bp2[[x]]$llk})

#### Looks good! ###
plot(bp2_candidates,bp2_all_llk,
     type = 'l',
     xlim = c(300,330),
     ylim = c(1013,1018))

#### optimal likelihood - 1.92 ####
abline(h = max(bp2_all_llk) - 1.92,col = 'red')
abline(v = 304.3,col = 'grey')
abline(v = 320.6,col = 'grey')


####### Profile on val1 ########
val1_cand = seq(0.087,0.487,0.1)

profile_llk_val1 = list()
i = 1

for(val1 in val1_cand) {
    
    profile_llk_val1[[i]] = profile_llk_2_stepfunc(data = ili_plus_21d_df$ili_plus_smth_7d,
                                                   cycle = 365,
                                                   param_to_profile = 'val1',
                                                   profile_param_value = val1,
                                                   bp1_lower_limit = 50,
                                                   bp1_upper_limit = 60,
                                                   bp2_lower_limit = 290,
                                                   bp2_upper_limit = 300,
                                                   val1_lower_limit = 0.25,
                                                   val1_upper_limit = 0.3,
                                                   val2_lower_limit = 0.15,
                                                   val2_upper_limit = 0.2,
                                                   repeats = 100)
    i = i + 1
}
save(profile_llk_val1,file = 'Rdata/profile_lk_val1_2-stepfunc_cyc365_iliplus_holidayNA.Rdata')
profile_llk_val1

val1_all_llk = sapply(1:length(profile_llk_val1),function(x) {profile_llk_val1[[x]]$llk})

plot(val1_cand,val1_all_llk,
     xlim = c(0.286,0.29),
     ylim = c(20,30),
     type = 'l')

abline(h = max(val1_all_llk) - 1.92,col = 'red')
abline(v = 0.2865,col = 'grey')
abline(v = 0.28845,col = 'grey')


####### Profile on val2 ########
val2_cand = seq(0.085,0.285,0.05)

profile_llk_val2 = list()
i = 1

for(val2 in val2_cand) {
    
    profile_llk_val2[[i]] = profile_llk_2_stepfunc(data = ili_plus_21d_df$ili_plus_smth_7d,
                                                   cycle = 365,
                                                   param_to_profile = 'val2',
                                                   profile_param_value = val2,
                                                   bp1_lower_limit = 50,
                                                   bp1_upper_limit = 60,
                                                   bp2_lower_limit = 290,
                                                   bp2_upper_limit = 300,
                                                   val1_lower_limit = 0.25,
                                                   val1_upper_limit = 0.3,
                                                   val2_lower_limit = 0.15,
                                                   val2_upper_limit = 0.2,
                                                   repeats = 100)
    i = i + 1
}
save(profile_llk_val2,file = 'Rdata/profile_lk_val2_2-stepfunc_cyc365_iliplus_holidayNA.Rdata')
profile_llk_val2

val2_all_llk = sapply(1:length(profile_llk_val2),function(x) {profile_llk_val2[[x]]$llk})

plot(val2_cand,val2_all_llk,
    xlim = c(0.183,0.186),
    ylim = c(20,30),
     type = 'l')

abline(h = max(val2_all_llk) - 1.92,col = 'red')
abline(v = 0.1844,col = 'grey')
abline(v = 0.18542,col = 'grey')



############## Cycle = 385 ###############
######### Profile likelihood on stepfunction on zeta ###########################
rm(list = ls())
load('Rdata/Overall.ILIplus.21-day.aggregated.pcr.RData')
source('functions/profile likelihood function on 2-stepfunction.R')

#### Profile on startpoint, bp1 #######
bp1_candidates = seq(28.6,68.6,10)

profile_llk_cyc385_bp1 = list()
i = 1

for(bp1 in bp1_candidates) {
    
    profile_llk_cyc385_bp1[[i]] = profile_llk_2_stepfunc(data = ili_plus_21d_df$ili_plus_smth_7d,
                                                         cycle = 385,
                                                         param_to_profile = 'bp1',
                                                         profile_param_value = bp1,
                                                         bp1_lower_limit = 40,
                                                         bp1_upper_limit = 50,
                                                         bp2_lower_limit = 270,
                                                         bp2_upper_limit = 280,
                                                         val1_lower_limit = 0.3,
                                                         val1_upper_limit = 0.35,
                                                         val2_lower_limit = 0.13,
                                                         val2_upper_limit = 0.5,
                                                         repeats = 100)
    i = i + 1
}
save(profile_llk_cyc385_bp1,file = 'Rdata/profile_lk_bp1_2-stepfunc_cyc385_iliplus_holiday0_sigma0.1.Rdata')
load('Rdata/profile_lk_bp1_2-stepfunc_cyc385_iliplus_holiday0_sigma0.1.Rdata')

bp1_all_llk = sapply(1:length(profile_llk_cyc385_bp1),function(x) {profile_llk_cyc385_bp1[[x]]$llk})

profile_llk_cyc385_bp1
#### Looks good! ###
plot(bp1_candidates,bp1_all_llk,
     type = 'l')

#### optimal likelihood - 1.92 ####
abline(h = max(bp1_all_llk) - 1.92,col = 'red')
abline(v = 48.33,col = 'grey')
abline(v = 49.12,col = 'grey')


#### Profile on endpoint, bp2 #######
bp2_candidates = seq(256,296,10)

profile_llk_cyc385_bp2 = list()
i = 1

for(bp2 in bp2_candidates) {
    
    profile_llk_cyc385_bp2[[i]] = profile_llk_2_stepfunc(data = ili_plus_21d_df$ili_plus_smth_7d,
                                                         cycle = 385,
                                                         param_to_profile = 'bp2',
                                                         profile_param_value = bp2,
                                                         bp1_lower_limit = 40,
                                                         bp1_upper_limit = 50,
                                                         bp2_lower_limit = 270,
                                                         bp2_upper_limit = 280,
                                                         val1_lower_limit = 0.3,
                                                         val1_upper_limit = 0.35,
                                                         val2_lower_limit = 0.13,
                                                         val2_upper_limit = 0.5,
                                                         repeats = 100)
    i = i + 1
}
save(profile_llk_cyc385_bp2,file = 'Rdata/profile_lk_bp2_2-stepfunc_cyc385_iliplus_holiday0_sigma0.1.Rdata')
load('Rdata/profile_lk_bp2_2-stepfunc_cyc365_iliplus_holiday0_sigma0.1.Rdata')

profile_llk_cyc385_bp2
bp2_all_llk = sapply(1:length(profile_llk_cyc385_bp2),function(x) {profile_llk_cyc385_bp2[[x]]$llk})

#### Looks good! ###
plot(bp2_candidates,bp2_all_llk,type = 'l',
     ylim = c(585,590),
     xlim = c(275,280))

#### optimal likelihood - 1.92 ####
abline(h = max(bp2_all_llk) - 1.92,col = 'red')
abline(v = 275.05,col = 'grey')
abline(v = 276.9,col = 'grey')


####### Profile on val1 ########
val1_cand = seq(0.102,0.502,0.1)

profile_llk_cyc385_val1 = list()
i = 1

for(val1 in val1_cand) {
    
    profile_llk_cyc385_val1[[i]] = profile_llk_2_stepfunc(data = ili_plus_21d_df$ili_plus_smth_7d,
                                                          cycle = 385,
                                                          param_to_profile = 'val1',
                                                          profile_param_value = val1,
                                                          bp1_lower_limit = 40,
                                                          bp1_upper_limit = 50,
                                                          bp2_lower_limit = 270,
                                                          bp2_upper_limit = 280,
                                                          val1_lower_limit = 0.3,
                                                          val1_upper_limit = 0.35,
                                                          val2_lower_limit = 0.13,
                                                          val2_upper_limit = 0.5,
                                                          repeats = 100)
    i = i + 1
}
save(profile_llk_cyc385_val1,file = 'Rdata/profile_lk_val1_2-stepfunc_cyc385_iliplus_holiday0_sigma0.1.Rdata')
profile_llk_cyc385_val1

val1_all_llk = sapply(1:length(profile_llk_cyc385_val1),function(x) {profile_llk_cyc385_val1[[x]]$llk})

plot(val1_cand,val1_all_llk,type = 'l',
     xlim = c(0.3,0.305),
     ylim = c(585,590))

abline(h = max(val1_all_llk) - 1.92,col = 'red')
abline(v = 0.3016,col = 'grey')
abline(v = 0.3025,col = 'grey')


####### Profile on val2 ########
val2_cand = seq(0.035,0.235,0.05)

profile_llk_cyc385_val2 = list()
i = 1

for(val2 in val2_cand) {
    
    profile_llk_cyc385_val2[[i]] = profile_llk_2_stepfunc(data = ili_plus_21d_df$ili_plus_smth_7d,
                                                          cycle = 385,
                                                          param_to_profile = 'val2',
                                                          profile_param_value = val2,
                                                          bp1_lower_limit = 40,
                                                          bp1_upper_limit = 50,
                                                          bp2_lower_limit = 270,
                                                          bp2_upper_limit = 280,
                                                          val1_lower_limit = 0.3,
                                                          val1_upper_limit = 0.35,
                                                          val2_lower_limit = 0.13,
                                                          val2_upper_limit = 0.5,
                                                          repeats = 100)
    i = i + 1
}
save(profile_llk_cyc385_val2,file = 'Rdata/profile_lk_val2_2-stepfunc_cyc385_iliplus_holiday0_sigma0.1.Rdata')
profile_llk_cyc385_val2

val2_all_llk = sapply(1:length(profile_llk_cyc385_val2),function(x) {profile_llk_cyc385_val2[[x]]$llk})

plot(val2_cand,val2_all_llk,type = 'l',
     xlim = c(0.134,0.136),ylim = c(585,590))

abline(h = max(val2_all_llk) - 1.92,col = 'red')
abline(v = 0.13447,col = 'grey')
abline(v = 0.1355,col = 'grey')





