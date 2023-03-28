rm(list = ls())
library(ggplot2)
library(RColorBrewer)
source('functions/Multiple_stepfunc_v5.R')
load('Rdata/7d.smth.zeta.ilipercent.Rdata')

####### Read all the fits #######
all_fits_pt1 = list()
for(i in 1:100) {
    
    load(paste0('Rdata/Results from Roar/Roar_zeta_minAvailDay_150_sensitivity_analysis_1_100/best_fit_',i,'.RData'))
    all_fits_pt1[[i]] = best_fit
}

all_fits_pt2 = list()
for(i in 1:95) {
    
    load(paste0('Rdata/Results from Roar/Roar_zeta_minAvailDay_150_sensitivity_analysis_101_195/best_fit_',i,'.RData'))
    all_fits_pt2[[i]] = best_fit
}

all_fits_pt3 = list()
for(i in 1:40) {
    
    load(paste0('Rdata/Results from Roar/Roar_zeta_minAvailDay_150_sensitivity_analysis_cyc150_189/best_fit_',i,'.RData'))
    all_fits_pt3[[i]] = best_fit
}

all_fits_pt4 = list()
for(i in 1:70) {
    
    load(paste0('Rdata/Results from Roar/Roar_zeta_minAvailDay_150_sensitivity_analysis_cyc385_450/best_fit_',i,'.RData'))
    all_fits_pt4[[i]] = best_fit
}

all_fits = c(all_fits_pt1,all_fits_pt2,all_fits_pt3,all_fits_pt4)

all_N = sapply(1:length(all_fits),function(x){all_fits[[x]]$N})
all_aic = sapply(1:length(all_fits),function(x){all_fits[[x]]$aic})
all_cycle = sapply(1:length(all_fits), function(x){all_fits[[x]]$cycle})
cycles = unique(all_cycle)

all_results = data.frame(cycle = all_cycle, N = all_N,aic = all_aic)
all_results[which.min(all_results$aic),]
all_results = all_results[order(all_results$aic),]


selected_results = all_results[all_results$aic < 0,]
ggplot(selected_results) +
    geom_raster(aes(x = cycle,y = N,fill = aic)) +
    scale_fill_viridis_c(option = "H") +
    scale_x_continuous(breaks = seq(min(cycles),max(cycles),15)) +
    labs(fill = 'AIC',tag = 'A') +
    xlab('Cycle') +
    ylab('N') +
    theme_bw() +
    theme(text = element_text(size = 25),
          axis.text = element_text(size = 20),
          plot.margin = margin(t = 55,b = 10,r = 5,l = 0),
          plot.tag = element_text(size = 35),
          plot.tag.position = c(0.05,1.1))
# ggsave(filename = '../plots/zeta step function/annotated_grid_search_stepfunc_zeta_sigma_unfixed.jpg',
#            width = 5000, height = 2000, units = c('px'))


### Get the best step function for each cycle ###

# convert df to matrix 
m_results = reshape(all_results,v.names = 'aic',idvar = 'cycle',
                    timevar = 'N',direction = 'wide')

best_aic_per_cycle = apply(m_results[,-1],1,min)

best_aic_per_cycle = data.frame(cycle = m_results$cycle,
                                aic = best_aic_per_cycle)

best_aic_per_cycle = best_aic_per_cycle[order(best_aic_per_cycle$aic),]

plot(m_results$cycle,best_aic_per_cycle,pch = 20,xaxt = 'n')
axis(1,at = unique(all_cycle),labels = unique(all_cycle))

#### Does any N is best fitting the data ? ###
all_best_N_cycle = data.frame()

for(cycle in unique(all_cycle)) {
    
    results_per_cycle = all_results[all_results$cycle == cycle,]
    N_lowest_aic = results_per_cycle$N[which.min(results_per_cycle$aic)]
    
    all_best_N_cycle = rbind(all_best_N_cycle,
                             data.frame(cycle = cycle,
                                        best_N = N_lowest_aic))
}
hist(all_best_N_cycle$best_N)
plot(all_best_N_cycle$cycle,all_best_N_cycle$best_N,pch = 20,cex = 0.5)

### Extract the best step function from each cycle ###
all_sim_ts = c()

for(cycle in unique(all_cycle)) {
    
    indices_aic = data.frame(index = which(all_results$cycle == cycle),
                             aic = all_results$aic[all_results$cycle == cycle])
    index = indices_aic$index[which.min(indices_aic$aic)]
    
    sim_ts = sapply(1:nrow(zeta_perc_avg_df),multi_stepfunc,cycle = cycle,
                    all_fits[[index]]$all_results)
    
    all_sim_ts = rbind(all_sim_ts,sim_ts)
}

dim(all_sim_ts)
rownames(all_sim_ts) = unique(all_cycle)

all_zeta_sim_ts = all_sim_ts
save(all_zeta_sim_ts,file = 'Rdata/zeta_best_stepfunc_for_sensitivity_analysis_cyc150_450.Rdata')
