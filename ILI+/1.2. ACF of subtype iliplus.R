#### ACF of subtype ILI+ ####
rm(list = ls())
setwd("~/Dropbox/Influenza and Respiratory Disease in Vietnam/code")
load('Rdata/Subtype.iliplus.21d.aggregated.pcr.Rdata')


jpeg(filename = 'plots/SuppFig8.jpg',width = 10,height = 8,
     units = 'in',res = 300)
par(mfcol = c(3,2),mar = c(2,3,2,2))
ilip_names = c('H1N1 ILI+','H3N2 ILI+','B ILI+')
for(i in 1:length(all_subtype_ilip)) {
    
    ilip_spec = spectrum(all_subtype_ilip[[i]]$ilip_smth_7d,
                         plot = F,log = 'no',
                         pad = 0.051,
                         fast = F)
    
    cycle_spec = data.frame(cycle = 1/ilip_spec$freq,
                            spec = ilip_spec$spec)
    cycle_spec = cycle_spec[order(cycle_spec$spec,decreasing = T),]
    
    plot(cycle_spec$cycle,cycle_spec$spec,type = 'h',
         xlim = c(0,600),
         xlab = "",ylab = "",main = "",
         lwd = 2,mgp = c(2,1,0),cex.axis = 1.5)
    title(main = ilip_names[i],line = 0.5,
          cex.main = 2)
    # points(x = cycle_spec$cycle[1],y = cycle_spec$spec[1],
    #        pch = 20,cex = 1.5)
}


all_subtype_acf = list()
for(i in 1:length(all_subtype_ilip)) {
        
        all_subtype_acf[[i]] = acf(all_subtype_ilip[[i]]$ilip_smth_7d,
                                   lag.max = 365*3,
                                   ylim = c(-0.2,0.5),
                                   col = 'grey',
                                   main = '',
                                   xaxt = 'n',
                                   cex.axis = 1.4)
        title(main = ilip_names[i],line = 0.5,
              cex.main = 2)
        max_lag = all_subtype_acf[[i]]$lag[which.max(all_subtype_acf[[i]]$acf[150:450]) + 149]
        max_cycle = seq(max_lag,max_lag*4,max_lag)
        axis(1, at = max_cycle, labels = max_cycle)
        abline(v = max_cycle, lty = 'dashed',col = 'black')
        points(x = seq(365,1095,365),
               y = c(all_subtype_acf[[i]]$acf[365],
                     all_subtype_acf[[i]]$acf[730],
                     all_subtype_acf[[i]]$acf[1095]),
               pch = 20,
               cex = 1.5)
        
}
dev.off()

