plot_data_by_year = function(data_df,data_name,time_name,
                             ymin,ymax,smoothed_data_name,
                             Nrow,Ncol,tag,
                             save = F,
                             filename = filename,
                             width = width,
                             height = height) {
    
    data_df$year = format(data_df[,time_name],"%Y")
    years = unique(data_df$year)
    
    if(save) { jpeg(filename = filename,width = width, height = height,
                         units = 'in',res = 300) }
    
    par(mfcol = c(Nrow,Ncol),mar = c(0.5,3.5,0.5,0.8),mgp = c(2,0.5,0),oma = c(1,0.5,1,0.5))
    for(year in years) {
        
        data_each_year = data_df[data_df$year == year,]
        
        plot(data_each_year[,time_name],data_each_year[,data_name],
             type = "l",lwd = 1.5,ylim = c(ymin,ymax),
             xaxt = "n",xaxs = "i",main = '',xlab = '',
             ylab = year,
             yaxt = 'n',
             cex.lab = 2,
             col = 'grey')
        abline(h = mean(data_each_year[,data_name],na.rm = T),col = 'blue',lty = 'dashed')
        if(year %in% c('2014','2019')) {
        axis(1,at = seq(as.Date(paste0(year,"-01-01")),as.Date(paste0(year,"-12-31")),"2 months"),
             labels = format(seq(as.Date(paste0(year,"-01-01")),as.Date(paste0(year,"-12-31")),"2 months"),"%b"),
             cex.axis = 1.7)
        }
        axis(2,at = seq(ymin,ymax,0.2),labels = seq(ymin,ymax,0.2),cex.axis = 1.7)
        abline(v = seq(as.Date(paste0(year,"-01-01")),as.Date(paste0(year,"-12-31")),"1 month"),lty = "dashed",col = "grey")
        
        if(!is.null(smoothed_data_name)) {
            
            lines(x = data_each_year[,time_name],
                  y = data_each_year[,smoothed_data_name],col = "black",lwd = 2)
        }
        
    }
     mtext(tag,side = 3,outer = T,adj = 0,cex = 2)
     
     if(save) {dev.off()}
    
    
    
}


plot_data_with_stepfunc_by_year_fixed_cycle = function(data_df,data_name,time_name,
                                                       stepfunc_params,cycle,
                                                       ymin,ymax,
                                                       Nrow,Ncol) {
    
    data_df$step_ts = sapply(1:nrow(data_df),multi_stepfunc,
                             cycle = cycle, 
                             stepfunc_params)
    data_df$year = format(data_df[,time_name],"%Y")
    years = unique(data_df$year)
    
    par(mfcol = c(Nrow,Ncol),mar = c(2,2,2,0.8),mgp = c(2,0.5,0))
    for(year in years) {
        
        data_each_year = data_df[data_df$year == year,]
        
        plot(data_each_year[,time_name],data_each_year[,data_name],
             type = "l",lwd = 1.5,
             xlim = c(as.Date(paste0(year,"-01-01")),as.Date(paste0(year,"-12-31"))),
             ylim = c(ymin,ymax),
             main = '',
             xaxs = "i",
             xaxt = "n",yaxt = 'n')
        title(main = year,xlab = '',ylab = '',cex.main = 2,line = 0.5)
        axis(1,at = seq(as.Date(paste0(year,"-01-01")),as.Date(paste0(year,"-12-31")),"1 month"),
             labels = format(seq(as.Date(paste0(year,"-01-01")),as.Date(paste0(year,"-12-31")),"1 month"),"%m"),
             cex.axis = 1.4)
        axis(2, at = seq(ymin,ymax,0.2),
             labels = seq(ymin,ymax,0.2),
             cex.axis = 1.4)
        abline(v = seq(as.Date(paste0(year,"-01-01")),as.Date(paste0(year,"-12-31")),"1 month"),
               lty = "dashed",
               col = "grey")
        lines(data_each_year[,time_name],data_each_year$step_ts,col = "red",lwd = 1.5)
        
    }
    
}


plot_data_with_fitted_value_by_year_fixed_cycle = function(data_df,data_name,fitted_value,
                                                           ymin,ymax,
                                                           time_name,Nrow,Ncol) {
    
    data_df$fitted_value = fitted_value
    data_df$year = format(data_df[,time_name],"%Y")
    years = unique(data_df$year)
    
    par(mfcol = c(Nrow,Ncol),mar = c(2,2,2,0.8),mgp = c(2,0.5,0))
    for(year in years) {
        
        data_each_year = data_df[data_df$year == year,]
        
        plot(data_each_year[,time_name],data_each_year[,data_name],
             type = "l",lwd = 1.5,
             xlim = c(as.Date(paste0(year,"-01-01")),as.Date(paste0(year,"-12-31"))),
             ylim = c(ymin,ymax),
             main = '',
             xaxs = "i",
             xaxt = "n",yaxt = 'n')
        title(main = year,xlab = '',ylab = '',cex.main = 2,line = 0.5)
        axis(1,at = seq(as.Date(paste0(year,"-01-01")),as.Date(paste0(year,"-12-31")),"1 month"),
             labels = format(seq(as.Date(paste0(year,"-01-01")),as.Date(paste0(year,"-12-31")),"1 month"),"%m"),
             cex.axis = 1.4)
        axis(2, at = seq(ymin,ymax,0.2),
             labels = seq(ymin,ymax,0.2),
             cex.axis = 1.4)
        abline(v = seq(as.Date(paste0(year,"-01-01")),as.Date(paste0(year,"-12-31")),"1 month"),
               lty = "dashed",
               col = "grey")
        lines(data_each_year[,time_name],data_each_year$fitted_value,col = "forestgreen",lwd = 1.5)
        
    }
    
    
    
}

### ggplot version ####
library(ggplot2)
library(reshape2)

ggplot_data_fitted_value_by_year = function(data_df,data_name,time_name,fitted_value) {
    plot_data_df = data.frame(date = data_df[,time_name],
                              data = data_df[,data_name],
                              fitted_value = fitted_value)
    molten_plot_data_df = melt(plot_data_df,id.vars = 'date')
    
    molten_plot_data_df$year = format(molten_plot_data_df$date,'%Y')
    molten_plot_data_df$date = as.Date(paste0('2020-',format(molten_plot_data_df$date,'%m-%d')))
    molten_plot_data_df$fake_year = factor(molten_plot_data_df$year,
                                           levels = seq(2019,2010,-1),
                                           labels = seq(2019,2010,-1))
    
    ggplot(molten_plot_data_df) +
        geom_line(aes(x = date,y = value,color = variable)) + 
        scale_color_manual(values = c(data = 'black',fitted_value = 'red'),
                           labels = c(data = data_name,
                                      fitted_value = 'fitted value')) +
        scale_x_date(date_breaks = '1 month',date_labels = '%m',expand = c(0,0)) +
        xlab('') +
        ylab('') +
        facet_grid(fake_year~.) +
        theme(plot.title = element_text(hjust = 0.5),
              aspect.ratio = 2/10,
              panel.background = element_blank(),
              panel.border = element_rect(color = "grey",inherit.blank = T,fill = NA),
              panel.grid.major.x = element_line(linetype = "dashed",color = "grey"),
              legend.position = "bottom",
              panel.spacing = unit(0.5,"lines")) +
        guides(color = guide_legend(title = "",nrow = 2))
    
}


