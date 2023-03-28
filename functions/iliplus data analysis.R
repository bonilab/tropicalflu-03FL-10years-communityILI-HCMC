##### Plot clinical iliplus and sentinel iliplus #####

plot_sentinel_iliplus_with_clinical_iliplus = function(sentinel_data) {
    
    sentinel_data = sentinel_data[sentinel_data$firstday %in% common_period,]
    
    par(mfcol = c(2,2),mar = c(2,2,2,2))
    for(i in 2:ncol(clinical_weekly_iliplus)) {
        
        plot(clinical_weekly_iliplus$firstday,clinical_weekly_iliplus[,i],
             type = 'l',
             main = colnames(clinical_weekly_iliplus)[i],
             ylim = c(min(cbind(clinical_weekly_iliplus[,i],sentinel_data[,i]),na.rm = T),
                      max(cbind(clinical_weekly_iliplus[,i],sentinel_data[,i]),na.rm = T)))
        lines(sentinel_data$firstday,sentinel_data[,i],
              col = "red")
    }
}


##### Interpolate daily value within 1 week ####
weekly_to_daily_iliplus = function(iliplus_data, daily_period_in_data) {
    
    # remove overall iliplus 
    plot_iliplus_data = iliplus_data[,-2]
    
    ## To ensure the stack plot is continuous, interpolate daily value between weeks ##
    
    plot_daily_iliplus_data = data.frame()
    
    for(i in 1:nrow(plot_iliplus_data)) {
        
        week_start = plot_iliplus_data$firstday[i]
        
        if(i == nrow(plot_iliplus_data)) {
            
            week_end = plot_iliplus_data$firstday[i] + 6
        } else {
            week_end = plot_iliplus_data$firstday[i+1] - 1
        }
        
        week_df = data.frame(date = daily_period_in_data[daily_period_in_data %in% seq(week_start,week_end,'1 day')],
                             H1plus = plot_iliplus_data$H1plus[i],
                             H3plus = plot_iliplus_data$H3plus[i],
                             Bplus = plot_iliplus_data$Bplus[i])
        
        plot_daily_iliplus_data = rbind(plot_daily_iliplus_data,week_df)
        
    }
    
    return(plot_daily_iliplus_data)
    
}


#### Plot stacked iliplus by year #####

library(ggplot2)
library(reshape2)
library(grid)

plot_stacked_iliplus_by_year = function(iliplus_data,
                                        Bplus_label,
                                        H1plus_label,
                                        H3plus_label,
                                        tag = '',
                                        change_tag_position = F,
                                        margins = NULL,
                                        tag_coordiantes = NULL) {
    
    molten_iliplus_data = melt(iliplus_data,id.vars = 'date')
    molten_iliplus_data$year = format(molten_iliplus_data$date,'%Y')
    years = unique(molten_iliplus_data$year)
    molten_iliplus_data$fake_year = factor(molten_iliplus_data$year,levels = seq(max(years),min(years),-1))
    molten_iliplus_data$date = as.Date(paste0('2020',format(molten_iliplus_data$date,'-%m-%d')))
    
    p = ggplot(molten_iliplus_data) +
        geom_col(aes(x = date,y = value,fill = variable),
                 position = "stack") +
        scale_fill_manual(values = c(Bplus = "skyblue",H1plus = "orange",H3plus = "red3"),
                          labels = c(Bplus = Bplus_label,
                                     H1plus = H1plus_label,
                                     H3plus = H3plus_label)) + 
        scale_x_date(date_breaks = "1 month",date_labels = "%m",expand = c(0,0))+
        xlab("Date") +
        ylab("ILI+") +
        labs(tag = tag) +  ## Add annotation ###
        facet_grid(fake_year~.) +
        guides(fill = guide_legend(title = "",nrow = 3)) +
        theme(plot.title = element_text(hjust = 0.5),
              aspect.ratio = 2/10,
              panel.background = element_blank(),
              panel.border = element_rect(color = "grey",inherit.blank = T,fill = NA),
              panel.grid.major.x = element_line(linetype = "dashed",color = "grey"),
              legend.position = "bottom",
              panel.spacing = unit(0.8,"lines"),
              text = element_text(size = 15))
    
    if(change_tag_position) {
        
        p = p + theme(plot.margin = margin(margins),
                      plot.tag.position = tag_coordiantes)
    }
    p
}


##### Calculate correlation ######

correlation_all_iliplus = function(iliplus_data_1, iliplus_data_2) {
    
    cor_test_list = list()
    j = 1
    for(i in 2:ncol(iliplus_data_1)) {
    
        
        cor_test_list[[j]] = cor.test(iliplus_data_1[,i],iliplus_data_2[,i])
        j = j + 1
        
    }
    
    names(cor_test_list) = colnames(iliplus_data_1)[2:ncol(iliplus_data_1)]
    
    return(cor_test_list)
}
