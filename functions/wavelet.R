library(forecast)
library(Rwave)
library(WaveletComp)
library(ggplot2)
library(metR)
library(reshape2)


### function of getting Morlet wavelet ###
wavelet_transform = function(no,nv,start_cycle,end_cycle,
                             time_name,data_name,df) {
    
    #### Find the cycle resolution ###

    all_periods = 2^seq(1,no +1-1/nv,1/nv)
    max_period = 2^(no+1-1/nv)
    
    if(start_cycle < 2 | end_cycle > max_period) {
        
        warning('End cycle is larger than max period in wavelet transform, 
                higher no is suggested.')
    }
    wfit = cwt(df[,data_name],
               no,nv,plot = F)
    
    wspec = Mod(wfit)
    rownames(wspec) = df[,time_name]
    colnames(wspec) = all_periods
    
    return(list(wspec = wspec,all_periods = all_periods))
}


### function of getting cone of influence 
get_coi = function(df,time_name,all_periods) {
    
    ### cone of influence ###
    coi = COI(start = 1,
              dt = 1,
              nc = nrow(df),
              Period = all_periods)
    
    ### Get the coordinates of coi ###
    coi_coords_x = coi$x
    coi_coords_y = coi$y
    
    ### remove the padding zeros at the edges ###
    num_pad_0s_x = (length(coi_coords_x) - nrow(df))/2
    coi_coords_x = coi_coords_x[(1 + num_pad_0s_x):(length(coi_coords_x) - num_pad_0s_x)]
    
    num_pad_0s_y = (length(coi_coords_y) - nrow(df))/2
    coi_coords_y = coi_coords_y[(1 + num_pad_0s_y):(length(coi_coords_y) - num_pad_0s_y)]
    
    coi_coords = data.frame(x = coi_coords_x,
                            y = coi_coords_y)
    
    ### change coi_coords_x to date to accommodate with date x-axis ###
    ### y is the exponent of 2, calculate the actual period from y ###
    coi_coords$x = df[,time_name]
    coi_coords$y = 2^coi_coords$y
    
    return(coi_coords)
}

#### plot wavelet and coi ###
plot_wavelet = function(wspec,
                        coi_coords,
                        time_name,
                        tile_width,
                        tile_height,
                        contour_binwidth,
                        contour_threshold,
                        date_breaks = date_breaks,
                        cycle_breaks,
                        annual_cycle,
                        start_cycle,
                        end_cycle,
                        tag,
                        signal_range,
                        filename,
                        change_tag_margin,
                        margins,
                        tag.coordiantes,
                        pic_width,
                        pic_height) {
    
    molten_wspec = melt(wspec)
    colnames(molten_wspec) = c('date','cycle','coefficients')
    selected_molten_wspec = molten_wspec[molten_wspec$cycle >= start_cycle & molten_wspec$cycle <= end_cycle,]
    selected_molten_wspec$date = as.Date(selected_molten_wspec$date,origin = '1970-01-01')
    
    p = ggplot() +
        geom_tile(aes(x = date,y = cycle, fill = coefficients),
                  width = tile_width,height = tile_height,
                  data = selected_molten_wspec) +
        geom_line(aes(x = x,y = y),
                  col = 'black',
                  data = coi_coords) +
        geom_contour(aes(x = date,y = cycle, z = coefficients),
                     color = 'grey',
                     show.legend = T,
                     binwidth = contour_binwidth,
                     data = selected_molten_wspec[selected_molten_wspec$coefficients > contour_threshold,]) +
        # geom_text_contour(aes(x = date, y = cycle, z = coefficients),
        #                   size = 0.5,color = 'grey',
        #                   data = selected_molten_wspec[selected_molten_wspec$coefficients > contour_threshold,]) +
        geom_hline(yintercept = annual_cycle,linetype = 'dashed') +
        coord_cartesian(xlim = range(selected_molten_wspec$date),ylim = range(selected_molten_wspec$cycle)) +
        scale_fill_viridis_c(option = "H",limit = signal_range) +
        scale_x_date(date_breaks = date_breaks,date_labels = '%Y',expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0),n.breaks = cycle_breaks) + 
        labs(tag = tag) +
        theme(text = element_text(size = 20),
              legend.title = element_blank()) 
    
    if(change_tag_margin) {
        
        p = p + 
            theme(plot.margin = margin(margins),
                  plot.tag.position = tag.coordiantes)
        
    }
    ggsave(filename = paste0('../plots/wavelet/',filename,'.jpg'),
           plot = p,
           units = 'px',
           width = pic_width,height = pic_height)
}

##### wrapping all functions ###
wavelet_wrapper = function(df,time_name,data_name,
                           start_cycle,
                           end_cycle,
                           no,
                           nv,
                           tile_width,
                           tile_height,
                           contour_binwidth,
                           contour_threshold,
                           date_breaks = '1 year',
                           cycle_breaks = 10,
                           pic_width = 2000,
                           pic_height = 1200,
                           annual_cycle,
                           signal_range,
                           tag,
                           filename,
                           change_tag_margin = F,
                           margins = NULL,
                           tag.coordiantes = NULL) {
    
    ### get Morlet wavelet ###
    results = wavelet_transform(no = no,
                                nv = nv,
                                start_cycle = start_cycle,
                                end_cycle = end_cycle,
                                time_name = time_name,
                                data_name = data_name,
                                df = df)
    
    ### get COI of zeta ##
    coi_coords = get_coi(df = df,
                       time_name = time_name,
                       all_periods = results$all_periods)
    
    
    ### plot wavelet ###
    plot_wavelet(wspec = results$wspec,
                 coi_coords = coi_coords,
                 time_name = time_name,
                 tile_width = tile_width,
                 tile_height = tile_height,
                 contour_binwidth = contour_binwidth,
                 contour_threshold = contour_threshold,
                 annual_cycle = annual_cycle,
                 date_breaks = date_breaks,
                 cycle_breaks = cycle_breaks,
                 start_cycle = start_cycle,
                 end_cycle = end_cycle,
                 tag = tag,
                 signal_range = signal_range,
                 filename = filename,
                 change_tag_margin = change_tag_margin,
                 margins = margins,
                 tag.coordiantes = tag.coordiantes,
                 pic_width = pic_width,
                 pic_height = pic_height)
    
}
