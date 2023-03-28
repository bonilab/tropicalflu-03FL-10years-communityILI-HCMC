# The following packages are for parallelization
library(future)
library(listenv)
library(stringr)
library(future.apply)
library(future.batchtools)

# Set up the environment
load("data.Rdata")
source("func.R")

# Remove stale logs
stale_log_files <- list.files('log_files', '*.log', full.names=T)

if (length(stale_log_files) > 0) {
    print('Remove stale log files in log_files/')
    file.remove(stale_log_files)
}

to_log <- function(msg, row_index) {
    con <- file(paste0('log_files/', str_pad(row_index, 6, pad=0), '.log'), 'a')
    writeLines(msg, con)
    close(con)
}

# Set up the parameter lists for parallelization
all_cycle = seq(385,450,5)
all_N = seq(4,8,1)
iters = 120

#all_cycle = 210
#all_N = c(4,5)
#iters = 24

# These are the arguments to be distributed to cores (compute nodes)
node_args = expand.grid(cycle = all_cycle, N = all_N)

# Specify resources
plan(list(
    tweak(batchtools_torque, 
          resources = list(account = "open", 
                           walltime = 'walltime=08:00:00', 
                           cores = 'select=1:ncpus=24:basic',
                           memory = "pmem=4gb")),
    # multisession,
    multisession
))

# Initialize an empty list to store best results from iterations
# all_best_fits = listenv()

print(paste0('There are ', nrow(node_args), ' arg combinations.'))
print('Start the parallel computation ...')

# Manually set seed
seed = 1

ret = listenv()

# Iterate through all node arguments
for (row_index in c(63,65,69)) {
    
    print('Requesting a new node ...')
    
    ret[[row_index]] %<-% {
        
        to_log(paste('Process started with row_index', row_index, '...'), row_index)
        
        all_fits = listenv()
        
        for (i in 1:iters) {
            all_fits[[i]] %<-% fit_multi_stepfunc(data = zeta_perc_avg_df$smoothed_zeta_score,
                                                  N = node_args[row_index, 'N'],
                                                  cycle = node_args[row_index, 'cycle'],
                                                  bps_lower_limit = 1,
                                                  bps_upper_limit = node_args[row_index, 'cycle'],
                                                  vals_lower_limit = 0.92,
                                                  vals_upper_limit = 1.02,
                                                  sigma_lower_limit = 0.11,
                                                  sigma_upper_limit = 0.13,
                                                  seed = seed)
            
            # print(paste(i,"iters is finished! And seed number is",seed))    
            seed = seed + 1
        }
        
        to_log('Loop finished!', row_index)
        to_log('Converting environment to list ...', row_index)
        
        all_fits = tryCatch({as.list(all_fits)},
                            error = function(err) {
                                to_log(c('Errored during conversion',
                                         as.character(err),
                                         '**** Terminating ****'), row_index)
                                stop()
                            })
        
        to_log('Select the best fit ...', row_index)
        
        best_fit = select_best_fit(all_fits=all_fits,
                                   N=node_args[row_index, 'N'],
                                   cycle=node_args[row_index, 'cycle'],
                                   iters=iters)
        
        ## Save result to an RData file
        out_file = paste0('best_fit_', row_index, '.RData')
        
        to_log(paste('Saving best fit to', out_file, '...'), row_index)
        save(best_fit, file = out_file)
        to_log('Save complete!', row_index)
        
        TRUE
    }
}

print('Computation has been distributed!')


