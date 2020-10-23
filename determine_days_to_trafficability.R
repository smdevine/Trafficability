laptop <- FALSE
if (laptop) {
  resultsDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/HYDRUS_runs'} else {
    resultsDir <- 'D:/PostDoc/Trafficability/Oct2020test'
    summaryDir <- file.path(resultsDir, 'summary')
}




#function to summarize trafficability based on 100 day model where start time is 24 hours before
#flood_start for a number of days described by flood_duration
#trafficability assumed to be equal to 90% of field capacity
# flood_duration <- 4
# f_path <- 'Delhi/Delhi'
# FC_factor <- 0.9
specific_depth <- function(depth, vec_length, FC_factor, flood_end, df) {
  df$time_days[head(which(df[[paste0('theta_', depth, 'cm')]][2:vec_length] < df[[paste0('theta_', depth, 'cm')]][1] * FC_factor), n=1)+1] - flood_end
}

determine_trafficability <- function(f_path, resultsDir, flood_duration) {
  print(f_path)
  # flood_end <- as.integer(format(as.Date(flood_start), format='%j')) + flood_duration - 1
  flood_end <- flood_duration + 1
  if(!file.exists(file.path(resultsDir, f_path, f_path, 'Obs_Node.out'))) {
    print(paste('No results for', f_path))
    next
  }
  file.remove(file.path(resultsDir, f_path, f_path, 'Nod_Inf.out'), showWarnings=FALSE)
  file.remove(file.path(resultsDir, f_path, f_path, 'Nod_Inf_V.out'), showWarnings=FALSE)
  obs <- read.table(file.path(resultsDir, f_path, f_path, 'Obs_Node.out'), col.names = c('time_days', 'h_0cm', 'theta_0cm', 'flux_0cm', 'h_1cm', 'theta_1cm', 'flux_1cm', 'h_2cm', 'theta_2cm', 'flux_2cm', 'h_3cm', 'theta_3cm', 'flux_3cm', 'h_4cm', 'theta_4cm', 'flux_4cm', 'h_5cm', 'theta_5cm', 'flux_5cm', 'h_6cm', 'theta_6cm', 'flux_6cm', 'h_7cm', 'theta_7cm', 'flux_7cm', 'h_8cm', 'theta_8cm', 'flux_8cm', 'h_9cm', 'theta_9cm', 'flux_9cm', 'h_10cm', 'theta_10cm', 'flux_10cm', 'h_15cm', 'theta_15cm', 'flux_15cm', 'h_20cm', 'theta_20cm', 'flux_20cm', 'h_30cm', 'theta_30cm', 'flux_30cm', 'h_50cm', 'theta_50cm', 'flux_50cm'), header = FALSE, nrow=948, skip=11)
  obs_10cm_theta <- obs[,c(1,3,6,9,12,15,18,21,24,27,30,33)]
  vec_length <- length(obs_10cm_theta$theta_5cm)
  obs_0_10cm_theta_avg <- apply(obs_10cm_theta[2:nrow(obs_10cm_theta), 2:12], 1, mean)
  
  depth_results_0.95 <- sapply(c(2,5,10), specific_depth, vec_length = vec_length, FC_factor= 0.95, flood_end = flood_end, df = obs_10cm_theta)
  depth_results_0.9 <- sapply(c(2,5,10), specific_depth, vec_length = vec_length, FC_factor= 0.9, flood_end = flood_end, df =obs_10cm_theta)
  depth_results_0.85 <- sapply(c(2,5,10), specific_depth, vec_length = vec_length, FC_factor= 0.85, flood_end = flood_end, df =obs_10cm_theta)
  depth_results_0.8 <- sapply(c(2,5,10), specific_depth, vec_length = vec_length, FC_factor= 0.8, flood_end = flood_end, df =obs_10cm_theta)
  depth_results <- sapply(c(depth_results_0.95, depth_results_0.9, depth_results_0.85, depth_results_0.8), function(x) if(length(x)==0){NA} else{x})
  
  #0-10 cm avg
  days_to_traff_0_10cm <- sapply(c(0.95, 0.9, 0.85, 0.8), function(x) {
    obs_10cm_theta$time_days[head(which(obs_0_10cm_theta_avg < apply(obs_10cm_theta[1, 3:12], 1, mean) * x), n=1)+1] - flood_end}) #FC estimate ignores upper layer
  
  results <- matrix(data=depth_results, nrow = 3, ncol = 4, byrow = FALSE, dimnames = list(c('days_to_traff_2cm', 'days_to_traff_5cm', 'days_to_traff_10cm'), c('FC_0.95', 'FC_0.9', 'FC_0.85', 'FC_0.8')))
  results <- rbind(results, days_to_traff_0_10cm)
  write.csv(results, file = file.path(summaryDir, paste0(f_path, '_results.csv')), row.names=TRUE)
}
# determine_trafficability(f_path = 'Delhi/Delhi', resultsDir = resultsDir, FC_factor=0.9, flood_duration=4)
# determine_trafficability(f_path = 'Hanford/Hanford', resultsDir = resultsDir, FC_factor=0.9, flood_duration=4)
# determine_trafficability(f_path = 'Yolo/Yolo', resultsDir = resultsDir, FC_factor=0.9, flood_duration=4)
# determine_trafficability(f_path = 'Panoche_4/Panoche_4', resultsDir = resultsDir, FC_factor=0.9, flood_duration=4)
# determine_trafficability(f_path = 'Yolo_2/Yolo_2', resultsDir = resultsDir, FC_factor=0.9, flood_duration=4)
# determine_trafficability(f_path = 'soil_166', resultsDir = resultsDir, flood_duration=4)

path_names <- list.dirs(resultsDir, full.names = FALSE, recursive = FALSE)
head(path_names)
tail(path_names)
path_names <- path_names[2:(length(path_names)-2)]

error_runs <- sapply(path_names, function(x) file.exists(file.path(resultsDir, x, x, 'Error.msg')))

sum(error_runs) #15 of 5509
error_names <- path_names[error_runs]
path_names <- path_names[!error_runs]
length(path_names)
sapply(path_names, determine_trafficability, resultsDir=resultsDir, flood_duration=4)
