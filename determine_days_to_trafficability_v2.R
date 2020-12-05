laptop <- FALSE
#function to summarize trafficability based on 100 day model where start time is 24 hours before
#flood_start for a number of days described by flood_duration
determine_trafficability_v2 <- function(f_path, resultsDir, flood_duration, removeFiles, subDir, print_f_path, summaryDir) {
  if(print_f_path) {
    print(f_path)
  }
  # flood_end <- as.integer(format(as.Date(flood_start), format='%j')) + flood_duration - 1
  flood_end <- flood_duration + 1
  if(!file.exists(file.path(resultsDir, subDir, f_path, f_path, 'Obs_Node.out'))) {
    print(paste('No results for', f_path))
    next
  }
  if (removeFiles) {
    file.remove(file.path(resultsDir, subDir, f_path, f_path, 'Nod_Inf.out'), showWarnings=FALSE)
    file.remove(file.path(resultsDir, subDir, f_path, f_path, 'Nod_Inf_V.out'), showWarnings=FALSE)
  }
  obs <- read.table(file.path(resultsDir, subDir, f_path, f_path, 'Obs_Node.out'), col.names = c('time_days', 'h_0cm', 'theta_0cm', 'flux_0cm', 'h_1cm', 'theta_1cm', 'flux_1cm', 'h_2cm', 'theta_2cm', 'flux_2cm', 'h_3cm', 'theta_3cm', 'flux_3cm', 'h_4cm', 'theta_4cm', 'flux_4cm', 'h_5cm', 'theta_5cm', 'flux_5cm', 'h_6cm', 'theta_6cm', 'flux_6cm', 'h_7cm', 'theta_7cm', 'flux_7cm', 'h_8cm', 'theta_8cm', 'flux_8cm', 'h_9cm', 'theta_9cm', 'flux_9cm', 'h_10cm', 'theta_10cm', 'flux_10cm', 'h_15cm', 'theta_15cm', 'flux_15cm', 'h_20cm', 'theta_20cm', 'flux_20cm', 'h_30cm', 'theta_30cm', 'flux_30cm', 'h_50cm', 'theta_50cm', 'flux_50cm'), header = FALSE, nrow=948, skip=11)
  obs_10cm_theta <- obs[,c(1,3,6,9,12,15,18,21,24,27,30,33)]
  obs_0_10cm_theta_avg <- apply(obs_10cm_theta[2:nrow(obs_10cm_theta), 2:12], 1, mean)
  
  #0-10 cm avg
  days_to_traff_0_10cm <- sapply(seq(0.95, 0.7, -0.01), function(x) {
    obs_10cm_theta$time_days[head(which(obs_0_10cm_theta_avg < apply(obs_10cm_theta[1, 3:12], 1, mean) * x), n=1)+1] - flood_end}) #FC estimate ignores upper layer
  
  results <- matrix(data=days_to_traff_0_10cm, nrow = length(days_to_traff_0_10cm), ncol = 1, byrow = TRUE, dimnames = list(paste0('FC_', seq(0.95, 0.7, -0.01)), c('days_to_traff_0_10cm')))
  write.csv(results, file = file.path(summaryDir, paste0(f_path, '_results.csv')), row.names=TRUE)
}
dirnames <- list.dirs('D:/PostDoc/Trafficability/climate_runs', full.names = FALSE, recursive = FALSE)
dirnames <- dirnames[grepl(glob2rx('cell*2005-02-15'), dirnames)]
dirnames
# subDir <- dirnames #143 GB before cleaning; 4.3 GB after)
condense_results <- function(folder) {
  subDir <- folder
  print(subDir)
  if (laptop) {
    resultsDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/HYDRUS_runs'
    summaryDir <- file.path(resultsDir, 'summary')
    summaryDir2 <- file.path(resultsDir, 'summary')} else {
      resultsDir <- 'D:/PostDoc/Trafficability/climate_runs'
      summaryDir <- file.path(resultsDir, subDir, 'summary')
  }
  path_names <- list.dirs(file.path(resultsDir, subDir), full.names = FALSE, recursive = FALSE)
# head(path_names)
# tail(path_names)
  path_names <- path_names[2:length(path_names)]
  if(length(path_names)!=2911){
    stop(print('Incorrect number of directories to parse!'))
  }
  error_runs <- sapply(path_names, function(x) file.exists(file.path(resultsDir, subDir, x, x, 'Error.msg')))
  print(paste('There are', sum(error_runs), 'with error messages.'))
# error_runs <- rep(TRUE, length(path_names))
#if error with atmos.in
# for (i in seq_along(error_runs)) {
#   test <- readLines(file.path(resultsDir, subDir, path_names[i], path_names[i], 'Error.msg'), n=1, skipNul = TRUE)
#   if(test == " Error when reading from an input file Atmosph.in !                             ") {
#     error_runs[i] <- FALSE
#   }
# }
# 
# sum(error_runs) #4
  error_names <- path_names[error_runs]
  path_names <- path_names[!error_runs]
  #write soil names with error.msgs to file
  if(sum(error_runs) > 0) {
    write.csv(data.frame(soilnames=error_names), file.path(resultsDir, 'errors', paste0('H1D_errors_', subDir, '.csv')), row.names = FALSE)
  }
  if(!dir.exists(file.path(resultsDir, subDir, 'summary'))) {
    dir.create(file.path(resultsDir, subDir, 'summary'))
  }
# if(!dir.exists(file.path(resultsDir, 'errors'))) {
#   dir.create(file.path(resultsDir, 'errors'))
# }
  #loop through directories applying function to condense results
  sapply(path_names, determine_trafficability_v2, resultsDir=resultsDir, flood_duration=4, removeFiles=TRUE, subDir=subDir, print_f_path=FALSE, summaryDir=summaryDir)
  print(paste(subDir, 'has been completed.'))
}
condense_results(folder = dirnames[1])
condense_results(folder = dirnames[2])
condense_results(folder = dirnames[3])
condense_results(folder = dirnames[4])
condense_results(folder = dirnames[5])
condense_results(folder = dirnames[6])
which(path_names=='soil_10789752')
# determine_trafficability_v2(f_path = 'Panoche_plants', resultsDir = resultsDir, flood_duration = 4)
