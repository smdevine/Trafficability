dataDir <- 'C:/Users/smdevine/Desktop/PostDoc/trafficability'
resultsDir <- 'D:/PostDoc/Trafficability/climate_runs'
errorDir <- file.path(resultsDir, 'climate_reruns')
saveDir <- file.path(resultsDir, 'overall_summaries')
errorDir2 <- file.path(resultsDir, 'errors_take2')
no_data <- read.csv(file.path(errorDir2, 'no_data.csv'), row.names = 1, stringsAsFactors = FALSE)
traffDefs <- read.csv(file.path(dataDir, 'trafficability_defs_11_3_20.csv'), stringsAsFactors = FALSE)
theta_seq <- as.character(seq(0.95, 0.7, -0.01)) ##how results are arranged by decreasing theta_fc thresholds
extract_0_10cm_trafficability <- function(x, dirPath) {
  results <- read.csv(file.path(dirPath, x), row.names=1, na.strings = 'numeric(0)')
}
soil_df <- read.csv(file.path(dataDir, 'soils_of_interest', 'soils_modeled_revised_QCpass_Oct2020.csv'), stringsAsFactors = FALSE)
months <- c('01', '02', '03', '04')
years <- c('2005', '2009', '2015')
subDirs <- expand.grid(years, months)
subDirs <- apply(subDirs, 1, function(x) paste(x, collapse = '-'))
subDirs <- paste0(subDirs, '-15')

i <- 1
for (i in 2:length(subDirs)) {
  subDirs2 <- list.dirs(file.path(resultsDir, subDirs[i]), full.names = FALSE, recursive = FALSE)
  results_template <- soil_df[,c('cokey', 'textural_class', 'traff_def_opt2')]
  colnames(results_template)[3] <- 'traff_def' 
  results_template$cokey <- as.character(results_template$cokey)
# j <- 2
  for (j in 1:length(subDirs2)) {
    print(j)
    summaryDir <- file.path(resultsDir, subDirs[i], subDirs2[j], 'summary') 
    fnames_results <- list.files(summaryDir, full.names = FALSE, recursive = FALSE)
    cokeys_modeled <- sub(pattern = 'soil_', replacement = '', x = fnames_results)
    cokeys_modeled <- sub(pattern = '_results.csv', replacement = '', x = cokeys_modeled)
    trafficability <- lapply(fnames_results, extract_0_10cm_trafficability, dirPath=summaryDir)
    names(trafficability) <- cokeys_modeled
    if(length(cokeys_modeled)!=2911){
      error_fnames_results <- list.files(file.path(errorDir, subDirs2[j], 'summary'), full.names = FALSE, recursive = FALSE)
      error_cokeys_modeled <- sub(pattern = 'soil_', replacement = '', x = error_fnames_results)
      error_cokeys_modeled <- sub(pattern = '_results.csv', replacement = '', x = error_cokeys_modeled)
      trafficability_reruns <- lapply(error_fnames_results, extract_0_10cm_trafficability, dirPath=file.path(errorDir, subDirs2[j], 'summary'))
      names(trafficability_reruns) <- error_cokeys_modeled
      trafficability <- c(trafficability, trafficability_reruns)
      if(length(trafficability)!=2911) {
        errors_fnames_final <- list.files(errorDir2)
        errors_fnames_final <- errors_fnames_final[grepl(subDirs2[j], errors_fnames_final)]
        errors_soils_final <- read.csv(file.path(errorDir2, errors_fnames_final))
        errors_soils_final <- errors_soils_final$soilnames
        errors_soils_final <- sub(pattern = 'soil_', replacement = '', x = errors_soils_final)
        for (k in length(errors_soils_final)) {
          trafficability <- c(trafficability, list(no_data))
          names(trafficability)[length(trafficability)] <- errors_soils_final[k]
        }
        if(length(trafficability)!=2911) {
          stop(print('Missing cokeys from results'))
        }
      }
      trafficability <- trafficability[order(names(trafficability))]
    }
    match_test <- sapply(seq_along(results_template$cokey), function(a) names(trafficability)[a] == results_template$cokey[a])
    if(!all(match_test)){
      stop(print('Data out of order'))
    }
    cellname <- gsub(paste0('_', subDirs[i]), '', subDirs2[j])
    results_template[[cellname]] <- sapply(seq_along(results_template$traff_def), function(a) {trafficability[[a]][as.character(results_template$traff_def[a])==theta_seq,]})
  }
  subDirs[i]
  write.csv(results_template, file.path(saveDir, paste0('trafficability_results_', subDirs[i], '.csv')), row.names=FALSE)
}
subDirs[i]
subDirs2[j]
