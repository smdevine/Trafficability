dataDir <- 'C:/Users/smdevine/Desktop/PostDoc/trafficability'
resultsDir <- 'D:/PostDoc/Trafficability/climate_runs'
saveDir <- file.path(resultsDir, 'overall_summaries')
CIMISdir <- file.path(resultsDir, 'CIMIS_cell_selection')
fnames <- list.files(saveDir)
fnames <- fnames[grepl(glob2rx('*.csv'), fnames)]
dates <- gsub('trafficability_results_', '', fnames)
dates <- gsub('.csv', '', dates)
results_template <- read.csv(file.path(saveDir, fnames[1]), stringsAsFactors = FALSE)
cellnames <- colnames(results_template)[4:14]
results_template <- results_template[,1:3]
results_all <- lapply(seq_along(fnames), function(i) read.csv(file.path(saveDir, fnames[i]), stringsAsFactors = FALSE))
names(results_all) <- dates
for (i in seq_along(cellnames)) {
  results <- results_template
  for (j in seq_along(results_all)) {
    results[[paste0('fd_', dates[j])]] <- results_all[[j]][[cellnames[i]]]
  }
  write.csv(results, file.path(saveDir, 'by_cell', paste0(cellnames[i], '_results.csv')), row.names = FALSE)
}
list.files(file.path(saveDir, 'by_cell'))
median_climate_results <- read.csv(file.path(saveDir, 'by_cell', 'cell_130212_results.csv'), stringsAsFactors = FALSE)
colnames(median_climate_results)
median_climate_stats_by_date <- do.call(cbind, lapply(median_climate_results[,4:15], function(x) summary(x)))

write.csv(median_climate_stats_by_date, file.path(saveDir, 'median_climate_summaries', 'overall_stats_by_flood_date.csv'), row.names=TRUE)
textural_classes <- unique(median_climate_results$textural_class)
textural_classes
sapply(textural_classes, function(x) sum(median_climate_results$textural_class==x))
textural_classes <- textural_classes[textural_classes != 'sandy clay']
textural_classes <- textural_classes[c(6,3,2,7,5,8,1,9,10,4)]
textural_classes
for (i in seq_along(textural_classes)) {
  climate_stats_by_date <- do.call(cbind, lapply(median_climate_results[median_climate_results$textural_class==textural_classes[i],4:15], function(x) { 
    result <- summary(x) 
    if(length(result)==7) {
      result
    } else {c(result, "NA's"=0)}}))
  write.csv(climate_stats_by_date, file.path(saveDir, 'median_climate_summaries', paste0(textural_classes[i],'_stats_by_flood_date.csv')), row.names=TRUE)
}

list.files(file.path(saveDir, 'by_cell'), pattern = glob2rx('*.csv'), recursive = FALSE)

monthly_mean_summarize <- function(fname, month, month_name) {
  cellname <- gsub('_results.csv', '', fname)
  df <- read.csv(file.path(saveDir, 'by_cell', fname), stringsAsFactors = FALSE)
  monthly_avg <- apply(df[,grepl(month, colnames(median_climate_results))], 1, mean, na.rm=FALSE)
  df_monthly_avg <- cbind(df[,1:3], monthly_avg)
  colnames(df_monthly_avg)[4] <- paste0(month_name, '_mean_days_to_traff')
  write.csv(df_monthly_avg, file.path(saveDir, 'by_cell', 'monthly_means', paste0(cellname, '_', month_name, '.csv')), row.names = FALSE)
  climate_stats_by_month <- do.call(cbind, lapply(seq_along(textural_classes), function(i) {
    result <- summary(df_monthly_avg[df_monthly_avg$textural_class==textural_classes[i],4])
    if(length(result)==7) {
      result
    } else {result <- c(result, "NA's"=0)}
    result <- c(result, 'n'=length(df_monthly_avg[df_monthly_avg$textural_class==textural_classes[i],4]))}))
  colnames(climate_stats_by_month) <- textural_classes
  climate_stats_by_month <- cbind(all_textures=c(summary(df_monthly_avg[,4]), 'n'=length(df_monthly_avg[,4])), climate_stats_by_month)
  write.csv(climate_stats_by_month, file.path(saveDir, 'by_cell', 'monthly_means', 'summaries', paste0(cellname, '_', month_name, '_summary.csv')), row.names = TRUE)
}
sapply(list.files(file.path(saveDir, 'by_cell'), pattern = glob2rx('*.csv'), recursive = FALSE), monthly_mean_summarize, month='01.15', month_name='Jan')
sapply(list.files(file.path(saveDir, 'by_cell'), pattern = glob2rx('*.csv'), recursive = FALSE), monthly_mean_summarize, month='02.15', month_name='Feb')
sapply(list.files(file.path(saveDir, 'by_cell'), pattern = glob2rx('*.csv'), recursive = FALSE), monthly_mean_summarize, month='03.15', month_name='Mar')
sapply(list.files(file.path(saveDir, 'by_cell'), pattern = glob2rx('*.csv'), recursive = FALSE), monthly_mean_summarize, month='04.15', month_name='Apr')

summarize_by_texture <- function(fname, texture) {
  cellname <- gsub('_results.csv', '', fname)
  df <- read.csv(file.path(saveDir, 'by_cell', fname), stringsAsFactors = FALSE)
  climate_stats_by_month <- do.call(cbind, lapply(df[,4:ncol(df)], function(x) {
    result <- summary(x[df$textural_class==texture])
    if(length(result)==7) {
      result
    } else {result <- c(result, "NA's"=0)}
    result <- c(result, 'n'=length(x[df$textural_class==texture]))}))
  colnames(climate_stats_by_month) <- colnames(df)[4:ncol(df)]
  if(!dir.exists(file.path(saveDir, 'by_cell', 'by_texture', texture))) {
    dir.create(file.path(saveDir, 'by_cell', 'by_texture', texture))
  } 
  write.csv(climate_stats_by_month, file.path(saveDir, 'by_cell', 'by_texture', texture, paste0(cellname, '_', texture, '_summary.csv')), row.names = TRUE)
}
lapply(list.files(file.path(saveDir, 'by_cell'), pattern = glob2rx('*.csv'), recursive = FALSE), summarize_by_texture, texture='sandy loam')
lapply(list.files(file.path(saveDir, 'by_cell'), pattern = glob2rx('*.csv'), recursive = FALSE), summarize_by_texture, texture='clay')
lapply(list.files(file.path(saveDir, 'by_cell'), pattern = glob2rx('*.csv'), recursive = FALSE), summarize_by_texture, texture='loam')

#now bind to climate data 
ETo_cells_of_interest <- read.csv(file.path(CIMISdir, 'ETo_cells_of_interest.csv'), stringsAsFactors = FALSE)
getETo_daily_mn <- function(cellname, start_date, days, days_to_flood=4) {
  flood_yr <- unlist(strsplit(start_date, '-'))[1]
  df <- ETo_cells_of_interest[, c('dates', cellname)]
  df <- df[which(df$dates==paste0('01_01_', flood_yr)):which(df$dates==paste0('12_31_', flood_yr)),2] / 10
  mean(df[(as.integer(format(as.Date(start_date), format='%j'))+days_to_flood):(as.integer(format(as.Date(start_date), format='%j'))+days+days_to_flood)])
}

# texture <- 'sandy loam'
# stat <- 'Median'
bind_CIMIS_to_days_to_traffic <- function(texture, stat) {
  fnames <- list.files(file.path(saveDir, 'by_cell', 'by_texture', texture))
  cellnames <- gsub(paste0('_', texture, '_summary.csv'), '', fnames)
  results <- data.frame(cellname=as.character(sapply(cellnames, function(x) rep(x, times=12), simplify = TRUE)), date=NA, days_to_traffic=NA, meanET=NA, stringsAsFactors = FALSE)
  results$days_to_traffic <- do.call(c, sapply(seq_along(fnames), function(i) {
    df <- read.csv(file.path(saveDir, 'by_cell', 'by_texture', texture, fnames[i]), stringsAsFactors=FALSE, row.names=1)
    df[stat,]}))
  results$date <- do.call(c, lapply(seq_along(fnames), function(i) {
    df <- read.csv(file.path(saveDir, 'by_cell', 'by_texture', texture, fnames[i]), stringsAsFactors=FALSE, row.names=1)
    colnames(df)}))
  results$date <- gsub('fd_', '', results$date)
  results$date <- gsub('[.]', '-', results$date)
  results$meanET <- mapply(FUN=getETo_daily_mn, cellname=results$cellname, start_date=results$date, days=ceiling(results$days_to_traffic))
  
  write.csv(results, file.path(saveDir, 'by_cell', 'by_texture', 'summaries', paste0('ETo_vs_days_to_traffic_', texture, '_', stat, '.csv')))
  results
}
sandy_loam_summary <- bind_CIMIS_to_days_to_traffic(texture = 'sandy loam', stat = 'Median')
clay_summary <- bind_CIMIS_to_days_to_traffic(texture = 'clay', stat = 'Median')
loam_summary <- bind_CIMIS_to_days_to_traffic(texture = 'loam', stat = 'Median')
plot(sandy_loam_summary$meanET, sandy_loam_summary$days_to_traffic)
plot(clay_summary$meanET, clay_summary$days_to_traffic)
plot(loam_summary$meanET, loam_summary$days_to_traffic)

lapply(results_df[,4:14], function(x) tapply(x, results_df$textural_class, summary))

lapply(results_df[,4:14], function(x) summary(x))
lapply(results_df[,4:14], function(x) tapply(x, results_df$textural_class, summary))



#old functions from summarize_results_v2.R draft analysis
summarize_by_texture <- function(traff_def) {
  textural_class_summary <- do.call(cbind, lapply(unique(soils_modeled_10cm$textural_class), function(x) as.matrix(summary(soils_modeled_10cm[[traff_def]][soils_modeled_10cm$textural_class==x]))[1:6,]))
  colnames(textural_class_summary) <- unique(soils_modeled_10cm$textural_class)
  textural_class_summary <- rbind(textural_class_summary, n=as.integer(sapply(unique(soils_modeled_10cm$textural_class), function(x) sum(!is.na(soils_modeled_10cm[[traff_def]][soils_modeled_10cm$textural_class==x])))))
  textural_class_summary <- rbind(textural_class_summary, mean_clay=sapply(unique(soils_modeled_10cm$textural_class), function(x) mean(soils_modeled_10cm$clay_10cm[soils_modeled_10cm$textural_class==x], na.rm=TRUE)))
  textural_class_summary <- rbind(textural_class_summary, Not_Determined=as.integer(sapply(unique(soils_modeled_10cm$textural_class), function(x) sum(is.na(soils_modeled_10cm[[traff_def]][soils_modeled_10cm$textural_class==x])))))
  colnames(textural_class_summary)[order(textural_class_summary[8,])]
  textural_class_summary <- textural_class_summary[,order(textural_class_summary[8,])]
  print(textural_class_summary)
  write.csv(textural_class_summary, file.path(tablesDir, paste0('summary_by_textural_class_0_10cm_', traff_def, '_unique.csv')), row.names = TRUE)
}
summarize_by_texture('result_opt1')


summarize_by_names <- function(traff_def) {
  soilname_summary <- do.call(cbind, lapply(names(soil_name_mas30), function(x) as.matrix(summary(soils_modeled_10cm[[traff_def]][soils_modeled_10cm$soil_name==x]))[1:6,]))
  colnames(soilname_summary) <- names(soil_name_mas30)
  soilname_summary <- rbind(soilname_summary, n=as.integer(sapply(names(soil_name_mas30), function(x) sum(!is.na(soils_modeled_10cm[[traff_def]][soils_modeled_10cm$soil_name==x])))))
  soilname_summary <- rbind(soilname_summary, mean_clay=sapply(names(soil_name_mas30), function(x) mean(soils_modeled_10cm$clay_10cm[soils_modeled_10cm$soil_name==x], na.rm=TRUE)))
  soilname_summary <- rbind(soilname_summary, texture=sapply(names(soil_name_mas30), function(x) paste(unique(soils_modeled_10cm$textural_class[soils_modeled_10cm$soil_name==x]), collapse = ', ')))
  print(soilname_summary)
  write.csv(soilname_summary, file.path(tablesDir, paste0('summary_by_soilname_0_10cm_properties_', traff_def, '_rosetta5.csv')), row.names = TRUE)
}
summarize_by_names('result_opt1')
