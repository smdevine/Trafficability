library(extrafont)
library(extrafontdb)
#font_import() #only needs to be done one time after updating and re-installing R and moving and updating packages
loadfonts(device = 'win')
dataDir <- 'C:/Users/smdevine/Desktop/PostDoc/trafficability'
FiguresDir <- file.path(dataDir, 'Figures', 'final')
resultsDir <- 'D:/PostDoc/Trafficability/climate_runs'
saveDir <- file.path(resultsDir, 'overall_summaries_final')
CIMISdir <- file.path(resultsDir, 'CIMIS_cell_selection')
if (!dir.exists(file.path(saveDir, 'by_cell'))) {
  dir.create(file.path(saveDir, 'by_cell'))
  dir.create(file.path(saveDir, 'by_cell', 'monthly_means'))
  dir.create(file.path(saveDir, 'by_cell', 'monthly_means', 'summaries'))
}
if (!dir.exists(file.path(saveDir, 'median_climate_summaries'))) {
  dir.create(file.path(saveDir, 'median_climate_summaries'))
}

#power function & exponential functions to be tested
power_func <- function(x, b, z, a) {b*x^z + a}
exp_func <- function(x, b, z, a) {b*exp(z*x) + a}
# power_simp_func <- function(x, b, z) {b*x^z}

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
  if(!dir.exists(file.path(saveDir, 'by_cell', 'by_texture'))) {
    dir.create(file.path(saveDir, 'by_cell', 'by_texture'))
  } 
  if(!dir.exists(file.path(saveDir, 'by_cell', 'by_texture', texture))) {
    dir.create(file.path(saveDir, 'by_cell', 'by_texture', texture))
  } 
  write.csv(climate_stats_by_month, file.path(saveDir, 'by_cell', 'by_texture', texture, paste0(cellname, '_', texture, '_summary.csv')), row.names = TRUE)
}
textural_classes
lapply(list.files(file.path(saveDir, 'by_cell'), pattern = glob2rx('*.csv'), recursive = FALSE), summarize_by_texture, texture='sand')
lapply(list.files(file.path(saveDir, 'by_cell'), pattern = glob2rx('*.csv'), recursive = FALSE), summarize_by_texture, texture='loamy sand')
lapply(list.files(file.path(saveDir, 'by_cell'), pattern = glob2rx('*.csv'), recursive = FALSE), summarize_by_texture, texture='sandy loam')
lapply(list.files(file.path(saveDir, 'by_cell'), pattern = glob2rx('*.csv'), recursive = FALSE), summarize_by_texture, texture='silt loam')
lapply(list.files(file.path(saveDir, 'by_cell'), pattern = glob2rx('*.csv'), recursive = FALSE), summarize_by_texture, texture='loam')
lapply(list.files(file.path(saveDir, 'by_cell'), pattern = glob2rx('*.csv'), recursive = FALSE), summarize_by_texture, texture='sandy clay loam')
lapply(list.files(file.path(saveDir, 'by_cell'), pattern = glob2rx('*.csv'), recursive = FALSE), summarize_by_texture, texture='clay loam')
lapply(list.files(file.path(saveDir, 'by_cell'), pattern = glob2rx('*.csv'), recursive = FALSE), summarize_by_texture, texture='silty clay loam')
lapply(list.files(file.path(saveDir, 'by_cell'), pattern = glob2rx('*.csv'), recursive = FALSE), summarize_by_texture, texture='silty clay')
lapply(list.files(file.path(saveDir, 'by_cell'), pattern = glob2rx('*.csv'), recursive = FALSE), summarize_by_texture, texture='clay')

#now bind to climate data 
ETo_cells_of_interest <- read.csv(file.path(CIMISdir, 'ETo_cells_of_interest.csv'), stringsAsFactors = FALSE)
getETo_daily_mn <- function(cellname, start_date, days, days_to_flood=4) {
  flood_yr <- unlist(strsplit(start_date, '-'))[1]
  df <- ETo_cells_of_interest[, c('dates', cellname)]
  df <- df[which(df$dates==paste0('01_01_', flood_yr)):which(df$dates==paste0('12_31_', flood_yr)),2] #keep as mm per day
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
  results <- results[,c(1:2,4,3)]
  if(!dir.exists(file.path(saveDir, 'by_cell', 'by_texture', 'summaries'))) {
    dir.create(file.path(saveDir, 'by_cell', 'by_texture', 'summaries'))
  }
  write.csv(results, file.path(saveDir, 'by_cell', 'by_texture', 'summaries', paste0('ETo_vs_days_to_traffic_', texture, '_', stat, '.csv')))
  results
}
sand_summary <- bind_CIMIS_to_days_to_traffic(texture = 'sand', stat = 'Median')
loamy_sand_summary <- bind_CIMIS_to_days_to_traffic(texture = 'loamy sand', stat = 'Median')
sandy_loam_summary <- bind_CIMIS_to_days_to_traffic(texture = 'sandy loam', stat = 'Median')
silt_loam_summary <- bind_CIMIS_to_days_to_traffic(texture = 'silt loam', stat = 'Median')
loam_summary <- bind_CIMIS_to_days_to_traffic(texture = 'loam', stat = 'Median')
sandy_clay_loam_summary <- bind_CIMIS_to_days_to_traffic(texture = 'sandy clay loam', stat = 'Median')
clay_loam_summary <- bind_CIMIS_to_days_to_traffic(texture = 'clay loam', stat = 'Median')
silty_clay_loam_summary <- bind_CIMIS_to_days_to_traffic(texture = 'silty clay loam', stat = 'Median')
silty_clay_summary <- bind_CIMIS_to_days_to_traffic(texture = 'silty clay', stat = 'Median')
clay_summary <- bind_CIMIS_to_days_to_traffic(texture = 'clay', stat = 'Median')

#Q1 summaries
sand_q1 <- bind_CIMIS_to_days_to_traffic(texture = 'sand', stat = '1st Qu.')
loamy_sand_q1 <- bind_CIMIS_to_days_to_traffic(texture = 'loamy sand', stat = '1st Qu.')
sandy_loam_q1 <- bind_CIMIS_to_days_to_traffic(texture = 'sandy loam', stat = '1st Qu.')
silt_loam_q1 <- bind_CIMIS_to_days_to_traffic(texture = 'silt loam', stat = '1st Qu.')
loam_q1 <- bind_CIMIS_to_days_to_traffic(texture = 'loam', stat = '1st Qu.')
sandy_clay_loam_q1 <- bind_CIMIS_to_days_to_traffic(texture = 'sandy clay loam', stat = '1st Qu.')
clay_loam_q1 <- bind_CIMIS_to_days_to_traffic(texture = 'clay loam', stat = '1st Qu.')
silty_clay_loam_q1 <- bind_CIMIS_to_days_to_traffic(texture = 'silty clay loam', stat = '1st Qu.')
silty_clay_q1 <- bind_CIMIS_to_days_to_traffic(texture = 'silty clay', stat = '1st Qu.')
clay_q1 <- bind_CIMIS_to_days_to_traffic(texture = 'clay', stat = '1st Qu.')

#3rd Qu.
sand_q3 <- bind_CIMIS_to_days_to_traffic(texture = 'sand', stat = '3rd Qu.')
loamy_sand_q3 <- bind_CIMIS_to_days_to_traffic(texture = 'loamy sand', stat = '3rd Qu.')
sandy_loam_q3 <- bind_CIMIS_to_days_to_traffic(texture = 'sandy loam', stat = '3rd Qu.')
silt_loam_q3 <- bind_CIMIS_to_days_to_traffic(texture = 'silt loam', stat = '3rd Qu.')
loam_q3 <- bind_CIMIS_to_days_to_traffic(texture = 'loam', stat = '3rd Qu.')
sandy_clay_loam_q3 <- bind_CIMIS_to_days_to_traffic(texture = 'sandy clay loam', stat = '3rd Qu.')
clay_loam_q3 <- bind_CIMIS_to_days_to_traffic(texture = 'clay loam', stat = '3rd Qu.')
silty_clay_loam_q3 <- bind_CIMIS_to_days_to_traffic(texture = 'silty clay loam', stat = '3rd Qu.')
silty_clay_q3 <- bind_CIMIS_to_days_to_traffic(texture = 'silty clay', stat = '3rd Qu.')
clay_q3 <- bind_CIMIS_to_days_to_traffic(texture = 'clay', stat = '3rd Qu.')

#function to estimate power curve parameters to fit traffic vs. daily mean ET by textural class, documenting R^2
power_curve_fit <- function(df, texture) {
  nls_fit <- nls(days_to_traffic ~ b*meanET^z+a, data = df, start = list(b=25, z=-0.8, a=3), control=list(warnOnly=TRUE))
  nls_fit_summary <- summary(nls_fit)
  plot(df$meanET, df$days_to_traffic, main=texture)
  curve(power_func(x, b=nls_fit_summary$coefficients[1,1], z=nls_fit_summary$coefficients[2,1], a=nls_fit_summary$coefficients[3,1]), from=min(df$meanET), to=max(df$meanET), lty=2, add=TRUE)
  data.frame(texture_class=texture, b=nls_fit_summary$coefficients[1,1], z=nls_fit_summary$coefficients[2,1], a=nls_fit_summary$coefficients[3,1], R2=cor(predict(nls_fit), df$days_to_traffic)^2, stringsAsFactors = FALSE)
}

power_curve_results <- do.call(rbind, mapply(FUN=power_curve_fit, list(sand_summary, loamy_sand_summary, sandy_loam_summary, silt_loam_summary, loam_summary, sandy_clay_loam_summary, clay_loam_summary, silty_clay_loam_summary, silty_clay_summary, clay_summary), textural_classes, SIMPLIFY = FALSE))
power_curve_results

if(!dir.exists(file.path(saveDir, 'by_cell', 'by_texture', 'summaries', 'math_models'))) {
  dir.create(file.path(saveDir, 'by_cell', 'by_texture', 'summaries', 'math_models'))
}
write.csv(power_curve_results, file.path(saveDir, 'by_cell', 'by_texture', 'summaries', 'math_models', 'power_curve_median.csv'), row.names = FALSE)
write.csv(power_curve_results, file.path(dataDir, 'climate_soil', 'final results', 'power_curve_median.csv'), row.names = FALSE)

power_curve_q1 <- do.call(rbind, mapply(FUN=power_curve_fit, list(sand_q1, loamy_sand_q1, sandy_loam_q1, silt_loam_q1, loam_q1, sandy_clay_loam_q1, clay_loam_q1, silty_clay_loam_q1, silty_clay_q1, clay_q1), textural_classes, SIMPLIFY = FALSE))
write.csv(power_curve_q1, file.path(saveDir, 'by_cell', 'by_texture', 'summaries', 'math_models', 'power_curve_q1.csv'), row.names = FALSE)
write.csv(power_curve_q1, file.path(dataDir, 'climate_soil', 'final results', 'power_curve_q1.csv'), row.names = FALSE)

power_curve_q3 <- do.call(rbind, mapply(FUN=power_curve_fit, list(sand_q3, loamy_sand_q3, sandy_loam_q3, silt_loam_q3, loam_q3, sandy_clay_loam_q3, clay_loam_q3, silty_clay_loam_q3, silty_clay_q3, clay_q3), textural_classes, SIMPLIFY = FALSE))
write.csv(power_curve_q3, file.path(saveDir, 'by_cell', 'by_texture', 'summaries', 'math_models', 'power_curve_q3.csv'), row.names = FALSE)
write.csv(power_curve_q3, file.path(dataDir, 'climate_soil', 'final results', 'power_curve_q3.csv'), row.names = FALSE)

# df <- sandy_loam_summary
# texture <- 'sandy loam'
exp_decay_fit <- function(df, texture) {
  nls_fit <- nls(days_to_traffic ~ b*exp(meanET*z)+a, data = df, start = list(b=50, z=-0.8, a=4), control=list(warnOnly=TRUE))
  nls_fit_summary <- summary(nls_fit)
  plot(df$meanET, df$days_to_traffic, main=texture)
  curve(exp_func(x, b=nls_fit_summary$coefficients[1,1], z=nls_fit_summary$coefficients[2,1], a=nls_fit_summary$coefficients[3,1]), from=min(df$meanET), to=max(df$meanET), lty=2, add=TRUE)
  data.frame(texture_class=texture, b=nls_fit_summary$coefficients[1,1], z=nls_fit_summary$coefficients[2,1], a=nls_fit_summary$coefficients[3,1], R2=cor(predict(nls_fit), df$days_to_traffic)^2, stringsAsFactors = FALSE)
}
exp_decay_fit(loam_summary, 'loam')
exp_decay_results <- do.call(rbind, mapply(FUN=exp_decay_fit, list(sand_summary, loamy_sand_summary, sandy_loam_summary, silt_loam_summary, loam_summary, sandy_clay_loam_summary, clay_loam_summary, silty_clay_loam_summary, silty_clay_summary, clay_summary), textural_classes, SIMPLIFY = FALSE))
write.csv(exp_decay_results, file.path(saveDir, 'by_cell', 'by_texture', 'summaries', 'math_models', 'exp_decay_median.csv'), row.names = FALSE)
write.csv(exp_decay_results, file.path(dataDir, 'climate_soil', 'final results', 'exp_decay_median.csv'), row.names = FALSE)

exp_decay_q1 <- do.call(rbind, mapply(FUN=exp_decay_fit, list(sand_q1, loamy_sand_q1, sandy_loam_q1, silt_loam_q1, loam_q1, sandy_clay_loam_q1, clay_loam_q1, silty_clay_loam_q1, silty_clay_q1, clay_q1), textural_classes, SIMPLIFY = FALSE))
write.csv(exp_decay_q1, file.path(saveDir, 'by_cell', 'by_texture', 'summaries', 'math_models', 'exp_decay_q1.csv'), row.names = FALSE)
write.csv(exp_decay_q1, file.path(dataDir, 'climate_soil', 'final results', 'exp_decay_q1.csv'), row.names = FALSE)

exp_decay_q3 <- do.call(rbind, mapply(FUN=exp_decay_fit, list(sand_q3, loamy_sand_q3, sandy_loam_q3, silt_loam_q3, loam_q3, sandy_clay_loam_q3, clay_loam_q3, silty_clay_loam_q3, silty_clay_q3, clay_q3), textural_classes, SIMPLIFY = FALSE))
write.csv(exp_decay_q3, file.path(saveDir, 'by_cell', 'by_texture', 'summaries', 'math_models', 'exp_decay_q3.csv'), row.names = FALSE)
write.csv(exp_decay_q3, file.path(dataDir, 'climate_soil', 'final results', 'exp_decay_q3.csv'), row.names = FALSE)


#make a figure showing clay, clay loam, sandy loam, and loam

texture_colors <- data.frame(textures=c('clay', 'silty clay' , 'silty clay loam', 'clay loam', 'silt loam', 'sandy clay loam', 'loam', 'sandy loam', 'loamy sand', 'sand'), red=c(169, 0, 0, 223, 0, 170, 255, 230, 115, 255), green=c(0, 112, 197, 115, 168, 255, 0, 152, 76, 255), blue=c(230, 255, 255, 255, 132, 0, 0, 0, 0, 0), stringsAsFactors = FALSE)

power_curve_results <- read.csv()

tiff(file = file.path(FiguresDir, 'days_to_traffic_vs_ETo_all_textures.tif'), family = 'Times New Roman', width = 6.5, height = 3.75, pointsize = 11, units = 'in', res=800, compression = 'lzw')
par(mar=c(3.5, 4.5, 0.5, 0.5))
curve(power_func(x, b=power_curve_results$b[power_curve_results$texture_class=='clay'], z=power_curve_results$z[power_curve_results$texture_class=='clay'], a=power_curve_results$a[power_curve_results$texture_class=='clay']), from=1, to=7, lwd=1.5, col=rgb(red = 169/255, green = 0, blue = 230/255), xlab='', ylab='Days to trafficability after flooding', xlim = c(0.95, 6.95), ylim=c(0,43), xaxt='n')
mtext(text=expression(paste('ETo (mm ', day^-1, ')')), side=1, line=2.5)
axis(side=1, at=1:7,labels = 1:7)
curve(exp_func(x, b=exp_decay_results$b[exp_decay_results$texture_class=='silty clay'], z=exp_decay_results$z[exp_decay_results$texture_class=='silty clay'], a=exp_decay_results$a[exp_decay_results$texture_class=='silty clay']), from=1, to=7, lwd=1.5, col=rgb(red = 0/255, green = 112/255, blue = 255/255), add=TRUE)
for (i in 3:10) {
  curve(power_func(x, b=power_curve_results$b[power_curve_results$texture_class==texture_colors$textures[i]], z=power_curve_results$z[power_curve_results$texture_class==texture_colors$textures[i]], a=power_curve_results$a[power_curve_results$texture_class==texture_colors$textures[i]]), from=1, to=7, lwd=1.5, col=rgb(red = texture_colors$red[i]/255, green = texture_colors$green[i]/255, blue = texture_colors$blue[i]/255), add=TRUE)
}
legend('topright', legend=texture_colors$textures, col=rgb(red=texture_colors$red/255, green = texture_colors$green/255, blue = texture_colors$blue/255), lty=1, lwd=1.5, ncol=2)
<<<<<<< HEAD
legend('topleft', legend='a', bty='n')
=======
legend('topleft', legend='a', bty='n', inset = c(-0.025, 0))
>>>>>>> 89925259eb3b1234bc8c35d172fa0e8e5a3637b7
# text(x=1.35, y=40, 'A')
dev.off()

#make a figure showing IQRs for clay, clay loam, loam, and sandy loam
tiff(file = file.path(FiguresDir, 'days_to_traffic_vs_ETo_all_IQRs_4textures.tif'), family = 'Times New Roman', width = 6.5, height = 3.75, pointsize = 11, units = 'in', res=800, compression = 'lzw')
par(mar=c(3.5, 4.5, 0.5, 0.5))
curve(power_func(x, b=power_curve_results$b[power_curve_results$texture_class=='clay loam'], z=power_curve_results$z[power_curve_results$texture_class=='clay loam'], a=power_curve_results$a[power_curve_results$texture_class=='clay loam']), from=1.1, to=6.8, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='clay loam',2:4]/255), xlab='', ylab='Days to trafficability after flooding', xlim = c(1.1, 6.8), ylim=c(0,43), xaxt='n')
mtext(text=expression(paste('ETo (mm ', day^-1, ')')), side=1, line=2.5)
axis(side=1, at=1:7,labels = 1:7)
curve(power_func(x, b=power_curve_q1$b[power_curve_q1$texture_class=='clay loam'], z=power_curve_q1$z[power_curve_q1$texture_class=='clay loam'], a=power_curve_q1$a[power_curve_q1$texture_class=='clay loam']), from=1.1, to=6.8, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='clay loam',2:4]/255), add=TRUE, lty=2)
curve(power_func(x, b=power_curve_q3$b[power_curve_q3$texture_class=='clay loam'], z=power_curve_q3$z[power_curve_q3$texture_class=='clay loam'], a=power_curve_q3$a[power_curve_q3$texture_class=='clay loam']), from=1.1, to=6.8, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='clay loam',2:4]/255), add=TRUE, lty=2)

curve(power_func(x, b=power_curve_results$b[power_curve_results$texture_class=='loam'], z=power_curve_results$z[power_curve_results$texture_class=='loam'], a=power_curve_results$a[power_curve_results$texture_class=='loam']), from=1.1, to=6.8, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='loam',2:4]/255), add=TRUE)
curve(power_func(x, b=power_curve_q1$b[power_curve_q1$texture_class=='loam'], z=power_curve_q1$z[power_curve_q1$texture_class=='loam'], a=power_curve_q1$a[power_curve_q1$texture_class=='loam']), from=1.1, to=6.8, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='loam',2:4]/255), add=TRUE, lty=2)
curve(power_func(x, b=power_curve_q3$b[power_curve_q3$texture_class=='loam'], z=power_curve_q3$z[power_curve_q3$texture_class=='loam'], a=power_curve_q3$a[power_curve_q3$texture_class=='loam']), from=1.1, to=6.8, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='loam',2:4]/255), add=TRUE, lty=2)

curve(power_func(x, b=power_curve_results$b[power_curve_results$texture_class=='sandy loam'], z=power_curve_results$z[power_curve_results$texture_class=='sandy loam'], a=power_curve_results$a[power_curve_results$texture_class=='sandy loam']), from=1.1, to=6.8, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='sandy loam',2:4]/255), add=TRUE)
curve(power_func(x, b=power_curve_q1$b[power_curve_q1$texture_class=='sandy loam'], z=power_curve_q1$z[power_curve_q1$texture_class=='sandy loam'], a=power_curve_q1$a[power_curve_q1$texture_class=='sandy loam']), from=1.1, to=6.8, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='sandy loam',2:4]/255), add=TRUE, lty=2)
curve(power_func(x, b=power_curve_q3$b[power_curve_q3$texture_class=='sandy loam'], z=power_curve_q3$z[power_curve_q3$texture_class=='sandy loam'], a=power_curve_q3$a[power_curve_q3$texture_class=='sandy loam']), from=1.1, to=6.8, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='sandy loam',2:4]/255), add=TRUE, lty=2)

curve(power_func(x, b=power_curve_results$b[power_curve_results$texture_class=='clay'], z=power_curve_results$z[power_curve_results$texture_class=='clay'], a=power_curve_results$a[power_curve_results$texture_class=='clay']), from=1.1, to=6.8, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='clay',2:4]/255), add=TRUE)
curve(power_func(x, b=power_curve_q1$b[power_curve_q1$texture_class=='clay'], z=power_curve_q1$z[power_curve_q1$texture_class=='clay'], a=power_curve_q1$a[power_curve_q1$texture_class=='clay']), from=1.1, to=6.8, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='clay',2:4]/255), add=TRUE, lty=2)
curve(exp_func(x, b=exp_decay_q3$b[exp_decay_q3$texture_class=='clay'], z=exp_decay_q3$z[exp_decay_q3$texture_class=='clay'], a=exp_decay_q3$a[exp_decay_q3$texture_class=='clay']), from=1.1, to=6.8, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='clay',2:4]/255), add=TRUE, lty=2)
legend('topright', legend=c('clay', 'clay loam', 'loam', 'sandy loam'), col=rgb(texture_colors[which(texture_colors$textures %in% c('clay', 'clay loam', 'loam', 'sandy loam')),2:4]/255), lty=1, lwd=1.5, ncol=2)
legend('topleft', legend='b', bty='n', inset = c(-0.025, 0))
# text(x=1.5, y=43, 'B')
dev.off()

#4 texture plot with points
#make a figure showing IQRs for clay, clay loam, loam, and sandy loam
tiff(file = file.path(FiguresDir, 'days_to_traffic_vs_ETo_all_IQRs_4textures_pts.tif'), family = 'Times New Roman', width = 6.5, height = 3.75, pointsize = 11, units = 'in', res=800, compression = 'lzw')
par(mar=c(3.5, 4.5, 0.5, 0.5))
curve(power_func(x, b=power_curve_results$b[power_curve_results$texture_class=='clay loam'], z=power_curve_results$z[power_curve_results$texture_class=='clay loam'], a=power_curve_results$a[power_curve_results$texture_class=='clay loam']), from=0.89, to=6.37, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='clay loam',2:4]/255), xlab='', ylab='Days to trafficability after flooding', xlim = c(0.8, 6.8), ylim=c(0,45), xaxt='n')
points(clay_loam_summary$meanET, clay_loam_summary$days_to_traffic, pch=1, col=rgb(texture_colors[texture_colors$textures=='clay loam',2:4]/255), cex=0.7)
mtext(text=expression(paste('ETo (mm ', day^-1, ')')), side=1, line=2.5)
axis(side=1, at=1:7,labels = 1:7)
curve(power_func(x, b=power_curve_q1$b[power_curve_q1$texture_class=='clay loam'], z=power_curve_q1$z[power_curve_q1$texture_class=='clay loam'], a=power_curve_q1$a[power_curve_q1$texture_class=='clay loam']), from=0.89, to=6.37, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='clay loam',2:4]/255), add=TRUE, lty=2)
curve(power_func(x, b=power_curve_q3$b[power_curve_q3$texture_class=='clay loam'], z=power_curve_q3$z[power_curve_q3$texture_class=='clay loam'], a=power_curve_q3$a[power_curve_q3$texture_class=='clay loam']), from=0.89, to=6.37, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='clay loam',2:4]/255), add=TRUE, lty=2)

curve(power_func(x, b=power_curve_results$b[power_curve_results$texture_class=='loam'], z=power_curve_results$z[power_curve_results$texture_class=='loam'], a=power_curve_results$a[power_curve_results$texture_class=='loam']), from=0.96, to=6.83, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='loam',2:4]/255), add=TRUE)
points(loam_summary$meanET, loam_summary$days_to_traffic, pch=1, col=rgb(texture_colors[texture_colors$textures=='loam',2:4]/255), cex=0.7)
curve(power_func(x, b=power_curve_q1$b[power_curve_q1$texture_class=='loam'], z=power_curve_q1$z[power_curve_q1$texture_class=='loam'], a=power_curve_q1$a[power_curve_q1$texture_class=='loam']), from=0.96, to=6.83, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='loam',2:4]/255), add=TRUE, lty=2)
curve(power_func(x, b=power_curve_q3$b[power_curve_q3$texture_class=='loam'], z=power_curve_q3$z[power_curve_q3$texture_class=='loam'], a=power_curve_q3$a[power_curve_q3$texture_class=='loam']), from=0.96, to=6.83, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='loam',2:4]/255), add=TRUE, lty=2)

curve(power_func(x, b=power_curve_results$b[power_curve_results$texture_class=='sandy loam'], z=power_curve_results$z[power_curve_results$texture_class=='sandy loam'], a=power_curve_results$a[power_curve_results$texture_class=='sandy loam']), from=0.799, to=6.74, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='sandy loam',2:4]/255), add=TRUE)
points(sandy_loam_summary$meanET, sandy_loam_summary$days_to_traffic, pch=1, col=rgb(texture_colors[texture_colors$textures=='sandy loam',2:4]/255), cex=0.7)
curve(power_func(x, b=power_curve_q1$b[power_curve_q1$texture_class=='sandy loam'], z=power_curve_q1$z[power_curve_q1$texture_class=='sandy loam'], a=power_curve_q1$a[power_curve_q1$texture_class=='sandy loam']), from=0.799, to=6.74, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='sandy loam',2:4]/255), add=TRUE, lty=2)
curve(power_func(x, b=power_curve_q3$b[power_curve_q3$texture_class=='sandy loam'], z=power_curve_q3$z[power_curve_q3$texture_class=='sandy loam'], a=power_curve_q3$a[power_curve_q3$texture_class=='sandy loam']), from=0.799, to=6.74, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='sandy loam',2:4]/255), add=TRUE, lty=2)

curve(power_func(x, b=power_curve_results$b[power_curve_results$texture_class=='clay'], z=power_curve_results$z[power_curve_results$texture_class=='clay'], a=power_curve_results$a[power_curve_results$texture_class=='clay']), from=1.07, to=6.82, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='clay',2:4]/255), add=TRUE)
points(clay_summary$meanET, clay_summary$days_to_traffic, pch=1, col=rgb(texture_colors[texture_colors$textures=='clay',2:4]/255), cex=0.7)
curve(power_func(x, b=power_curve_q1$b[power_curve_q1$texture_class=='clay'], z=power_curve_q1$z[power_curve_q1$texture_class=='clay'], a=power_curve_q1$a[power_curve_q1$texture_class=='clay']), from=1.07, to=6.82, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='clay',2:4]/255), add=TRUE, lty=2)
curve(power_func(x, b=power_curve_q3$b[power_curve_q3$texture_class=='clay'], z=power_curve_q3$z[power_curve_q3$texture_class=='clay'], a=power_curve_q3$a[power_curve_q3$texture_class=='clay']), from=1.07, to=6.82, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='clay',2:4]/255), add=TRUE, lty=2)
legend('topright', legend=c('clay', 'clay loam', 'loam', 'sandy loam'), col=rgb(texture_colors[which(texture_colors$textures %in% c('clay', 'clay loam', 'loam', 'sandy loam')),2:4]/255), lty=1, lwd=1.5, ncol=2)
text(x=1.35, y=43, 'B')
dev.off()

#3 texture plot with points
#make a figure showing IQRs for clay, loam, and sandy loam, top 3 surface textures across study area
tiff(file = file.path(FiguresDir, 'days_to_traffic_vs_ETo_all_IQRs_3textures_pts.tif'), family = 'Times New Roman', width = 6.5, height = 3.75, pointsize = 11, units = 'in', res=800, compression = 'lzw')
par(mar=c(3.5, 4.5, 0.5, 0.5))
curve(power_func(x, b=power_curve_results$b[power_curve_results$texture_class=='loam'], z=power_curve_results$z[power_curve_results$texture_class=='loam'], a=power_curve_results$a[power_curve_results$texture_class=='loam']), from=0.96, to=6.98, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='loam',2:4]/255),  xlab='', ylab='Days to trafficability after flooding', xlim = c(0.95, 6.95), ylim=c(0,43), xaxt='n')
points(loam_summary$meanET, loam_summary$days_to_traffic, pch=1, col=rgb(texture_colors[texture_colors$textures=='loam',2:4]/255), cex=0.7)
mtext(text=expression(paste('ETo (mm ', day^-1, ')')), side=1, line=2.5)
axis(side=1, at=1:7,labels = 1:7)
curve(power_func(x, b=power_curve_q1$b[power_curve_q1$texture_class=='loam'], z=power_curve_q1$z[power_curve_q1$texture_class=='loam'], a=power_curve_q1$a[power_curve_q1$texture_class=='loam']), from=0.96, to=6.98, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='loam',2:4]/255), add=TRUE, lty=2)
curve(power_func(x, b=power_curve_q3$b[power_curve_q3$texture_class=='loam'], z=power_curve_q3$z[power_curve_q3$texture_class=='loam'], a=power_curve_q3$a[power_curve_q3$texture_class=='loam']), from=0.96, to=6.98, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='loam',2:4]/255), add=TRUE, lty=2)

curve(power_func(x, b=power_curve_results$b[power_curve_results$texture_class=='sandy loam'], z=power_curve_results$z[power_curve_results$texture_class=='sandy loam'], a=power_curve_results$a[power_curve_results$texture_class=='sandy loam']), from=0.778, to=6.74, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='sandy loam',2:4]/255), add=TRUE)
points(sandy_loam_summary$meanET, sandy_loam_summary$days_to_traffic, pch=1, col=rgb(texture_colors[texture_colors$textures=='sandy loam',2:4]/255), cex=0.7)
curve(power_func(x, b=power_curve_q1$b[power_curve_q1$texture_class=='sandy loam'], z=power_curve_q1$z[power_curve_q1$texture_class=='sandy loam'], a=power_curve_q1$a[power_curve_q1$texture_class=='sandy loam']), from=0.778, to=6.74, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='sandy loam',2:4]/255), add=TRUE, lty=2)
curve(power_func(x, b=power_curve_q3$b[power_curve_q3$texture_class=='sandy loam'], z=power_curve_q3$z[power_curve_q3$texture_class=='sandy loam'], a=power_curve_q3$a[power_curve_q3$texture_class=='sandy loam']), from=0.778, to=6.74, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='sandy loam',2:4]/255), add=TRUE, lty=2)

curve(power_func(x, b=power_curve_results$b[power_curve_results$texture_class=='clay'], z=power_curve_results$z[power_curve_results$texture_class=='clay'], a=power_curve_results$a[power_curve_results$texture_class=='clay']), from=0.94, to=6.36, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='clay',2:4]/255), add=TRUE)
points(clay_summary$meanET, clay_summary$days_to_traffic, pch=1, col=rgb(texture_colors[texture_colors$textures=='clay',2:4]/255), cex=0.7)
curve(power_func(x, b=power_curve_q1$b[power_curve_q1$texture_class=='clay'], z=power_curve_q1$z[power_curve_q1$texture_class=='clay'], a=power_curve_q1$a[power_curve_q1$texture_class=='clay']), from=0.94, to=6.36, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='clay',2:4]/255), add=TRUE, lty=2)
curve(power_func(x, b=power_curve_q3$b[power_curve_q3$texture_class=='clay'], z=power_curve_q3$z[power_curve_q3$texture_class=='clay'], a=power_curve_q3$a[power_curve_q3$texture_class=='clay']), from=0.94, to=6.36
, lwd=1.5, col=rgb(texture_colors[texture_colors$textures=='clay',2:4]/255), add=TRUE, lty=2)
legend('topright', legend=c('clay','loam','sandy loam'), col=rgb(texture_colors[which(texture_colors$textures %in% c('clay','loam', 'sandy loam')),2:4]/255), lty=1, lwd=1.5, ncol=2)
legend('topleft', legend='b', bty='n', inset = c(-0.025, 0))
# text(x=1.35, y=40, 'B')
dev.off()