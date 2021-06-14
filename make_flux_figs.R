#calculate normalized drawdown for each soil
#check silt loam median profile
#check silty clay median profile
#get April 2009 results as example
library(extrafont)
library(extrafontdb)
# font_import() #only needs to be done one time after updating and re-installing R and moving and updating packages
loadfonts(device = 'win')
FiguresDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/Figures'
FebResultsDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/climate_runs/typical profiles/cell_181807_2009-02-15'
AprResultsDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/climate_runs/typical profiles/cell_181807_2009-04-15'
trafficDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/climate_run_summaries_final/by_cell/'
soilsDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/soils_of_interest'
textures <- c('clay', 'silty clay' , 'silty clay loam', 'clay loam', 'silt loam', 'sandy clay loam', 'loam', 'sandy loam', 'loamy sand', 'sand')
texture_colors <- data.frame(textures=c('clay', 'silty clay' , 'silty clay loam', 'clay loam', 'silt loam', 'sandy clay loam', 'loam', 'sandy loam', 'loamy sand', 'sand'), texture_labs=c('clay', 'silty\nclay' , 'silty\nclay\nloam', 'clay\nloam', 'silt\nloam', 'sandy\nclay\nloam', 'loam', 'sandy\nloam', 'loamy\nsand', 'sand'), red=c(169, 0, 0, 223, 0, 170, 255, 230, 115, 255), green=c(0, 112, 197, 115, 168, 255, 0, 152, 76, 255), blue=c(230, 255, 255, 255, 132, 0, 0, 0, 0, 0), median_cokeys=c(11694392, 11681, 12763923, 11731790, 12769121, 11693996, 11855610, 11834554, 12350998, 10807468), stringsAsFactors = FALSE)
soilsDF <- read.csv(file.path(soilsDir, 'soils_modeled_revised_QCpass_Oct2020.csv'), stringsAsFactors = FALSE)
trafficDF <- read.csv(file.path(trafficDir, 'cell_181807_results.csv'), stringsAsFactors = FALSE)
trafficDF[trafficDF$cokey==11855610,]

# texture <- 'Sandy loam'
# directory <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/HYDRUS_runs/zero_flux_tests/sandy loam'
theta_fc_est <- function(theta_s, theta_r, Ks, n) {
  (n^(-0.6*(2+log10(Ks))))*(theta_s - theta_r) + theta_r
}
read_in_flux_files <- function(texture, directory, climateRuns=TRUE, columnNames=c('time_days', 'h_0cm', 'theta_0cm', 'flux_0cm', 'h_1cm', 'theta_1cm', 'flux_1cm', 'h_2cm', 'theta_2cm', 'flux_2cm', 'h_3cm', 'theta_3cm', 'flux_3cm', 'h_4cm', 'theta_4cm', 'flux_4cm', 'h_5cm', 'theta_5cm', 'flux_5cm', 'h_6cm', 'theta_6cm', 'flux_6cm', 'h_7cm', 'theta_7cm', 'flux_7cm', 'h_8cm', 'theta_8cm', 'flux_8cm', 'h_9cm', 'theta_9cm', 'flux_9cm', 'h_10cm', 'theta_10cm', 'flux_10cm', 'h_15cm', 'theta_15cm', 'flux_15cm', 'h_20cm', 'theta_20cm', 'flux_20cm', 'h_30cm', 'theta_30cm', 'flux_30cm', 'h_50cm', 'theta_50cm', 'flux_50cm'), profile_1cm=FALSE) {
  subDir <- texture
  texture <- tolower(texture)
  if(climateRuns) {
    f_path <- list.dirs(file.path(directory, subDir), recursive = FALSE, full.names = FALSE)
    obs <- read.table(file.path(directory, subDir, f_path, f_path, 'Obs_Node.out'), col.names = c('time_days', 'h_0cm', 'theta_0cm', 'flux_0cm', 'h_1cm', 'theta_1cm', 'flux_1cm', 'h_2cm', 'theta_2cm', 'flux_2cm', 'h_3cm', 'theta_3cm', 'flux_3cm', 'h_4cm', 'theta_4cm', 'flux_4cm', 'h_5cm', 'theta_5cm', 'flux_5cm', 'h_6cm', 'theta_6cm', 'flux_6cm', 'h_7cm', 'theta_7cm', 'flux_7cm', 'h_8cm', 'theta_8cm', 'flux_8cm', 'h_9cm', 'theta_9cm', 'flux_9cm', 'h_10cm', 'theta_10cm', 'flux_10cm', 'h_15cm', 'theta_15cm', 'flux_15cm', 'h_20cm', 'theta_20cm', 'flux_20cm', 'h_30cm', 'theta_30cm', 'flux_30cm', 'h_50cm', 'theta_50cm', 'flux_50cm'), header = FALSE, nrow=948, skip=11)
  } else {
    obs <- read.table(file.path(directory, subDir, 'Obs_Node.out'), col.names = columnNames, header = FALSE, nrow=948, skip=11)
  }
  if (profile_1cm) {
    obs$theta_0_1cm_avg <- apply(obs[,c(3,6,9)], 1, mean)
    obs <- obs[2:nrow(obs),] #ignore time 0.001 when flooding starts
    obs$time_days <- obs$time_days - 5 #time 0 means when flooding ends
    theta_fc <- soilsDF$theta_fc_10cm[match(texture_colors$median_cokeys[texture_colors$textures==texture], soilsDF$cokey)]
    obs$theta_0_1cm_rel <- obs$theta_0_1cm_avg / theta_fc
  } else {
      obs$theta_0_10cm_avg <- apply(obs[,c(3,6,9,12,15,18,21,24,27,30,33)], 1, mean)
      obs <- obs[2:nrow(obs),] #ignore time 0.001 when flooding starts
      obs$time_days <- obs$time_days - 5 #time 0 means when flooding ends
      obs$time_diff <- c(0, obs$time_days[2:length(obs$time_days)] - obs$time_days[1:(length(obs$time_days)-1)])
      obs$cum_flux_0cm <- cumsum(obs$flux_0cm * obs$time_diff)
      obs$cum_flux_10cm <- cumsum(obs$flux_10cm * obs$time_diff)
  # obs$cum_flux_30cm <- cumsum(obs$flux_30cm * obs$time_diff)
      obs$net_flux_0_10cm <- obs$flux_0cm - obs$flux_10cm
  # obs$net_flux_0_30cm <- obs$flux_0cm - obs$flux_30cm
      theta_fc <- soilsDF$theta_fc_10cm[match(texture_colors$median_cokeys[texture_colors$textures==texture], soilsDF$cokey)]
      obs$theta_0_10cm_rel <- obs$theta_0_10cm_avg / theta_fc
  }
    obs
}
results_files <- lapply(textures, read_in_flux_files, directory=FebResultsDir)
# summary(obs$flux_0cm)
head(results_files[[1]], 10)

results_files_Apr <- lapply(textures, read_in_flux_files, directory=AprResultsDir)

sandy_loam_no_evap <- read_in_flux_files(texture = 'Sandy loam', directory = 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/HYDRUS_runs/zero_flux_tests/sandy loam', climateRuns = FALSE, columnNames = c('time_days', 'h_0cm', 'theta_0cm', 'flux_0cm', 'h_1cm', 'theta_1cm', 'flux_1cm', 'h_2cm', 'theta_2cm', 'flux_2cm', 'h_3cm', 'theta_3cm', 'flux_3cm', 'h_4cm', 'theta_4cm', 'flux_4cm', 'h_5cm', 'theta_5cm', 'flux_5cm', 'h_6cm', 'theta_6cm', 'flux_6cm', 'h_7cm', 'theta_7cm', 'flux_7cm', 'h_8cm', 'theta_8cm', 'flux_8cm', 'h_9cm', 'theta_9cm', 'flux_9cm', 'h_10cm', 'theta_10cm', 'flux_10cm', 'h_15cm', 'theta_15cm', 'flux_15cm', 'h_20cm', 'theta_20cm', 'flux_20cm', 'h_30cm', 'theta_30cm', 'flux_30cm', 'h_50cm', 'theta_50cm', 'flux_50cm', 'h_200cm', 'theta_200cm', 'flux_200cm'))
colnames(sandy_loam_no_evap)
sandy_loam_no_evap[,c('time_days', 'theta_0_10cm_avg', 'flux_10cm')]
theta_fc_est(theta_s = 0.3639, theta_r = 0.03333, Ks = 48.21478, n = 1.36604)

sandy_loam_no_evap_100cm <- read_in_flux_files(texture = 'Sandy loam', directory = 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/HYDRUS_runs/zero_flux_tests/sandy loam 100cm', climateRuns = FALSE, columnNames = c('time_days', 'h_0cm', 'theta_0cm', 'flux_0cm', 'h_1cm', 'theta_1cm', 'flux_1cm', 'h_2cm', 'theta_2cm', 'flux_2cm', 'h_3cm', 'theta_3cm', 'flux_3cm', 'h_4cm', 'theta_4cm', 'flux_4cm', 'h_5cm', 'theta_5cm', 'flux_5cm', 'h_6cm', 'theta_6cm', 'flux_6cm', 'h_7cm', 'theta_7cm', 'flux_7cm', 'h_8cm', 'theta_8cm', 'flux_8cm', 'h_9cm', 'theta_9cm', 'flux_9cm', 'h_10cm', 'theta_10cm', 'flux_10cm', 'h_15cm', 'theta_15cm', 'flux_15cm', 'h_20cm', 'theta_20cm', 'flux_20cm', 'h_30cm', 'theta_30cm', 'flux_30cm', 'h_50cm', 'theta_50cm', 'flux_50cm', 'h_100cm', 'theta_100cm', 'flux_100cm'))

sandy_loam_no_evap_10cm <- read_in_flux_files(texture = 'Sandy loam', directory = 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/HYDRUS_runs/zero_flux_tests/sandy loam 10cm', climateRuns = FALSE, columnNames = c('time_days', 'h_0cm', 'theta_0cm', 'flux_0cm', 'h_1cm', 'theta_1cm', 'flux_1cm', 'h_2cm', 'theta_2cm', 'flux_2cm', 'h_3cm', 'theta_3cm', 'flux_3cm', 'h_4cm', 'theta_4cm', 'flux_4cm', 'h_5cm', 'theta_5cm', 'flux_5cm', 'h_6cm', 'theta_6cm', 'flux_6cm', 'h_7cm', 'theta_7cm', 'flux_7cm', 'h_8cm', 'theta_8cm', 'flux_8cm', 'h_9cm', 'theta_9cm', 'flux_9cm', 'h_10cm', 'theta_10cm', 'flux_10cm'))

sandy_loam_no_evap_1cm <- read_in_flux_files(texture = 'Sandy loam', directory = 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/HYDRUS_runs/zero_flux_tests/sandy loam 1cm', climateRuns = FALSE, columnNames = c('time_days', 'h_0cm', 'theta_0cm', 'flux_0cm', 'h_0.5cm', 'theta_0.5cm', 'flux_0.5cm', 'h_1cm', 'theta_1cm', 'flux_1cm'), profile_1cm = TRUE)

#plot cumulative fluxes by textural class
for (i in seq_along(textures)) {
  plot(results_files[[i]]$time_days, results_files[[i]]$cum_flux_0cm, type = 'l', ylim=c(-2,5), xlim = c(0,30), col='black', ylab='Cumulative water flux (cm)', xlab = 'Days post-flooding')
  lines(results_files[[i]]$time_days, results_files[[i]]$cum_flux_10cm, col='grey')
  lines(results_files[[i]]$time_days, results_files[[i]]$cum_flux_30cm, col='grey', lty=2)
  abline(h=0, lty=3)
  text(x=3, y=4.5, textures[i])
}

#plot theta drawdown
for (i in seq_along(textures)) {
  plot(results_files[[i]]$time_days, results_files[[i]]$theta_0_10cm_avg, type = 'l', ylim=c(0,0.5), xlim = c(0,30), col='black', ylab=expression('0-10 cm water content (cm H'[2]*'O cm soil'^-1*')'), xlab = 'Days post-flooding')
  text(x=5, y=0.49, textures[i])
}


#clay
plot(results_files[[1]]$time_days, results_files[[1]]$theta_0_10cm_avg, type = 'l', ylim=c(0,0.5), xlim = c(0,30), col='black', ylab=expression('0-10 cm water content (cm H'[2]*'O cm soil'^-1*')'), xlab = 'Days post-flooding')
text(x=5, y=0.49, textures[1])
abline(h=0.35415, lty=3)
abline(h=0.35415*0.81, lty=3)
results_files[[1]][ ,c('time_days', 'flux_0cm', 'flux_10cm', 'net_flux_0_10cm')]
soilsDF[soilsDF$cokey==texture_colors$median_cokeys[1],]

#clay loam
plot(results_files[[4]]$time_days, results_files[[4]]$theta_0_10cm_avg, type = 'l', ylim=c(0,0.5), xlim = c(0,30), col='black', ylab=expression('0-10 cm water content (cm H'[2]*'O cm soil'^-1*')'), xlab = 'Days post-flooding')
text(x=5, y=0.49, textures[4])
abline(h=0.30604, lty=3)
abline(h=0.30604*0.82, lty=3)
results_files[[4]][ ,c('time_days', 'flux_0cm', 'flux_10cm', 'net_flux_0_10cm')]
soilsDF[soilsDF$cokey==texture_colors$median_cokeys[4],]

#loam
plot(results_files[[7]]$time_days, results_files[[7]]$theta_0_10cm_avg, type = 'l', ylim=c(0,0.5), xlim = c(0,30), col='black', ylab=expression('0-10 cm water content (cm H'[2]*'O cm soil'^-1*')'), xlab = 'Days post-flooding')
text(x=5, y=0.49, textures[7])
abline(h=0.2523672, lty=3)
abline(h=0.2523672*0.89, lty=3)
results_files[[7]][ ,c('time_days', 'flux_0cm', 'flux_10cm', 'flux_30cm', 'net_flux_0_10cm', 'theta_0_10cm_avg')]
soilsDF[soilsDF$cokey==texture_colors$median_cokeys[7],]

#sandy loam theta fc exploration
tiff(file = file.path(FiguresDir, 'FC_exploration_sandy_loam.tif'), family = 'Times New Roman', width = 6.5, height = 4.5, pointsize = 12, units = 'in', res=800, compression = 'lzw')
par(mar=c(3,3.5,0.25,0.25))
plot(results_files[[8]]$time_days, results_files[[8]]$theta_0_10cm_avg, type = 'l', ylim=c(0.1,0.36), xlim = c(0,60), col=rgb(texture_colors$red[8]/255, texture_colors$green[8]/255, texture_colors$blue[8]/255), ylab='', xlab = '', lty=1)
mtext(expression('Soil water content (cm H'[2]*'O cm soil'^-1*')'), side=2, line=2.25)
mtext('Days post-flooding', side = 1, line=2)
lines(results_files_Apr[[8]]$time_days, results_files_Apr[[8]]$theta_0_10cm_avg, type = 'l', ylim=c(0.1,0.5), xlim = c(0,30), col=rgb(texture_colors$red[8]/255, texture_colors$green[8]/255, texture_colors$blue[8]/255), lty=2)
lines(sandy_loam_no_evap_1cm$time_days, sandy_loam_no_evap_1cm$theta_1cm, col='grey', lty=3)
lines(sandy_loam_no_evap_10cm$time_days, sandy_loam_no_evap_10cm$theta_10cm, col='dark grey', lty=4)
lines(sandy_loam_no_evap$time_days, sandy_loam_no_evap$theta_1cm, col='black', lty=3)
lines(sandy_loam_no_evap$time_days, sandy_loam_no_evap$theta_10cm, col='black', lty=4)
lines(sandy_loam_no_evap$time_days, sandy_loam_no_evap$theta_30cm, col='black', lty=5)
# text(x=5, y=0.49, textures[8])
abline(h=theta_fc_est(theta_s = 0.3639, theta_r = 0.03333, Ks = 48.21478, n = 1.36604), col='red', lty=2)
points(sandy_loam_no_evap_1cm$time_days[head(which(sandy_loam_no_evap_1cm$flux_1cm > -0.01), 1)], sandy_loam_no_evap_1cm$theta_1cm[head(which(sandy_loam_no_evap_1cm$flux_1cm > -0.01), 1)], pch=1, col='grey')
points(sandy_loam_no_evap_10cm$time_days[head(which(sandy_loam_no_evap_10cm$flux_10cm > -0.01), 1)], sandy_loam_no_evap_10cm$theta_10cm[head(which(sandy_loam_no_evap_10cm$flux_10cm > -0.01), 1)], pch=1, col='grey')
points(sandy_loam_no_evap$time_days[head(which(sandy_loam_no_evap$flux_1cm > -0.01), 1)], sandy_loam_no_evap$theta_1cm[head(which(sandy_loam_no_evap$flux_1cm > -0.01), 1)], pch=1)
points(sandy_loam_no_evap$time_days[head(which(sandy_loam_no_evap$flux_10cm > -0.01), 1)], sandy_loam_no_evap$theta_10cm[head(which(sandy_loam_no_evap$flux_10cm > -0.01), 1)], pch=1)
points(sandy_loam_no_evap$time_days[head(which(sandy_loam_no_evap$flux_30cm > -0.01), 1)], sandy_loam_no_evap$theta_30cm[head(which(sandy_loam_no_evap$flux_30cm > -0.01), 1)], pch=1)
legend('topright', legend=c('30-cm depth, 200-cm profile, PET=0', '10-cm depth, 200-cm profile, PET=0', '1-cm depth, 200-cm profile, PET=0', '0-10 cm average, 200-cm profile, PET=February', '0-10 cm average, 200-cm profile, PET=April', '10-cm depth, 10-cm profile, PET=0', '1-cm depth, 1-cm profile, PET=0', 'field capacity'), lty = c(5:3, 1:2, 4:2), col = c(rep('black', 3), rep(rgb(texture_colors$red[8]/255, texture_colors$green[8]/255, texture_colors$blue[8]/255), 2), 'dark grey', 'grey', 'red'), bty = 'n')
dev.off()

#at what water content does flux at 50 cm equal 0.01 cm day-1?
results_files[[8]]$theta_50cm[head(which(results_files[[8]]$flux_50cm > -0.01), 1)]
sandy_loam_no_evap$theta_10cm[head(which(sandy_loam_no_evap$flux_10cm > -0.01), 1)] #0.2232
sandy_loam_no_evap$theta_30cm[head(which(sandy_loam_no_evap$flux_30cm > -0.01), 1)] #0.2081
sandy_loam_no_evap$theta_50cm[head(which(sandy_loam_no_evap$flux_50cm > -0.01), 1)] #NA

#at what water content does flux at 10 cm equal 0.01 cm day-1 in 10 cm profile with zero evap?
sandy_loam_no_evap_1cm$theta_1cm[head(which(sandy_loam_no_evap_10cm$flux_1cm > -0.01), 1)] #0.1866
sandy_loam_no_evap_10cm$theta_10cm[head(which(sandy_loam_no_evap_10cm$flux_10cm > -0.01), 1)] #0.1909

sandy_loam_no_evap_100cm$theta_100cm[head(which(sandy_loam_no_evap_100cm$flux_100cm > -0.01), 1)] #NA

#plot Feb vwc drawdown
tiff(file = file.path(FiguresDir, '0_10cm_absdrawdown_cell_181807_2009-02-15.tif'), family = 'Times New Roman', width = 4.5, height = 3.25, pointsize = 12, units = 'in', res=800, compression = 'lzw')
par(mar=c(3,3.5,0.25,0.25))
plot(results_files[[1]]$time_days, results_files[[1]]$theta_0_10cm_avg, type = 'l', col=rgb(texture_colors$red[1]/255, texture_colors$green[1]/255, texture_colors$blue[1]/255), ylim=c(0.08,0.48), xlim = c(0,25), ylab='', xlab = '')
mtext('Days post-flooding (February)', side = 1, line=2)
mtext(expression('0-10 cm water content (cm H'[2]*'O cm soil'^-1*')'), side=2, line=2.25, at=0.25)
lines(results_files[[2]]$time_days, results_files[[2]]$theta_0_10cm_avg, type = 'l', col=rgb(texture_colors$red[2]/255, texture_colors$green[2]/255, texture_colors$blue[2]/255))
lines(results_files[[3]]$time_days, results_files[[3]]$theta_0_10cm_avg, type = 'l', col=rgb(texture_colors$red[3]/255, texture_colors$green[3]/255, texture_colors$blue[3]/255))
lines(results_files[[4]]$time_days, results_files[[4]]$theta_0_10cm_avg, type = 'l', col=rgb(texture_colors$red[4]/255, texture_colors$green[4]/255, texture_colors$blue[4]/255))
lines(results_files[[5]]$time_days, results_files[[5]]$theta_0_10cm_avg, type = 'l', col=rgb(texture_colors$red[5]/255, texture_colors$green[5]/255, texture_colors$blue[5]/255))
lines(results_files[[6]]$time_days, results_files[[6]]$theta_0_10cm_avg, type = 'l', col=rgb(texture_colors$red[6]/255, texture_colors$green[6]/255, texture_colors$blue[6]/255))
lines(results_files[[7]]$time_days, results_files[[7]]$theta_0_10cm_avg, type = 'l', col=rgb(texture_colors$red[7]/255, texture_colors$green[7]/255, texture_colors$blue[7]/255))
lines(results_files[[8]]$time_days, results_files[[8]]$theta_0_10cm_avg, type = 'l', col=rgb(texture_colors$red[8]/255, texture_colors$green[8]/255, texture_colors$blue[8]/255))
lines(results_files[[9]]$time_days, results_files[[9]]$theta_0_10cm_avg, type = 'l', col=rgb(texture_colors$red[9]/255, texture_colors$green[9]/255, texture_colors$blue[9]/255))
lines(results_files[[10]]$time_days, results_files[[10]]$theta_0_10cm_avg, type = 'l', col=rgb(texture_colors$red[10]/255, texture_colors$green[10]/255, texture_colors$blue[10]/255))
# legend('topright', legend=texture_colors$textures, col=rgb(red=texture_colors$red/255, green = texture_colors$green/255, blue = texture_colors$blue/255), lty=1, lwd=1.5, ncol=2, cex = 0.8)
for (i in 1:nrow(texture_colors)) {
    points(x=trafficDF$fd_2009.02.15[trafficDF$cokey==texture_colors$median_cokeys[i]], y=soilsDF$theta_fc_10cm[soilsDF$cokey==texture_colors$median_cokeys[i]]*trafficDF$traff_def[trafficDF$cokey==texture_colors$median_cokeys[i]], col=rgb(texture_colors$red[i]/255, texture_colors$green[i]/255, texture_colors$blue[i]/255))
}
dev.off()

##plot Apr vwc drawdown
tiff(file = file.path(FiguresDir, '0_10cm_absdrawdown_cell_181807_2009-04-15.tif'), family = 'Times New Roman', width = 4.5, height = 3.25, pointsize = 12, units = 'in', res=800, compression = 'lzw')
par(mar=c(3,3.5,0.25,0.25))
plot(results_files_Apr[[1]]$time_days, results_files_Apr[[1]]$theta_0_10cm_avg, type = 'l', col=rgb(texture_colors$red[1]/255, texture_colors$green[1]/255, texture_colors$blue[1]/255), ylim=c(0.08,0.48), xlim = c(0,15), ylab='', xlab = '')
mtext('Days post-flooding (April)', side = 1, line=2)
mtext(expression('0-10 cm water content (cm H'[2]*'O cm soil'^-1*')'), side=2, line=2.25, at=0.25)
lines(results_files_Apr[[2]]$time_days, results_files_Apr[[2]]$theta_0_10cm_avg, type = 'l', col=rgb(texture_colors$red[2]/255, texture_colors$green[2]/255, texture_colors$blue[2]/255))
lines(results_files_Apr[[3]]$time_days, results_files_Apr[[3]]$theta_0_10cm_avg, type = 'l', col=rgb(texture_colors$red[3]/255, texture_colors$green[3]/255, texture_colors$blue[3]/255))
lines(results_files_Apr[[4]]$time_days, results_files_Apr[[4]]$theta_0_10cm_avg, type = 'l', col=rgb(texture_colors$red[4]/255, texture_colors$green[4]/255, texture_colors$blue[4]/255))
lines(results_files_Apr[[5]]$time_days, results_files_Apr[[5]]$theta_0_10cm_avg, type = 'l', col=rgb(texture_colors$red[5]/255, texture_colors$green[5]/255, texture_colors$blue[5]/255))
lines(results_files_Apr[[6]]$time_days, results_files_Apr[[6]]$theta_0_10cm_avg, type = 'l', col=rgb(texture_colors$red[6]/255, texture_colors$green[6]/255, texture_colors$blue[6]/255))
lines(results_files_Apr[[7]]$time_days, results_files_Apr[[7]]$theta_0_10cm_avg, type = 'l', col=rgb(texture_colors$red[7]/255, texture_colors$green[7]/255, texture_colors$blue[7]/255))
lines(results_files_Apr[[8]]$time_days, results_files_Apr[[8]]$theta_0_10cm_avg, type = 'l', col=rgb(texture_colors$red[8]/255, texture_colors$green[8]/255, texture_colors$blue[8]/255))
lines(results_files_Apr[[9]]$time_days, results_files_Apr[[9]]$theta_0_10cm_avg, type = 'l', col=rgb(texture_colors$red[9]/255, texture_colors$green[9]/255, texture_colors$blue[9]/255))
lines(results_files_Apr[[10]]$time_days, results_files_Apr[[10]]$theta_0_10cm_avg, type = 'l', col=rgb(texture_colors$red[10]/255, texture_colors$green[10]/255, texture_colors$blue[10]/255))
# legend('topright', legend=texture_colors$textures, col=rgb(red=texture_colors$red/255, green = texture_colors$green/255, blue = texture_colors$blue/255), lty=1, lwd=1.5, ncol=2, cex = 0.8)
for (i in 1:nrow(texture_colors)) {
  points(x=trafficDF$fd_2009.04.15[trafficDF$cokey==texture_colors$median_cokeys[i]], y=soilsDF$theta_fc_10cm[soilsDF$cokey==texture_colors$median_cokeys[i]]*trafficDF$traff_def[trafficDF$cokey==texture_colors$median_cokeys[i]], col=rgb(texture_colors$red[i]/255, texture_colors$green[i]/255, texture_colors$blue[i]/255))
}
dev.off()

#plot 0-10 cm relative water content (February)
tiff(file = file.path(FiguresDir, '0_10cm_reldrawdown_cell_181807_2009-02-15.tif'), family = 'Times New Roman', width = 4.5, height = 3.25, pointsize = 12, units = 'in', res=800)
par(mar=c(3,3.5,0.25,0.25))
plot(results_files[[1]]$time_days, results_files[[1]]$theta_0_10cm_rel, type = 'l', col=rgb(texture_colors$red[1]/255, texture_colors$green[1]/255, texture_colors$blue[1]/255), ylim=c(0.6,1.9), xlim = c(0,25), ylab='', xlab = '')
mtext('Days post-flooding (February)', side=1, line=2)
mtext('0-10 cm water content (FC fraction)', side=2, line=2.25)
lines(results_files[[2]]$time_days, results_files[[2]]$theta_0_10cm_rel, type = 'l', col=rgb(texture_colors$red[2]/255, texture_colors$green[2]/255, texture_colors$blue[2]/255))
lines(results_files[[3]]$time_days, results_files[[3]]$theta_0_10cm_rel, type = 'l', col=rgb(texture_colors$red[3]/255, texture_colors$green[3]/255, texture_colors$blue[3]/255))
lines(results_files[[4]]$time_days, results_files[[4]]$theta_0_10cm_rel, type = 'l', col=rgb(texture_colors$red[4]/255, texture_colors$green[4]/255, texture_colors$blue[4]/255))
lines(results_files[[5]]$time_days, results_files[[5]]$theta_0_10cm_rel, type = 'l', col=rgb(texture_colors$red[5]/255, texture_colors$green[5]/255, texture_colors$blue[5]/255))
lines(results_files[[6]]$time_days, results_files[[6]]$theta_0_10cm_rel, type = 'l', col=rgb(texture_colors$red[6]/255, texture_colors$green[6]/255, texture_colors$blue[6]/255))
lines(results_files[[7]]$time_days, results_files[[7]]$theta_0_10cm_rel, type = 'l', col=rgb(texture_colors$red[7]/255, texture_colors$green[7]/255, texture_colors$blue[7]/255))
lines(results_files[[8]]$time_days, results_files[[8]]$theta_0_10cm_rel, type = 'l', col=rgb(texture_colors$red[8]/255, texture_colors$green[8]/255, texture_colors$blue[8]/255))
lines(results_files[[9]]$time_days, results_files[[9]]$theta_0_10cm_rel, type = 'l', col=rgb(texture_colors$red[9]/255, texture_colors$green[9]/255, texture_colors$blue[9]/255))
lines(results_files[[10]]$time_days, results_files[[10]]$theta_0_10cm_rel, type = 'l', col=rgb(texture_colors$red[10]/255, texture_colors$green[10]/255, texture_colors$blue[10]/255))
# abline(h=1, lty=2)
abline(h=0.95, lty=2)
abline(h=0.85, lty=2)
legend('topright', legend=texture_colors$textures, col=rgb(red=texture_colors$red/255, green = texture_colors$green/255, blue = texture_colors$blue/255), lty=1, lwd=1.5, ncol=2, cex=0.8, bty = 'n')
for (i in 1:nrow(texture_colors)) {
  points(x=trafficDF$fd_2009.02.15[trafficDF$cokey==texture_colors$median_cokeys[i]], y=trafficDF$traff_def[trafficDF$cokey==texture_colors$median_cokeys[i]], col=rgb(texture_colors$red[i]/255, texture_colors$green[i]/255, texture_colors$blue[i]/255), cex=0.8)
}
dev.off()

#plot 0-10 cm relative water content (April)
tiff(file = file.path(FiguresDir, '0_10cm_reldrawdown_cell_181807_2009-04-15.tif'), family = 'Times New Roman', width = 4.5, height = 3.25, pointsize = 12, units = 'in', res=800)
par(mar=c(3,3.5,0.25,0.25))
plot(results_files_Apr[[1]]$time_days, results_files_Apr[[1]]$theta_0_10cm_rel, type = 'l', col=rgb(texture_colors$red[1]/255, texture_colors$green[1]/255, texture_colors$blue[1]/255), ylim=c(0.6,1.9), xlim = c(0,15), ylab='', xlab = '')
mtext('Days post-flooding (April)', side=1, line=2)
mtext('0-10 cm water content (FC fraction)', side=2, line=2.25)
lines(results_files_Apr[[2]]$time_days, results_files_Apr[[2]]$theta_0_10cm_rel, type = 'l', col=rgb(texture_colors$red[2]/255, texture_colors$green[2]/255, texture_colors$blue[2]/255))
lines(results_files_Apr[[3]]$time_days, results_files_Apr[[3]]$theta_0_10cm_rel, type = 'l', col=rgb(texture_colors$red[3]/255, texture_colors$green[3]/255, texture_colors$blue[3]/255))
lines(results_files_Apr[[4]]$time_days, results_files_Apr[[4]]$theta_0_10cm_rel, type = 'l', col=rgb(texture_colors$red[4]/255, texture_colors$green[4]/255, texture_colors$blue[4]/255))
lines(results_files_Apr[[5]]$time_days, results_files_Apr[[5]]$theta_0_10cm_rel, type = 'l', col=rgb(texture_colors$red[5]/255, texture_colors$green[5]/255, texture_colors$blue[5]/255))
lines(results_files_Apr[[6]]$time_days, results_files_Apr[[6]]$theta_0_10cm_rel, type = 'l', col=rgb(texture_colors$red[6]/255, texture_colors$green[6]/255, texture_colors$blue[6]/255))
lines(results_files_Apr[[7]]$time_days, results_files_Apr[[7]]$theta_0_10cm_rel, type = 'l', col=rgb(texture_colors$red[7]/255, texture_colors$green[7]/255, texture_colors$blue[7]/255))
lines(results_files_Apr[[8]]$time_days, results_files_Apr[[8]]$theta_0_10cm_rel, type = 'l', col=rgb(texture_colors$red[8]/255, texture_colors$green[8]/255, texture_colors$blue[8]/255))
lines(results_files_Apr[[9]]$time_days, results_files_Apr[[9]]$theta_0_10cm_rel, type = 'l', col=rgb(texture_colors$red[9]/255, texture_colors$green[9]/255, texture_colors$blue[9]/255))
lines(results_files_Apr[[10]]$time_days, results_files_Apr[[10]]$theta_0_10cm_rel, type = 'l', col=rgb(texture_colors$red[10]/255, texture_colors$green[10]/255, texture_colors$blue[10]/255))
# abline(h=1, lty=2)
abline(h=0.95, lty=2)
abline(h=0.85, lty=2)
# legend('topright', legend=texture_colors$textures, col=rgb(red=texture_colors$red/255, green = texture_colors$green/255, blue = texture_colors$blue/255), lty=1, lwd=1.5, ncol=2, cex=0.8, bty = 'n')
for (i in 1:nrow(texture_colors)) {
  points(x=trafficDF$fd_2009.04.15[trafficDF$cokey==texture_colors$median_cokeys[i]], y=trafficDF$traff_def[trafficDF$cokey==texture_colors$median_cokeys[i]], col=rgb(texture_colors$red[i]/255, texture_colors$green[i]/255, texture_colors$blue[i]/255), cex=0.8)
}
dev.off()

#plot 0-10 cm net daily flux
plot(results_files[[1]]$time_days, results_files[[1]]$net_flux_0_10cm, type = 'l', col=rgb(texture_colors$red[1]/255, texture_colors$green[1]/255, texture_colors$blue[1]/255), ylim=c(0,0.5), xlim = c(0,30), ylab=expression('0-10 cm net water flux (cm H'[2]*'O day'^-1*')'), xlab = 'Days post-flooding')
lines(results_files[[2]]$time_days, results_files[[2]]$net_flux_0_10cm, type = 'l', col=rgb(texture_colors$red[2]/255, texture_colors$green[2]/255, texture_colors$blue[2]/255))
lines(results_files[[3]]$time_days, results_files[[3]]$net_flux_0_10cm, type = 'l', col=rgb(texture_colors$red[3]/255, texture_colors$green[3]/255, texture_colors$blue[3]/255))
lines(results_files[[4]]$time_days, results_files[[4]]$net_flux_0_10cm, type = 'l', col=rgb(texture_colors$red[4]/255, texture_colors$green[4]/255, texture_colors$blue[4]/255))
lines(results_files[[5]]$time_days, results_files[[5]]$net_flux_0_10cm, type = 'l', col=rgb(texture_colors$red[5]/255, texture_colors$green[5]/255, texture_colors$blue[5]/255))
lines(results_files[[6]]$time_days, results_files[[6]]$net_flux_0_10cm, type = 'l', col=rgb(texture_colors$red[6]/255, texture_colors$green[6]/255, texture_colors$blue[6]/255))
lines(results_files[[7]]$time_days, results_files[[7]]$net_flux_0_10cm, type = 'l', col=rgb(texture_colors$red[7]/255, texture_colors$green[7]/255, texture_colors$blue[7]/255))
lines(results_files[[8]]$time_days, results_files[[8]]$net_flux_0_10cm, type = 'l', col=rgb(texture_colors$red[8]/255, texture_colors$green[8]/255, texture_colors$blue[8]/255))
lines(results_files[[9]]$time_days, results_files[[9]]$net_flux_0_10cm, type = 'l', col=rgb(texture_colors$red[9]/255, texture_colors$green[9]/255, texture_colors$blue[9]/255))
lines(results_files[[10]]$time_days, results_files[[10]]$net_flux_0_10cm, type = 'l', col=rgb(texture_colors$red[10]/255, texture_colors$green[10]/255, texture_colors$blue[10]/255))

#30 cm depth water content
plot(results_files[[1]]$time_days, results_files[[1]]$theta_30cm, type = 'l', col=rgb(texture_colors$red[1]/255, texture_colors$green[1]/255, texture_colors$blue[1]/255), ylim=c(0,0.5), xlim = c(0,30), ylab=expression('0-30 cm water content (cm H'[2]*'O cm soil'^-1*')'), xlab = 'Days post-flooding')
lines(results_files[[2]]$time_days, results_files[[2]]$theta_30cm, type = 'l', col=rgb(texture_colors$red[2]/255, texture_colors$green[2]/255, texture_colors$blue[2]/255))
lines(results_files[[3]]$time_days, results_files[[3]]$theta_30cm, type = 'l', col=rgb(texture_colors$red[3]/255, texture_colors$green[3]/255, texture_colors$blue[3]/255))
lines(results_files[[4]]$time_days, results_files[[4]]$theta_30cm, type = 'l', col=rgb(texture_colors$red[4]/255, texture_colors$green[4]/255, texture_colors$blue[4]/255))
lines(results_files[[5]]$time_days, results_files[[5]]$theta_30cm, type = 'l', col=rgb(texture_colors$red[5]/255, texture_colors$green[5]/255, texture_colors$blue[5]/255))
lines(results_files[[6]]$time_days, results_files[[6]]$theta_30cm, type = 'l', col=rgb(texture_colors$red[6]/255, texture_colors$green[6]/255, texture_colors$blue[6]/255))
lines(results_files[[7]]$time_days, results_files[[7]]$theta_30cm, type = 'l', col=rgb(texture_colors$red[7]/255, texture_colors$green[7]/255, texture_colors$blue[7]/255))
lines(results_files[[8]]$time_days, results_files[[8]]$theta_30cm, type = 'l', col=rgb(texture_colors$red[8]/255, texture_colors$green[8]/255, texture_colors$blue[8]/255))
lines(results_files[[9]]$time_days, results_files[[9]]$theta_30cm, type = 'l', col=rgb(texture_colors$red[9]/255, texture_colors$green[9]/255, texture_colors$blue[9]/255))
lines(results_files[[10]]$time_days, results_files[[10]]$theta_30cm, type = 'l', col=rgb(texture_colors$red[10]/255, texture_colors$green[10]/255, texture_colors$blue[10]/255))

#plot 0-30 cm net daily flux
plot(results_files[[1]]$time_days, results_files[[1]]$net_flux_0_30cm, type = 'l', col=rgb(texture_colors$red[1]/255, texture_colors$green[1]/255, texture_colors$blue[1]/255), ylim=c(0,0.5), xlim = c(0,30), ylab=expression('0-30 cm water content (cm H'[2]*'O cm soil'^-1*')'), xlab = 'Days post-flooding')
lines(results_files[[2]]$time_days, results_files[[2]]$net_flux_0_30cm, type = 'l', col=rgb(texture_colors$red[2]/255, texture_colors$green[2]/255, texture_colors$blue[2]/255))
lines(results_files[[3]]$time_days, results_files[[3]]$net_flux_0_30cm, type = 'l', col=rgb(texture_colors$red[3]/255, texture_colors$green[3]/255, texture_colors$blue[3]/255))
lines(results_files[[4]]$time_days, results_files[[4]]$net_flux_0_30cm, type = 'l', col=rgb(texture_colors$red[4]/255, texture_colors$green[4]/255, texture_colors$blue[4]/255))
lines(results_files[[5]]$time_days, results_files[[5]]$net_flux_0_30cm, type = 'l', col=rgb(texture_colors$red[5]/255, texture_colors$green[5]/255, texture_colors$blue[5]/255))
lines(results_files[[6]]$time_days, results_files[[6]]$net_flux_0_30cm, type = 'l', col=rgb(texture_colors$red[6]/255, texture_colors$green[6]/255, texture_colors$blue[6]/255))
lines(results_files[[7]]$time_days, results_files[[7]]$net_flux_0_30cm, type = 'l', col=rgb(texture_colors$red[7]/255, texture_colors$green[7]/255, texture_colors$blue[7]/255))
lines(results_files[[8]]$time_days, results_files[[8]]$net_flux_0_30cm, type = 'l', col=rgb(texture_colors$red[8]/255, texture_colors$green[8]/255, texture_colors$blue[8]/255))
lines(results_files[[9]]$time_days, results_files[[9]]$net_flux_0_30cm, type = 'l', col=rgb(texture_colors$red[9]/255, texture_colors$green[9]/255, texture_colors$blue[9]/255))
lines(results_files[[10]]$time_days, results_files[[10]]$net_flux_0_30cm, type = 'l', col=rgb(texture_colors$red[10]/255, texture_colors$green[10]/255, texture_colors$blue[10]/255))

#plot 10 cm water flux
plot(results_files[[1]]$time_days, results_files[[1]]$flux_10cm, type = 'l', col=rgb(texture_colors$red[1]/255, texture_colors$green[1]/255, texture_colors$blue[1]/255), ylim=c(-0.25,0.25), xlim = c(0,30), ylab=expression('Water flux at 10 cm depth (cm H'[2]*'O day'^-1*')'), xlab = 'Days post-flooding')
lines(results_files[[2]]$time_days, results_files[[2]]$flux_10cm, type = 'l', col=rgb(texture_colors$red[2]/255, texture_colors$green[2]/255, texture_colors$blue[2]/255))
lines(results_files[[3]]$time_days, results_files[[3]]$flux_10cm, type = 'l', col=rgb(texture_colors$red[3]/255, texture_colors$green[3]/255, texture_colors$blue[3]/255))
lines(results_files[[4]]$time_days, results_files[[4]]$flux_10cm, type = 'l', col=rgb(texture_colors$red[4]/255, texture_colors$green[4]/255, texture_colors$blue[4]/255))
lines(results_files[[5]]$time_days, results_files[[5]]$flux_10cm, type = 'l', col=rgb(texture_colors$red[5]/255, texture_colors$green[5]/255, texture_colors$blue[5]/255))
lines(results_files[[6]]$time_days, results_files[[6]]$flux_10cm, type = 'l', col=rgb(texture_colors$red[6]/255, texture_colors$green[6]/255, texture_colors$blue[6]/255))
lines(results_files[[7]]$time_days, results_files[[7]]$flux_10cm, type = 'l', col=rgb(texture_colors$red[7]/255, texture_colors$green[7]/255, texture_colors$blue[7]/255))
lines(results_files[[8]]$time_days, results_files[[8]]$flux_10cm, type = 'l', col=rgb(texture_colors$red[8]/255, texture_colors$green[8]/255, texture_colors$blue[8]/255))
lines(results_files[[9]]$time_days, results_files[[9]]$flux_10cm, type = 'l', col=rgb(texture_colors$red[9]/255, texture_colors$green[9]/255, texture_colors$blue[9]/255))
lines(results_files[[10]]$time_days, results_files[[10]]$flux_10cm, type = 'l', col=rgb(texture_colors$red[10]/255, texture_colors$green[10]/255, texture_colors$blue[10]/255))
abline(h=0, lty=2)

#plot 30 cm water flux
plot(results_files[[1]]$time_days, results_files[[1]]$flux_30cm, type = 'l', col=rgb(texture_colors$red[1]/255, texture_colors$green[1]/255, texture_colors$blue[1]/255), ylim=c(-0.25,0.12), xlim = c(0,30), ylab=expression('Water flux at 30 cm depth (cm H'[2]*'O day'^-1*')'), xlab = 'Days post-flooding')
lines(results_files[[2]]$time_days, results_files[[2]]$flux_30cm, type = 'l', col=rgb(texture_colors$red[2]/255, texture_colors$green[2]/255, texture_colors$blue[2]/255))
lines(results_files[[3]]$time_days, results_files[[3]]$flux_30cm, type = 'l', col=rgb(texture_colors$red[3]/255, texture_colors$green[3]/255, texture_colors$blue[3]/255))
lines(results_files[[4]]$time_days, results_files[[4]]$flux_30cm, type = 'l', col=rgb(texture_colors$red[4]/255, texture_colors$green[4]/255, texture_colors$blue[4]/255))
lines(results_files[[5]]$time_days, results_files[[5]]$flux_30cm, type = 'l', col=rgb(texture_colors$red[5]/255, texture_colors$green[5]/255, texture_colors$blue[5]/255))
lines(results_files[[6]]$time_days, results_files[[6]]$flux_30cm, type = 'l', col=rgb(texture_colors$red[6]/255, texture_colors$green[6]/255, texture_colors$blue[6]/255))
lines(results_files[[7]]$time_days, results_files[[7]]$flux_30cm, type = 'l', col=rgb(texture_colors$red[7]/255, texture_colors$green[7]/255, texture_colors$blue[7]/255))
lines(results_files[[8]]$time_days, results_files[[8]]$flux_30cm, type = 'l', col=rgb(texture_colors$red[8]/255, texture_colors$green[8]/255, texture_colors$blue[8]/255))
lines(results_files[[9]]$time_days, results_files[[9]]$flux_30cm, type = 'l', col=rgb(texture_colors$red[9]/255, texture_colors$green[9]/255, texture_colors$blue[9]/255))
lines(results_files[[10]]$time_days, results_files[[10]]$flux_30cm, type = 'l', col=rgb(texture_colors$red[10]/255, texture_colors$green[10]/255, texture_colors$blue[10]/255))
abline(h=0, lty=2)
abline(h=-0.01, lty=2)

#water flux at 50 cm depth
plot(results_files[[1]]$time_days, results_files[[1]]$flux_50cm, type = 'l', col=rgb(texture_colors$red[1]/255, texture_colors$green[1]/255, texture_colors$blue[1]/255), ylim=c(-0.5,0.1), xlim = c(0,30), ylab=expression('Water flux at 50 cm depth (cm H'[2]*'O day'^-1*')'), xlab = 'Days post-flooding')
lines(results_files[[2]]$time_days, results_files[[2]]$flux_50cm, type = 'l', col=rgb(texture_colors$red[2]/255, texture_colors$green[2]/255, texture_colors$blue[2]/255))
lines(results_files[[3]]$time_days, results_files[[3]]$flux_50cm, type = 'l', col=rgb(texture_colors$red[3]/255, texture_colors$green[3]/255, texture_colors$blue[3]/255))
lines(results_files[[4]]$time_days, results_files[[4]]$flux_50cm, type = 'l', col=rgb(texture_colors$red[4]/255, texture_colors$green[4]/255, texture_colors$blue[4]/255))
lines(results_files[[5]]$time_days, results_files[[5]]$flux_50cm, type = 'l', col=rgb(texture_colors$red[5]/255, texture_colors$green[5]/255, texture_colors$blue[5]/255))
lines(results_files[[6]]$time_days, results_files[[6]]$flux_50cm, type = 'l', col=rgb(texture_colors$red[6]/255, texture_colors$green[6]/255, texture_colors$blue[6]/255))
lines(results_files[[7]]$time_days, results_files[[7]]$flux_50cm, type = 'l', col=rgb(texture_colors$red[7]/255, texture_colors$green[7]/255, texture_colors$blue[7]/255))
lines(results_files[[8]]$time_days, results_files[[8]]$flux_50cm, type = 'l', col=rgb(texture_colors$red[8]/255, texture_colors$green[8]/255, texture_colors$blue[8]/255))
lines(results_files[[9]]$time_days, results_files[[9]]$flux_50cm, type = 'l', col=rgb(texture_colors$red[9]/255, texture_colors$green[9]/255, texture_colors$blue[9]/255))
lines(results_files[[10]]$time_days, results_files[[10]]$flux_50cm, type = 'l', col=rgb(texture_colors$red[10]/255, texture_colors$green[10]/255, texture_colors$blue[10]/255))
abline(h=0, lty=2)

#cumulative flux at 10 cm
plot(results_files[[1]]$time_days, results_files[[1]]$cum_flux_10cm, type = 'l', col=rgb(texture_colors$red[1]/255, texture_colors$green[1]/255, texture_colors$blue[1]/255), ylim=c(-0.5,4), xlim = c(0,30), ylab=expression('Water flux at 10 cm depth (cm H'[2]*'O day'^-1*')'), xlab = 'Days post-flooding')
lines(results_files[[2]]$time_days, results_files[[2]]$cum_flux_10cm, type = 'l', col=rgb(texture_colors$red[2]/255, texture_colors$green[2]/255, texture_colors$blue[2]/255))
lines(results_files[[3]]$time_days, results_files[[3]]$cum_flux_10cm, type = 'l', col=rgb(texture_colors$red[3]/255, texture_colors$green[3]/255, texture_colors$blue[3]/255))
lines(results_files[[4]]$time_days, results_files[[4]]$cum_flux_10cm, type = 'l', col=rgb(texture_colors$red[4]/255, texture_colors$green[4]/255, texture_colors$blue[4]/255))
lines(results_files[[5]]$time_days, results_files[[5]]$cum_flux_10cm, type = 'l', col=rgb(texture_colors$red[5]/255, texture_colors$green[5]/255, texture_colors$blue[5]/255))
lines(results_files[[6]]$time_days, results_files[[6]]$cum_flux_10cm, type = 'l', col=rgb(texture_colors$red[6]/255, texture_colors$green[6]/255, texture_colors$blue[6]/255))
lines(results_files[[7]]$time_days, results_files[[7]]$cum_flux_10cm, type = 'l', col=rgb(texture_colors$red[7]/255, texture_colors$green[7]/255, texture_colors$blue[7]/255))
lines(results_files[[8]]$time_days, results_files[[8]]$cum_flux_10cm, type = 'l', col=rgb(texture_colors$red[8]/255, texture_colors$green[8]/255, texture_colors$blue[8]/255))
lines(results_files[[9]]$time_days, results_files[[9]]$cum_flux_10cm, type = 'l', col=rgb(texture_colors$red[9]/255, texture_colors$green[9]/255, texture_colors$blue[9]/255))
lines(results_files[[10]]$time_days, results_files[[10]]$cum_flux_10cm, type = 'l', col=rgb(texture_colors$red[10]/255, texture_colors$green[10]/255, texture_colors$blue[10]/255))
abline(h=0, lty=2)
