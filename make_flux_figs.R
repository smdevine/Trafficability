resultsDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/climate_runs/typical profiles/cell_181807_2009-02-15'
subDir <- 'clay'
f_path <- list.dirs(file.path(resultsDir, subDir), recursive = FALSE, full.names = FALSE)
obs <- read.table(file.path(resultsDir, subDir, f_path, f_path, 'Obs_Node.out'), col.names = c('time_days', 'h_0cm', 'theta_0cm', 'flux_0cm', 'h_1cm', 'theta_1cm', 'flux_1cm', 'h_2cm', 'theta_2cm', 'flux_2cm', 'h_3cm', 'theta_3cm', 'flux_3cm', 'h_4cm', 'theta_4cm', 'flux_4cm', 'h_5cm', 'theta_5cm', 'flux_5cm', 'h_6cm', 'theta_6cm', 'flux_6cm', 'h_7cm', 'theta_7cm', 'flux_7cm', 'h_8cm', 'theta_8cm', 'flux_8cm', 'h_9cm', 'theta_9cm', 'flux_9cm', 'h_10cm', 'theta_10cm', 'flux_10cm', 'h_15cm', 'theta_15cm', 'flux_15cm', 'h_20cm', 'theta_20cm', 'flux_20cm', 'h_30cm', 'theta_30cm', 'flux_30cm', 'h_50cm', 'theta_50cm', 'flux_50cm'), header = FALSE, nrow=948, skip=11)
obs <- obs[2:nrow(obs),]
obs$time_days <- obs$time_days - 5
# summary(obs$flux_0cm)
plot(obs$time_days, obs$flux_0cm, type = 'l', ylim=c(-0.35, 0.35), xlim = c(0,30), col='black', ylab=expression('Water flux (cm day'^-1*')'), xlab = 'Days post-flooding')
lines(obs$time_days, obs$flux_10cm, col='grey')
lines(obs$time_days, obs$flux_30cm, col='grey', lty=2)
abline(h=0, lty=3)

