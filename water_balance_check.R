# workDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/HYDRUS_runs/Delhi/Delhi'
# list.files(workDir)
laptop <- FALSE
if (laptop) {
  resultsDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/HYDRUS_runs'
  summaryDir <- file.path(resultsDir, 'summary')
  summaryDir2 <- file.path(resultsDir, 'summary')} else {
    resultsDir <- 'D:/PostDoc/Trafficability/Oct2020test'
    summaryDir <- file.path(resultsDir, 'summary')
    summaryDir2 <- file.path(resultsDir, 'summary_v2')
  }

check_water_blnc <- function(f_path, obs=948, workDir=resultsDir, floodDuration=4) {
  t_level.out <- read.table(file.path(workDir, f_path, f_path, 'T_Level.out'), col.names = c('Time', 'rTop', 'rRoot', 'vTop', 'vRoot', 'vBot', 'sum(rTop)', 'sum(rRoot)', 'sum(vTop)', 'sum(vRoot)', 'sum(vBot)', 'hTop', 'hRoot', 'hBot', 'RunOff', 'sum(RunOff)', 'Volume', 'sum(Infil)', 'sum(Evap)', 'TLevel', 'Cum(WTrans)', 'SnowLayer'), header = FALSE, nrow=obs, skip=9)
  runoff <- t_level.out$sum.RunOff.[948]
  infiltration <- t_level.out$sum.Infil.[948]
  # infiltration_pct <- round(100*(infiltration / (infiltration+runoff)), 2)
  # runoff_pct <- round(100 * (runoff / (infiltration+runoff)), 2)
  water_input_error <- round(100 * ((30.48 * floodDuration - runoff - infiltration) / (30.48 * floodDuration)), 2)
  thetaS_ck <- t_level.out$Volume[2]
  thetaS_max <- max(t_level.out$Volume, na.rm = TRUE)
  deltaS <- t_level.out$Volume[948] - t_level.out$Volume[1]
  drainage <- -t_level.out$sum.vBot.[948]
  evap <- t_level.out$sum.Evap.[948]
  water_balance <- drainage + evap + deltaS - infiltration
  water_balance_error <- round(100 * water_balance / infiltration, 2)
  deltaS_PF <- t_level.out$Volume[948] - t_level.out$Volume[2]
  drainage_PF <- -(t_level.out$sum.vBot.[948] - t_level.out$sum.vBot.[2])
  evap_PF <- t_level.out$sum.Evap.[948] - t_level.out$sum.Evap.[2]
  infiltration_PF <- t_level.out$sum.Infil.[948] - t_level.out$sum.Infil.[2]
  water_balance <- drainage + evap + deltaS - infiltration
  water_balance_error <- round(100 * (water_balance / infiltration), 2)
  water_balance_PF <- drainage_PF + evap_PF + deltaS_PF - infiltration_PF
  results <- data.frame(water_balance_cm = water_balance, water_balance_rel_pct = water_balance_error, water_balance_PF_cm=water_balance_PF, water_input_error = water_input_error, infiltration_cm = infiltration, runoff_cm = runoff, deltaS_cm = deltaS, drainage_cm = drainage, evaporation_cm = evap, deltaS_PF_cm = deltaS_PF, drainage_PF_cm = drainage_PF, evap_PF_cm = evap_PF, runoff_PF_cm = (t_level.out$sum.RunOff.[948] - t_level.out$sum.RunOff.[2]), infiltration_PF_cm = infiltration_PF, thetaS_ck = thetaS_ck, thetaS_max = thetaS_max) 
  results
}

path_names <- list.dirs(resultsDir, full.names = FALSE, recursive = FALSE)
head(path_names)
tail(path_names)
path_names <- path_names[grepl('soil_', path_names)]
length(path_names)
error_runs <- sapply(path_names, function(x) file.exists(file.path(resultsDir, x, x, 'Error.msg')))
sum(error_runs) #15 of 5509
error_names <- path_names[error_runs]
path_names <- path_names[!error_runs]
length(path_names)
path_names[1]

check_water_blnc(f_path = path_names[6])

water_blnc_results <- do.call(rbind, lapply(path_names, check_water_blnc))
head(water_blnc_results)
water_blnc_results$cokey <- gsub('soil_', '', path_names)
water_blnc_results$thetaS_test <- abs(water_blnc_results$thetaS_ck - water_blnc_results$thetaS_max)
sum(water_blnc_results$thetaS_test!=0, na.rm = TRUE) #7
write.csv(water_blnc_results, file.path(resultsDir, 'water_blnc_check', 'water_blnc_results.csv'), row.names = FALSE)
lapply(water_blnc_results, summary)
sum(abs(water_blnc_results$water_balance_rel_pct) > 5, na.rm = TRUE) #33 have > 5% rel. error; 133 are NA
sum(abs(water_blnc_results$water_balance_PF_cm > 0.5), na.rm = TRUE) #only 6 have > 0.5 cm post-flood water balance error
sum(water_blnc_results$infiltration_PF_cm > 0) #7 have infiltration continuing post-flood
sum(water_blnc_results$runoff_PF_cm > 0)
