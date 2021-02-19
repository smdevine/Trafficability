mainDir <- 'C:/Users/smdevine/Desktop/PostDoc'
trafficDir <- file.path(mainDir, 'trafficability')
cellsofinterestDir <- file.path(trafficDir, 'CIMIS')
medianDailyETo <- read.csv(file.path(cellsofinterestDir, 'SpatialCIMIS.ETo.dailymedian.csv'), stringsAsFactors = FALSE)
dim(medianDailyETo)
climsoil_df <- read.csv(file.path(trafficDir, 'climate_soil', 'climsoil_unique.csv'), stringsAsFactors = FALSE)
dim(climsoil_df)

cellnames <- paste0('cell_', climsoil_df$CIMIS)
sum(!cellnames %in% colnames(medianDailyETo)[2:ncol(medianDailyETo)]) #1 not in median daily ETo because that cell had almost no data
climsoil_df <- climsoil_df[cellnames %in% colnames(medianDailyETo)[2:ncol(medianDailyETo)],]
cellnames <- paste0('cell_', climsoil_df$CIMIS)
sum(!cellnames %in% colnames(medianDailyETo)[2:ncol(medianDailyETo)]) #0 is good

power_curv_params <- read.csv(file.path(trafficDir, 'climate_soil', 'power_curve_median.csv'), stringsAsFactors = FALSE) #this is based on cm day-1
exp_decay_params <- read.csv(file.path(trafficDir, 'climate_soil', 'exp_decay_median.csv'), stringsAsFactors = FALSE) #this is based on cm day-1

#define functions
power_func <- function(x, b, z, a) {b*x^z + a}
exp_func <- function(x, b, z, a) {b*exp(z*x) + a}
getETo_daily_mn <- function(df, cellname, start_month, start_day, days, days_to_flood=4) {
  df <- df[, c('month', 'day', cellname)]
  start_row <- which(df$month==start_month & df$day==start_day)
  days <- round(days, 0)
  mean(df[(start_row+days_to_flood):(start_row+days+days_to_flood),3])/10 #to convert from  mm to cm
}
getETo_daily_mn(df=medianDailyETo, cellname='cell_183344', start_month = 2, start_day=15, days=20)


cimis <- 182321
texture <- 'sandy loam'
initial_days <- 10
month <- 2
day <- 15
it_max <- 10

calc_days_to_traffic <- function(ET_df=medianDailyETo, texture, cimis, initial_days=10, month, day, it_max=10) {
  days_est <- initial_days
  test <- TRUE
  counter <- 0
  while(test) {
    ETo_est <- getETo_daily_mn(df=ET_df, cellname=paste0('cell_', cimis), start_month = month, start_day=day, days=days_est)
    days_est2 <- if(texture=='silty clay') {
      exp_func(ETo_est, b=exp_decay_params$b[exp_decay_params$texture_class=='silty clay'], z=exp_decay_params$z[exp_decay_params$texture_class=='silty clay'], a=exp_decay_params$a[exp_decay_params$texture_class=='silty clay'])
    } else {
        power_func(ETo_est, b=power_curv_params$b[power_curv_params$texture_class==texture], z=power_curv_params$z[power_curv_params$texture_class==texture], a=power_curv_params$a[power_curv_params$texture_class==texture])
  }
    # days_est
    # days_est2
    counter <- counter + 1
    test <- round(days_est, 2) != round(days_est2, 2)
    if(counter>it_max) {
      test <- FALSE
      days_est <- mean(c(days_est, days_est2))
    } else{days_est <- days_est2}
    # test
  }
  round(days_est, 2)
}
head(climsoil_df)
calc_days_to_traffic(texture = 'clay', cimis = 180281, month = 2, day = 15)

Jan_traffic_estimates <- mapply(FUN=calc_days_to_traffic, texture=climsoil_df$textural_class, cimis=climsoil_df$CIMIS, month=1, day=15)
Feb_traffic_estimates <- mapply(FUN=calc_days_to_traffic, texture=climsoil_df$textural_class, cimis=climsoil_df$CIMIS, month=2, day=15)
Mar_traffic_estimates <- mapply(FUN=calc_days_to_traffic, texture=climsoil_df$textural_class, cimis=climsoil_df$CIMIS, month=3, day=15)
Apr_traffic_estimates <- mapply(FUN=calc_days_to_traffic, texture=climsoil_df$textural_class, cimis=climsoil_df$CIMIS, month=4, day=15)
tapply(Jan_traffic_estimates, climsoil_df$textural_class, summary)
tapply(Feb_traffic_estimates, climsoil_df$textural_class, summary)
tapply(Mar_traffic_estimates, climsoil_df$textural_class, summary)
tapply(Apr_traffic_estimates, climsoil_df$textural_class, summary)

final_results <- cbind(climsoil_df, Jan_days=Jan_traffic_estimates, Feb_days=Feb_traffic_estimates, Mar_days=Mar_traffic_estimates, Apr_days=Apr_traffic_estimates)
dim(final_results)
lapply(final_results, class)
write.csv(final_results, file.path(trafficDir, 'climate_soil', 'Jan_Apr_traffic_est.csv'), row.names=FALSE)
#mapply(FUN=getETo_daily_mn, cellname=results$cellname, start_date=results$date, days=ceiling(results$days_to_traffic))