workDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/publication/climate vs. traffic math parameters'
resultsDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/publication/tables'
#function to estimate days to trafficability
power_func <- function(x, b, z, a) {b*x^z + a}
exp_func <- function(x, b, z, a) {b*exp(z*x) + a}


#mathematical constants by textural class
list.files(workDir)
power_curv_median <- read.csv(file.path(workDir, 'power_curve_median.csv'), stringsAsFactors = FALSE)
power_curv_Q1 <- read.csv(file.path(workDir, 'power_curve_q1.csv'), stringsAsFactors = FALSE)
power_curv_Q3 <- read.csv(file.path(workDir, 'power_curve_q3.csv'), stringsAsFactors = FALSE)
exp_curv_median <- read.csv(file.path(workDir, 'exp_decay_median.csv'), stringsAsFactors = FALSE)
exp_curv_Q1 <- read.csv(file.path(workDir, 'exp_decay_q1.csv'), stringsAsFactors = FALSE)
exp_curv_Q3 <- read.csv(file.path(workDir, 'exp_decay_q3.csv'), stringsAsFactors = FALSE)


days_to_traffic_est <- function(texture, exp_decay_params, power_curv_params, ETo) {
  if(texture=='none') {
    exp_func(ETo, b=exp_decay_params$b[exp_decay_params$texture_class=='silty clay'], z=exp_decay_params$z[exp_decay_params$texture_class=='silty clay'], a=exp_decay_params$a[exp_decay_params$texture_class=='silty clay'])
  } else {
    power_func(ETo, b=power_curv_params$b[power_curv_params$texture_class==texture], z=power_curv_params$z[power_curv_params$texture_class==texture], a=power_curv_params$a[power_curv_params$texture_class==texture])
  }
}

#test function
days_to_traffic_est('clay', exp_curv_median, power_curv_median, 1)
days_to_traffic_est('silty clay', exp_curv_median, power_curv_median, 6)

#median predictions
ETo_1mm_median <- sapply(power_curv_median$texture_class, days_to_traffic_est, exp_decay_params=exp_curv_median, power_curv_params=power_curv_median, ETo=1) 
ETo_1mm_median
ETo_3mm_median <- sapply(power_curv_median$texture_class, days_to_traffic_est, exp_decay_params=exp_curv_median, power_curv_params=power_curv_median, ETo=3) 
ETo_3mm_median
ETo_6mm_median <- sapply(power_curv_median$texture_class, days_to_traffic_est, exp_decay_params=exp_curv_median, power_curv_params=power_curv_median, ETo=6) 
ETo_6mm_median
ETo_median_df <- data.frame(texture=power_curv_median$texture_class, ETo_1mm=round(ETo_1mm_median, 1), ETo_3mm=round(ETo_3mm_median, 1), ETo_6mm=round(ETo_6mm_median, 1))
ETo_median_df
write.csv(ETo_median_df, file.path(resultsDir, 'days_to_traffic_median.csv'), row.names = FALSE)

#Q1 results
ETo_1mm_Q1 <- sapply(power_curv_Q1$texture_class, days_to_traffic_est, exp_decay_params=exp_curv_Q1, power_curv_params=power_curv_Q1, ETo=1) 
ETo_1mm_Q1
ETo_3mm_Q1 <- sapply(power_curv_Q1$texture_class, days_to_traffic_est, exp_decay_params=exp_curv_Q1, power_curv_params=power_curv_Q1, ETo=3) 
ETo_3mm_Q1
ETo_6mm_Q1 <- sapply(power_curv_Q1$texture_class, days_to_traffic_est, exp_decay_params=exp_curv_Q1, power_curv_params=power_curv_Q1, ETo=6) 
ETo_6mm_Q1
ETo_Q1_df <- data.frame(texture=power_curv_Q1$texture_class, ETo_1mm=round(ETo_1mm_Q1, 1), ETo_3mm=round(ETo_3mm_Q1, 1), ETo_6mm=round(ETo_6mm_Q1, 1))
ETo_Q1_df
write.csv(ETo_Q1_df, file.path(resultsDir, 'days_to_traffic_Q1.csv'), row.names = FALSE)

#Q3 results
ETo_1mm_Q3 <- sapply(power_curv_Q3$texture_class, days_to_traffic_est, exp_decay_params=exp_curv_Q3, power_curv_params=power_curv_Q3, ETo=1) 
ETo_1mm_Q3
ETo_3mm_Q3 <- sapply(power_curv_Q3$texture_class, days_to_traffic_est, exp_decay_params=exp_curv_Q3, power_curv_params=power_curv_Q3, ETo=3) 
ETo_3mm_Q3
ETo_6mm_Q3 <- sapply(power_curv_Q3$texture_class, days_to_traffic_est, exp_decay_params=exp_curv_Q3, power_curv_params=power_curv_Q3, ETo=6) 
ETo_6mm_Q3
ETo_Q3_df <- data.frame(texture=power_curv_Q3$texture_class, ETo_1mm=round(ETo_1mm_Q3, 1), ETo_3mm=round(ETo_3mm_Q3, 1), ETo_6mm=round(ETo_6mm_Q3, 1))
ETo_Q3_df
write.csv(ETo_Q3_df, file.path(resultsDir, 'days_to_traffic_Q3.csv'), row.names = FALSE)
