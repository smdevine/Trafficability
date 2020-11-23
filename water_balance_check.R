workDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/HYDRUS_runs/Delhi/Delhi'
list.files(workDir)
obs <- 948
t_level.out <- read.table(file.path(workDir, 'T_Level.out'), col.names = c('Time', 'rTop', 'rRoot', 'vTop', 'vRoot', 'vBot', 'sum(rTop)', 'sum(rRoot)', 'sum(vTop)', 'sum(vRoot)', 'sum(vBot)', 'hTop', 'hRoot', 'hBot', 'RunOff', 'sum(RunOff)', 'Volume', 'sum(Infil)', 'sum(Evap)', 'TLevel', 'Cum(WTrans)', 'SnowLayer'), header = FALSE, nrow=obs, skip=9)
dim(t_level.out)
runoff <- t_level.out$sum.RunOff.[948]
infiltration <- t_level.out$sum.Infil.[948]
infiltration_pct <- round(infiltration / (30.48 * 4) * 100, 2)
infiltration_pct
runoff_pct <- round(runoff / (infiltration+runoff), 1)
runoff_pct
water_input_error <- round(100 * (30.48 * 4 - runoff - infiltration) / (30.48 * 4), 2)
deltaS <- t_level.out$Volume[948] - t_level.out$Volume[1]
drainage <- -t_level.out$sum.vBot.[948]
evap <- t_level.out$sum.Evap.[948]
water_balance <- drainage + evap + deltaS - runoff - infiltration
water_balance_error <- round(100 * water_balance / (runoff + infiltration), 2)
results <- data.frame(water_balance_cm = water_balance, water_balance_rel_pct = water_balance_error, infiltration_cm = infiltration, runoff_cm = runoff, deltaS_cm = deltaS, drainage_cm = drainage, evaporation_cm = evap, infiltration_pct = infiltration_pct, runoff_pct = runoff_pct, water_input_error = water_input_error)
