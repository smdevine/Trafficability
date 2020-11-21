workDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/HYDRUS_runs/Delhi/Delhi'
list.files(workDir)
t_level.out <- read.table(file.path(workDir, 'T_Level.out'), col.names = c('Time', 'rTop', 'rRoot', 'vTop', 'vRoot', 'vBot', 'sum(rTop)', 'sum(rRoot)', 'sum(vTop)', 'sum(vRoot)', 'sum(vBot)', 'hTop', 'hRoot', 'hBot', 'RunOff', 'sum(RunOff)', 'Volume', 'sum(Infil)', 'sum(Evap)', 'TLevel', 'Cum(WTrans)', 'SnowLayer'), header = FALSE, nrow=948, skip=9)
dim(t_level.out)
runoff <- max(t_level.out$sum.RunOff.)
infiltration <- max(t_level.out$sum.Infil.)

drainage <- t_level.