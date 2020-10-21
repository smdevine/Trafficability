
resultsDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/HYDRUS_runs'



#HYDRUS data processing
f_path <- 'Panoche/Panoche'
obs <- read.table(file.path(resultsDir, f_path, 'Obs_Node.out'), col.names = c('time_days', 'h_0cm', 'theta_0cm', 'flux_0cm', 'h_1cm', 'theta_1cm', 'flux_1cm', 'h_2cm', 'theta_2cm', 'flux_2cm', 'h_3cm', 'theta_3cm', 'flux_3cm', 'h_4cm', 'theta_4cm', 'flux_4cm', 'h_5cm', 'theta_5cm', 'flux_5cm', 'h_6cm', 'theta_6cm', 'flux_6cm', 'h_7cm', 'theta_7cm', 'flux_7cm', 'h_8cm', 'theta_8cm', 'flux_8cm', 'h_9cm', 'theta_9cm', 'flux_9cm', 'h_10cm', 'theta_10cm', 'flux_10cm', 'h_15cm', 'theta_15cm', 'flux_15cm', 'h_20cm', 'theta_20cm', 'flux_20cm', 'h_30cm', 'theta_30cm', 'flux_30cm', 'h_50cm', 'theta_50cm', 'flux_50cm'), header = TRUE, nrow=199, skip=10)
colnames(obs)
obs_10cm_theta <- obs[,c(1,3,6,9,12,15,18,21,24,27,30, 33)]
head(obs_10cm_theta,80)
colnames(obs_10cm_theta)

#0-10 cm depth
obs_0_10cm_theta_avg <- apply(obs_10cm_theta[2:nrow(obs_10cm_theta), 2:12], 1, mean)
obs_10cm_theta$time_days[head(which(obs_0_10cm_theta_avg < obs_10cm_theta$theta_4cm[1] * 0.9), n=1)+1] - 35 #9.57

#2-10 cm depth
obs_2_10cm_theta_avg <- apply(obs_10cm_theta[2:nrow(obs_10cm_theta), 4:12], 1, mean)
obs_10cm_theta$time_days[head(which(obs_2_10cm_theta_avg < obs_10cm_theta$theta_4cm[1] * 0.9), n=1)+1] - 35 #11.88

#6-10 cm depth
obs_6_10cm_theta_avg <- apply(obs_10cm_theta[2:nrow(obs_10cm_theta), 8:12], 1, mean)
obs_10cm_theta$time_days[head(which(obs_6_10cm_theta_avg < obs_10cm_theta$theta_4cm[1] * 0.9), n=1)+1] - 35 #16.5 days


#surface
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_0cm[2:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_4cm[1] * 0.9), n=1)+1] - 35
#7.59 days

#1 cm depth
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_1cm[2:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_4cm[1] * 0.9), n=1)+1] - 35
#8.25 days

#2 cm depth
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_2cm[2:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_4cm[1] * 0.9), n=1)+1] - 35
#8.91 days

#5 cm depth
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_5cm[2:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_4cm[1] * 0.9), n=1)+1] - 35
#11.22 days

#10 cm depth
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_10cm[2:length(obs_10cm_theta$theta_10cm)] < obs_10cm_theta$theta_4cm[1] * 0.9), n=1)+1] - 35 #21.12 days

#Feb 15-18 flooding
f_path <- 'Panoche_2/Panoche_2'
obs <- read.table(file.path(resultsDir, f_path, 'Obs_Node.out'), col.names = c('time_days', 'h_0cm', 'theta_0cm', 'flux_0cm', 'h_1cm', 'theta_1cm', 'flux_1cm', 'h_2cm', 'theta_2cm', 'flux_2cm', 'h_3cm', 'theta_3cm', 'flux_3cm', 'h_4cm', 'theta_4cm', 'flux_4cm', 'h_5cm', 'theta_5cm', 'flux_5cm', 'h_6cm', 'theta_6cm', 'flux_6cm', 'h_7cm', 'theta_7cm', 'flux_7cm', 'h_8cm', 'theta_8cm', 'flux_8cm', 'h_9cm', 'theta_9cm', 'flux_9cm', 'h_10cm', 'theta_10cm', 'flux_10cm', 'h_15cm', 'theta_15cm', 'flux_15cm', 'h_20cm', 'theta_20cm', 'flux_20cm', 'h_30cm', 'theta_30cm', 'flux_30cm', 'h_50cm', 'theta_50cm', 'flux_50cm'), header = TRUE, nrow=199, skip=10)
colnames(obs)
obs_10cm_theta <- obs[,c(1,3,6,9,12,15,18,21,24,27,30, 33)]
head(obs_10cm_theta,80)

#0-10 cm depth
obs_0_10cm_theta_avg <- apply(obs_10cm_theta[33:nrow(obs_10cm_theta), 2:12], 1, mean)
obs_10cm_theta$time_days[head(which(obs_0_10cm_theta_avg < obs_10cm_theta$theta_4cm[1] * 0.9), n=1)+32] - 49  #8.11



#2-10 cm depth
obs_2_10cm_theta_avg <- apply(obs_10cm_theta[33:nrow(obs_10cm_theta), 4:12], 1, mean)
obs_10cm_theta$time_days[head(which(obs_2_10cm_theta_avg < obs_10cm_theta$theta_4cm[1] * 0.9), n=1)+32] - 49 #10.75

#6-10 cm depth
obs_6_10cm_theta_avg <- apply(obs_10cm_theta[33:nrow(obs_10cm_theta), 8:12], 1, mean)
obs_10cm_theta$time_days[head(which(obs_6_10cm_theta_avg < obs_10cm_theta$theta_4cm[1] * 0.9), n=1)+32] - 49 #16.03 days


#surface
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_0cm[33:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_4cm[1] * 0.9), n=1)+32] - 49
#6.46 days

#1 cm depth
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_1cm[33:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_4cm[1] * 0.9), n=1)+32] - 49
#6.79 days

#2 cm depth
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_2cm[33:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_4cm[1] * 0.9), n=1)+32] - 49
#7.45 days

#5 cm depth
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_5cm[33:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_4cm[1] * 0.9), n=1)+32] - 49
#10.42 days

#10 cm depth
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_10cm[33:length(obs_10cm_theta$theta_10cm)] < obs_10cm_theta$theta_4cm[1] * 0.9), n=1)+32] - 49 #20.32 days

#Mar 1-4 flooding
f_path <- 'Panoche_3/Panoche_3'
flood_end <- as.integer(format(as.Date('2005-03-01'), format='%j')) + flood_duration - 1
obs <- read.table(file.path(resultsDir, f_path, 'Obs_Node.out'), col.names = c('time_days', 'h_0cm', 'theta_0cm', 'flux_0cm', 'h_1cm', 'theta_1cm', 'flux_1cm', 'h_2cm', 'theta_2cm', 'flux_2cm', 'h_3cm', 'theta_3cm', 'flux_3cm', 'h_4cm', 'theta_4cm', 'flux_4cm', 'h_5cm', 'theta_5cm', 'flux_5cm', 'h_6cm', 'theta_6cm', 'flux_6cm', 'h_7cm', 'theta_7cm', 'flux_7cm', 'h_8cm', 'theta_8cm', 'flux_8cm', 'h_9cm', 'theta_9cm', 'flux_9cm', 'h_10cm', 'theta_10cm', 'flux_10cm', 'h_15cm', 'theta_15cm', 'flux_15cm', 'h_20cm', 'theta_20cm', 'flux_20cm', 'h_30cm', 'theta_30cm', 'flux_30cm', 'h_50cm', 'theta_50cm', 'flux_50cm'), header = TRUE, nrow=606, skip=10)
colnames(obs)
obs_10cm_theta <- obs[,c(1,3,6,9,12,15,18,21,24,27,30, 33)]
head(obs_10cm_theta,80)
colnames(obs_10cm_theta)

#0-10 cm depth
obs_0_10cm_theta_avg <- apply(obs_10cm_theta[2:nrow(obs_10cm_theta), 2:12], 1, mean)
obs_10cm_theta$time_days[head(which(obs_0_10cm_theta_avg < obs_10cm_theta$theta_5cm[1] * 0.9), n=1)+1] - flood_end #10.1

#2-10 cm depth
obs_2_10cm_theta_avg <- apply(obs_10cm_theta[2:nrow(obs_10cm_theta), 4:12], 1, mean)
obs_10cm_theta$time_days[head(which(obs_2_10cm_theta_avg < obs_10cm_theta$theta_6cm[1] * 0.9), n=1)+1] - flood_end #11.8

#6-10 cm depth
obs_6_10cm_theta_avg <- apply(obs_10cm_theta[2:nrow(obs_10cm_theta), 8:12], 1, mean)
obs_10cm_theta$time_days[head(which(obs_6_10cm_theta_avg < obs_10cm_theta$theta_8cm[1] * 0.9), n=1)+1] - flood_end #16.4 days


#surface
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_0cm[2:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_0cm[1] * 0.9), n=1)+1] - flood_end
#8 days

#1 cm depth
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_1cm[2:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_1cm[1] * 0.9), n=1)+1] - flood_end
#8.3 days

#2 cm depth
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_2cm[2:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_2cm[1] * 0.9), n=1)+1] - flood_end
#8.8 days

#5 cm depth
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_5cm[2:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_5cm[1] * 0.9), n=1)+1] - flood_end
#11.4 days

#10 cm depth
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_10cm[2:length(obs_10cm_theta$theta_10cm)] < obs_10cm_theta$theta_10cm[1] * 0.9), n=1)+1] - flood_end #20.7 days

#Mar 15-18 flooding
f_path <- 'Panoche_4/Panoche_4'
flood_end <- as.integer(format(as.Date('2005-03-15'), format='%j')) + flood_duration - 1
obs <- read.table(file.path(resultsDir, f_path, 'Obs_Node.out'), col.names = c('time_days', 'h_0cm', 'theta_0cm', 'flux_0cm', 'h_1cm', 'theta_1cm', 'flux_1cm', 'h_2cm', 'theta_2cm', 'flux_2cm', 'h_3cm', 'theta_3cm', 'flux_3cm', 'h_4cm', 'theta_4cm', 'flux_4cm', 'h_5cm', 'theta_5cm', 'flux_5cm', 'h_6cm', 'theta_6cm', 'flux_6cm', 'h_7cm', 'theta_7cm', 'flux_7cm', 'h_8cm', 'theta_8cm', 'flux_8cm', 'h_9cm', 'theta_9cm', 'flux_9cm', 'h_10cm', 'theta_10cm', 'flux_10cm', 'h_15cm', 'theta_15cm', 'flux_15cm', 'h_20cm', 'theta_20cm', 'flux_20cm', 'h_30cm', 'theta_30cm', 'flux_30cm', 'h_50cm', 'theta_50cm', 'flux_50cm'), header = TRUE, nrow=606, skip=10)
colnames(obs)
obs_10cm_theta <- obs[,c(1,3,6,9,12,15,18,21,24,27,30, 33)]
head(obs_10cm_theta,80)
colnames(obs_10cm_theta)

#0-10 cm depth
obs_0_10cm_theta_avg <- apply(obs_10cm_theta[2:nrow(obs_10cm_theta), 2:12], 1, mean)
obs_10cm_theta$time_days[head(which(obs_0_10cm_theta_avg < obs_10cm_theta$theta_5cm[1] * 0.9), n=1)+1] - flood_end #7.4

#2-10 cm depth
obs_2_10cm_theta_avg <- apply(obs_10cm_theta[2:nrow(obs_10cm_theta), 4:12], 1, mean)
obs_10cm_theta$time_days[head(which(obs_2_10cm_theta_avg < obs_10cm_theta$theta_6cm[1] * 0.9), n=1)+1] - flood_end #10

#6-10 cm depth
obs_6_10cm_theta_avg <- apply(obs_10cm_theta[2:nrow(obs_10cm_theta), 8:12], 1, mean)
obs_10cm_theta$time_days[head(which(obs_6_10cm_theta_avg < obs_10cm_theta$theta_8cm[1] * 0.9), n=1)+1] - flood_end #15.1 days


#surface
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_0cm[2:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_0cm[1] * 0.9), n=1)+1] - flood_end
#4.6 days

#1 cm depth
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_1cm[2:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_1cm[1] * 0.9), n=1)+1] - flood_end
#5 days

#2 cm depth
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_2cm[2:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_2cm[1] * 0.9), n=1)+1] - flood_end
#5.7 days

#5 cm depth
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_5cm[2:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_5cm[1] * 0.9), n=1)+1] - flood_end
#9.4 days

#10 cm depth
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_10cm[2:length(obs_10cm_theta$theta_10cm)] < obs_10cm_theta$theta_10cm[1] * 0.9), n=1)+1] - flood_end #19.6 days

#Mar 15-18 flooding on Hanford
f_path <- 'Hanford/Hanford'
flood_end <- as.integer(format(as.Date('2005-03-15'), format='%j')) + flood_duration - 1
obs <- read.table(file.path(resultsDir, f_path, 'Obs_Node.out'), col.names = c('time_days', 'h_0cm', 'theta_0cm', 'flux_0cm', 'h_1cm', 'theta_1cm', 'flux_1cm', 'h_2cm', 'theta_2cm', 'flux_2cm', 'h_3cm', 'theta_3cm', 'flux_3cm', 'h_4cm', 'theta_4cm', 'flux_4cm', 'h_5cm', 'theta_5cm', 'flux_5cm', 'h_6cm', 'theta_6cm', 'flux_6cm', 'h_7cm', 'theta_7cm', 'flux_7cm', 'h_8cm', 'theta_8cm', 'flux_8cm', 'h_9cm', 'theta_9cm', 'flux_9cm', 'h_10cm', 'theta_10cm', 'flux_10cm', 'h_15cm', 'theta_15cm', 'flux_15cm', 'h_20cm', 'theta_20cm', 'flux_20cm', 'h_30cm', 'theta_30cm', 'flux_30cm', 'h_50cm', 'theta_50cm', 'flux_50cm'), header = TRUE, nrow=606, skip=10)
colnames(obs)
obs_10cm_theta <- obs[,c(1,3,6,9,12,15,18,21,24,27,30, 33)]
head(obs_10cm_theta,80)
colnames(obs_10cm_theta)

#0-10 cm depth
obs_0_10cm_theta_avg <- apply(obs_10cm_theta[2:nrow(obs_10cm_theta), 2:12], 1, mean)
obs_10cm_theta$time_days[head(which(obs_0_10cm_theta_avg < obs_10cm_theta$theta_5cm[1] * 0.9), n=1)+1] - flood_end #16.5

#2-10 cm depth
obs_2_10cm_theta_avg <- apply(obs_10cm_theta[2:nrow(obs_10cm_theta), 4:12], 1, mean)
obs_10cm_theta$time_days[head(which(obs_2_10cm_theta_avg < obs_10cm_theta$theta_6cm[1] * 0.9), n=1)+1] - flood_end #17.7

#6-10 cm depth
obs_6_10cm_theta_avg <- apply(obs_10cm_theta[2:nrow(obs_10cm_theta), 8:12], 1, mean)
obs_10cm_theta$time_days[head(which(obs_6_10cm_theta_avg < obs_10cm_theta$theta_8cm[1] * 0.9), n=1)+1] - flood_end #20.3 days


#surface
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_0cm[2:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_0cm[1] * 0.9), n=1)+1] - flood_end
#13.2 days

#1 cm depth
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_1cm[2:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_1cm[1] * 0.9), n=1)+1] - flood_end
#13.6 days

#2 cm depth
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_2cm[2:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_2cm[1] * 0.9), n=1)+1] - flood_end
#14.1 days

#5 cm depth
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_5cm[2:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_5cm[1] * 0.9), n=1)+1] - flood_end
#17.1 days

#10 cm depth
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_10cm[2:length(obs_10cm_theta$theta_10cm)] < obs_10cm_theta$theta_10cm[1] * 0.9), n=1)+1] - flood_end #24.4 days

#Mar 15-18 flooding on Yolo
f_path <- 'Yolo/Yolo'
flood_end <- as.integer(format(as.Date('2005-03-15'), format='%j')) + flood_duration - 1
obs <- read.table(file.path(resultsDir, f_path, 'Obs_Node.out'), col.names = c('time_days', 'h_0cm', 'theta_0cm', 'flux_0cm', 'h_1cm', 'theta_1cm', 'flux_1cm', 'h_2cm', 'theta_2cm', 'flux_2cm', 'h_3cm', 'theta_3cm', 'flux_3cm', 'h_4cm', 'theta_4cm', 'flux_4cm', 'h_5cm', 'theta_5cm', 'flux_5cm', 'h_6cm', 'theta_6cm', 'flux_6cm', 'h_7cm', 'theta_7cm', 'flux_7cm', 'h_8cm', 'theta_8cm', 'flux_8cm', 'h_9cm', 'theta_9cm', 'flux_9cm', 'h_10cm', 'theta_10cm', 'flux_10cm', 'h_15cm', 'theta_15cm', 'flux_15cm', 'h_20cm', 'theta_20cm', 'flux_20cm', 'h_30cm', 'theta_30cm', 'flux_30cm', 'h_50cm', 'theta_50cm', 'flux_50cm'), header = TRUE, nrow=606, skip=10)
colnames(obs)
obs_10cm_theta <- obs[,c(1,3,6,9,12,15,18,21,24,27,30, 33)]
head(obs_10cm_theta,80)
colnames(obs_10cm_theta)

#0-10 cm depth
obs_0_10cm_theta_avg <- apply(obs_10cm_theta[2:nrow(obs_10cm_theta), 2:12], 1, mean)
obs_10cm_theta$time_days[head(which(obs_0_10cm_theta_avg < obs_10cm_theta$theta_5cm[1] * 0.9), n=1)+1] - flood_end #47.2

#2-10 cm depth
obs_2_10cm_theta_avg <- apply(obs_10cm_theta[2:nrow(obs_10cm_theta), 4:12], 1, mean)
obs_10cm_theta$time_days[head(which(obs_2_10cm_theta_avg < obs_10cm_theta$theta_6cm[1] * 0.9), n=1)+1] - flood_end #47.8

#6-10 cm depth
obs_6_10cm_theta_avg <- apply(obs_10cm_theta[2:nrow(obs_10cm_theta), 8:12], 1, mean)
obs_10cm_theta$time_days[head(which(obs_6_10cm_theta_avg < obs_10cm_theta$theta_8cm[1] * 0.9), n=1)+1] - flood_end #49.2 days


#surface
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_0cm[2:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_0cm[1] * 0.9), n=1)+1] - flood_end
#13.2 days

#1 cm depth
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_1cm[2:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_1cm[1] * 0.9), n=1)+1] - flood_end
#44.9 days

#2 cm depth
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_2cm[2:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_2cm[1] * 0.9), n=1)+1] - flood_end
#47.5

#5 cm depth
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_5cm[2:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_5cm[1] * 0.9), n=1)+1] - flood_end
#47.5 days

#10 cm depth
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_10cm[2:length(obs_10cm_theta$theta_10cm)] < obs_10cm_theta$theta_10cm[1] * 0.9), n=1)+1] - flood_end #51.8 days

#retest Panoche
#Mar 15-18 flooding on Panoche
f_path <- 'Panoche_4/Panoche_4'
# flood_end <- as.integer(format(as.Date('2005-03-15'), format='%j')) + flood_duration - 1
flood_end <- 5
obs <- read.table(file.path(resultsDir, f_path, 'Obs_Node.out'), col.names = c('time_days', 'h_0cm', 'theta_0cm', 'flux_0cm', 'h_1cm', 'theta_1cm', 'flux_1cm', 'h_2cm', 'theta_2cm', 'flux_2cm', 'h_3cm', 'theta_3cm', 'flux_3cm', 'h_4cm', 'theta_4cm', 'flux_4cm', 'h_5cm', 'theta_5cm', 'flux_5cm', 'h_6cm', 'theta_6cm', 'flux_6cm', 'h_7cm', 'theta_7cm', 'flux_7cm', 'h_8cm', 'theta_8cm', 'flux_8cm', 'h_9cm', 'theta_9cm', 'flux_9cm', 'h_10cm', 'theta_10cm', 'flux_10cm', 'h_15cm', 'theta_15cm', 'flux_15cm', 'h_20cm', 'theta_20cm', 'flux_20cm', 'h_30cm', 'theta_30cm', 'flux_30cm', 'h_50cm', 'theta_50cm', 'flux_50cm'), header = FALSE, nrow=948, skip=11)
colnames(obs)
obs_10cm_theta <- obs[,c(1,3,6,9,12,15,18,21,24,27,30, 33)]
head(obs_10cm_theta,80)
colnames(obs_10cm_theta)

#0-10 cm depth
obs_0_10cm_theta_avg <- apply(obs_10cm_theta[2:nrow(obs_10cm_theta), 2:12], 1, mean)
obs_10cm_theta$time_days[head(which(obs_0_10cm_theta_avg < obs_10cm_theta$theta_5cm[1] * 0.9), n=1)+1] - flood_end #7.3

#2-10 cm depth
obs_2_10cm_theta_avg <- apply(obs_10cm_theta[2:nrow(obs_10cm_theta), 4:12], 1, mean)
obs_10cm_theta$time_days[head(which(obs_2_10cm_theta_avg < obs_10cm_theta$theta_6cm[1] * 0.9), n=1)+1] - flood_end #9.9

#6-10 cm depth
obs_6_10cm_theta_avg <- apply(obs_10cm_theta[2:nrow(obs_10cm_theta), 8:12], 1, mean)
obs_10cm_theta$time_days[head(which(obs_6_10cm_theta_avg < obs_10cm_theta$theta_8cm[1] * 0.9), n=1)+1] - flood_end #15.1 days


#surface
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_0cm[2:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_0cm[1] * 0.9), n=1)+1] - flood_end
#4.6 days

#1 cm depth
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_1cm[2:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_1cm[1] * 0.9), n=1)+1] - flood_end
#5 days

#2 cm depth
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_2cm[2:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_2cm[1] * 0.9), n=1)+1] - flood_end
#5.7 days

#5 cm depth
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_5cm[2:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_5cm[1] * 0.9), n=1)+1] - flood_end
#9.4 days

#10 cm depth
obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_10cm[2:length(obs_10cm_theta$theta_10cm)] < obs_10cm_theta$theta_10cm[1] * 0.9), n=1)+1] - flood_end #19.5 days

#function to summarize trafficability based on 100 day model where start time is 24 hours before
#flood_start for a number of days described by flood_duration
#trafficability assumed to be equal to 90% of field capacity
# flood_duration <- 4
# f_path <- 'Delhi/Delhi'
# FC_factor <- 0.9
determine_trafficability <- function(f_path, resultsDir, FC_factor, flood_duration) {
  # flood_end <- as.integer(format(as.Date(flood_start), format='%j')) + flood_duration - 1
  flood_end <- flood_duration + 1
  obs <- read.table(file.path(resultsDir, f_path, 'Obs_Node.out'), col.names = c('time_days', 'h_0cm', 'theta_0cm', 'flux_0cm', 'h_1cm', 'theta_1cm', 'flux_1cm', 'h_2cm', 'theta_2cm', 'flux_2cm', 'h_3cm', 'theta_3cm', 'flux_3cm', 'h_4cm', 'theta_4cm', 'flux_4cm', 'h_5cm', 'theta_5cm', 'flux_5cm', 'h_6cm', 'theta_6cm', 'flux_6cm', 'h_7cm', 'theta_7cm', 'flux_7cm', 'h_8cm', 'theta_8cm', 'flux_8cm', 'h_9cm', 'theta_9cm', 'flux_9cm', 'h_10cm', 'theta_10cm', 'flux_10cm', 'h_15cm', 'theta_15cm', 'flux_15cm', 'h_20cm', 'theta_20cm', 'flux_20cm', 'h_30cm', 'theta_30cm', 'flux_30cm', 'h_50cm', 'theta_50cm', 'flux_50cm'), header = FALSE, nrow=948, skip=11)
  obs_10cm_theta <- obs[,c(1,3,6,9,12,15,18,21,24,27,30,33)]
  obs_0_10cm_theta_avg <- apply(obs_10cm_theta[2:nrow(obs_10cm_theta), 2:12], 1, mean)

  #surface
  result_0cm <- obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_0cm[2:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_2cm[1] * FC_factor), n=1)+1] - flood_end
  print(result_0cm)
  
  #2 cm depth
  result_2cm <- obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_2cm[2:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_2cm[1] * FC_factor), n=1)+1] - flood_end
  print(result_2cm)
  
  #5 cm depth
  result_5cm <- obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_5cm[2:length(obs_10cm_theta$theta_5cm)] < obs_10cm_theta$theta_5cm[1] * FC_factor), n=1)+1] - flood_end
  print(result_5cm)
  
  #10 cm depth
  result_10cm <- obs_10cm_theta$time_days[head(which(obs_10cm_theta$theta_10cm[2:length(obs_10cm_theta$theta_10cm)] < obs_10cm_theta$theta_10cm[1] * FC_factor), n=1)+1] - flood_end #19.5 days
  print(result_10cm)
  
  #0-10 cm avg
  result_0_10cm <- obs_10cm_theta$time_days[head(which(obs_0_10cm_theta_avg < apply(obs_10cm_theta[1, 3:12], 1, mean) * FC_factor), n=1)+1] - flood_end
  print(result_0_10cm)
  results <- data.frame(days_to_traff_0cm=result_0cm, days_to_traff_2cm=result_2cm, days_to_traff_5cm=result_5cm, days_to_traff_10cm=result_10cm, days_to_traff_0_10cm=result_0_10cm)
  results
}
determine_trafficability(f_path = 'Delhi/Delhi', resultsDir = resultsDir, FC_factor=0.9, flood_duration=4)
determine_trafficability(f_path = 'Hanford/Hanford', resultsDir = resultsDir, FC_factor=0.9, flood_duration=4)
determine_trafficability(f_path = 'Yolo/Yolo', resultsDir = resultsDir, FC_factor=0.9, flood_duration=4)
determine_trafficability(f_path = 'Panoche_4/Panoche_4', resultsDir = resultsDir, FC_factor=0.9, flood_duration=4)
determine_trafficability(f_path = 'Yolo_2/Yolo_2', resultsDir = resultsDir, FC_factor=0.9, flood_duration=4)
  
