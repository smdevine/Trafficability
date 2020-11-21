library(aqp)
options(max.print = 10000)
options(width=120)

workDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/data from Stathis'
templateDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/Nov2020test/Template/Template'
modelDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/Nov2020test'
climateDir <- 'C:/Users/smdevine/Desktop/Allowable_Depletion/model_scaffold/run_model/Mar2018'
summaryDir <- file.path(modelDir, 'summary')
list.files(templateDir)

#set flood date and duration
# flood_date <- '2005-03-15'
# flood_duration <- 4

#create PET vector for testing general soil differences in trafficability
ETo <- read.csv(file.path(climateDir, 'SpatialCIMIS.ETo.QCpass.csv'), stringsAsFactors = FALSE)
tapply(ETo$cell_170036, ETo$year, sum)
Panoche_ETo_170036 <- ETo[, c('dates', 'cell_170036')]
Panoche_ETo_2005 <- Panoche_ETo_170036[which(Panoche_ETo_170036$dates=='01_01_2005'):which(Panoche_ETo_170036$dates=='12_31_2005'),2] / 10 #to convert from cm to mm

#modify soil profile according to VG params
#test run with Hanford entry from modelling_database.R
mod_database <- read.csv(file.path(workDir, 'modelling_database.csv'), stringsAsFactors = FALSE)
mod_database <- mod_database[which(mod_database$hzn_top<200),] #delete horizons that start at or below 200 cm depth
# unique(mod_database$cokey[mod_database$taxonname=='Delhi'])
# unique(mod_database$cokey[mod_database$taxonname=='Yolo'])
# unique(mod_database$cokey[mod_database$taxonname=='Hanford'])
# unique(mod_database$cokey[mod_database$taxonname=='PANOCHE'])

# unique(mod_database$hzn_desgn[mod_database$Roseta.model==1])

#get data for soil name; only will work if soil has exactly one  instance
# select_soil <- function(soilname, cokey) {
#   soil <- mod_database[mod_database$taxonname==soilname & mod_database$cokey==cokey,]
#   mat_number <- nrow(soil)
#   VGs <- soil[,c(20:24)]
#   VGs$l <- 0.5
#   colnames(VGs) <- NULL
#   depths <- soil[,8:9]
#   depths$hzn_bot[length(depths$hzn_bot)] <- 201
#   list(soil=soil, VGs=VGs, depths=depths, mat_number=mat_number)
# }
# hanford <- select_soil('Hanford', 2793)
# delhi <- select_soil('Delhi', 7276)
# yolo <- select_soil('Yolo', 10997) #11067
# yolo_2 <- select_soil('Yolo', 11067)
# panoche <- select_soil('PANOCHE', 6829)

#get data for each cokey
soil_by_cokey <- function(cokey) {
  soil <- mod_database[mod_database$cokey==cokey,]
  mat_number <- nrow(soil)
  soil <- soil[order(soil$hzn_top),]
  VGs <- soil[,c(20:24)]
  VGs$l <- 0.5
  colnames(VGs) <- NULL
  depths <- soil[,8:9]
  depths$hzn_bot[length(depths$hzn_bot)] <- 201
  list(soil=soil, VGs=VGs, depths=depths, mat_number=mat_number)
}
unique_cokeys <- unique(mod_database$cokey)
names(unique_cokeys) <- paste0('soil_', unique_cokeys)
cemented_soils <- sapply(unique_cokeys, function(x) {ifelse(sum(mod_database$Roseta.model[mod_database$cokey==x]==1)>0, TRUE, FALSE)}, USE.NAMES = TRUE)

mod_database_aqp <- mod_database
depths(mod_database_aqp) <- cokey ~ hzn_top + hzn_bot
depth_logic <- checkHzDepthLogic(mod_database_aqp)
head(depth_logic)
sum(!depth_logic$valid) #23 are not valid
sum(depth_logic$missingDepth) #0
depth_logic_T <- depth_logic[!depth_logic$valid,]
depth_logic_T
for (i in 1:nrow(depth_logic_T)) {
  print(mod_database[mod_database$cokey==depth_logic_T$cokey[i], c('hzn_top', 'hzn_bot')])
}  #all have lower horizon top equal to bottom 

head(cemented_soils)
sum(cemented_soils) #173
soil_data <- lapply(unique_cokeys, function(x) {soil_by_cokey(x)})
length(soil_data) #5685
soil_data <- soil_data[!cemented_soils]
length(soil_data) #5512, now has cemented soils removed

#find soils with missing surface horizon
soils_w_topsoil <- sapply(unique_cokeys, function(x) {0 %in% mod_database$hzn_top[mod_database$cokey==x]}, USE.NAMES = TRUE)
sum(!soils_w_topsoil) #4 missing topsoil
names(soils_w_topsoil[!soils_w_topsoil])

soil_data <- soil_data[!(names(soil_data) %in% names(soils_w_topsoil[!soils_w_topsoil]))]
length(soil_data) #5509, now has soils with no topsoil removed

#order in selector.in file: thr     ths    Alfa      n         Ks
  
# Hanford <- mod_database[mod_database$taxonname=='Hanford',]
# nrow(Hanford)
# mat_number <- nrow(Hanford)
# colnames(Hanford)
# Hanford_VGs <- Hanford[,c(20:24)]
# colnames(Hanford_VGs)
# Hanford_VGs$l <- 0.5
# colnames(Hanford_VGs) <- NULL
# Hanford_depths <- Hanford[,8:9]
# Hanford_depths$hzn_bot[length(Hanford_depths$hzn_bot)] <- 201

30.48*4 #121.9
#modify ATMOSPH.IN file according to flood date and duration; will eventually need to modify PET vector
flood_inputs <- function(vector_input, days_to_flood, application=30.48) {
  # vector_input[as.integer(format(as.Date(start_date), format='%j')):(as.integer(format(as.Date(start_date), format='%j')) + days_to_flood - 1)] <- application #old version running annual model
  vector_input[2:(days_to_flood+1)] <- application
  vector_input
}

write_atmos <- function(modelDir, flood_date, flood_duration, soil_name, PET, climate_rows=100) {
  if(!dir.exists(file.path(modelDir, soil_name))) {
    dir.create(file.path(modelDir, soil_name))
    dir.create(file.path(modelDir, soil_name, soil_name))
  }
  atmos_header <- readLines(file.path(templateDir, 'ATMOSPH.IN'), n=9)
  atmos_data <- read.table(file.path(templateDir, 'ATMOSPH.IN'), header = FALSE, nrow=climate_rows, skip=9)
  atmos_data$V2 <- 0
  atmos_data$V2 <- flood_inputs(atmos_data$V2, days_to_flood=flood_duration)
  atmos_data$V3 <- PET[(as.integer(format(as.Date(flood_date), format='%j'))-1):(as.integer(format(as.Date(flood_date), format='%j'))+climate_rows-2)]
  atmos_colnames <- unlist(strsplit(atmos_header[9], ' '))
  atmos_colnames <- atmos_colnames[sapply(atmos_colnames, function(x) x!='')]
  atmos_data$V9 <- ""
  colnames(atmos_data) <- atmos_colnames
  atmos_header <- atmos_header[-9]
  atmos_output <- file(file.path(modelDir, soil_name, soil_name, 'ATMOSPH.IN'), 'w')
  writeLines(atmos_header, con = atmos_output)
  # write.table(atmos_data, file = atmos_output, append=TRUE, col.names = FALSE, row.names = FALSE, sep='    ')
  capture.output(print.data.frame(atmos_data, row.names = FALSE, print.gap=7), file = atmos_output, append = TRUE)
  writeLines("end*** END OF INPUT FILE 'ATMOSPH.IN' **********************************", con = atmos_output)
  close(atmos_output)
}
# cbind(Panoche_ETo_2005)
# write_atmos(modelDir = modelDir, flood_date = "2005-03-15", flood_duration = 4, soil_name='Delhi', PET = Panoche_ETo_2005)
# write_atmos(modelDir = modelDir, flood_date = "2005-03-15", flood_duration = 4, soil_name='Yolo', PET = Panoche_ETo_2005)
# write_atmos(modelDir = modelDir, flood_date = "2005-03-15", flood_duration = 4, soil_name='Yolo_2', PET = Panoche_ETo_2005)
# write_atmos(modelDir = modelDir, flood_date = "2005-03-15", flood_duration = 4, soil_name='Hanford', PET = Panoche_ETo_2005)
# write_atmos(modelDir = modelDir, flood_date = "2005-03-15", flood_duration = 4, soil_name='Panoche_4', PET = Panoche_ETo_2005)  
# write_atmos(modelDir = modelDir, flood_date = "2005-03-15", flood_duration = 4, soil_name=names(soil_data)[1], PET = Panoche_ETo_2005)


#SELECTOR.IN file creation
#print_start was 4
write_selector <- function(mat_number, flood_duration, modelDir, VGs, soil_name, MaxIT) { #flood_date,
  # if(!dir.exists(file.path(modelDir, soil_name))) {
  #   dir.create(file.path(modelDir, soil_name))
  #   dir.create(file.path(modelDir, soil_name, soil_name))
  # }
  modelDir <- file.path(modelDir, soil_name, soil_name)
  selector_header <- readLines(file.path(templateDir, 'SELECTOR.IN'), n=39)
  selector_header2 <- selector_header[32:39]
  selector_header <- selector_header[1:26]
  # length(unlist(strsplit(selector_header2[3], ''))) #to figure out where total number of print times is located in file with 5 material layers
  # substring(selector_header[14], 3, 3) #number of material layers
  substring(selector_header[14], ifelse(mat_number<10, 3, 2), 3) <- as.character(mat_number)
  substring(selector_header[17], 3, 4) <- as.character(MaxIT)
  # substring(selector_header2[3], 67, 69) #199 was original print times
  # print_start <- as.integer(format(as.Date(flood_date), format='%j')) + flood_duration - 1
  print_start <- 1 + flood_duration
  # print_times <- c(0.001, seq(from = print_start, to=99.5, by=0.1), 100)
  print_times <- c(0.001, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.2, 4.4, 4.6, 4.8, seq(from = print_start, to=99.5, by=0.1), 100)
  substring(selector_header2[3], 67, 69) <- as.character(length(print_times)) #inserts correct number of print times
  print_times_df <- as.data.frame(matrix(data=print_times, nrow=160, ncol=6, byrow = TRUE)) #nrow=158 if drop flood obs
  colnames(print_times_df) <- NULL
#write result to file
  selector_output <- file(file.path(modelDir, 'SELECTOR.IN'), 'w')
  writeLines(selector_header, con = selector_output)
  write.table(format(VGs, justify = 'right'), file = selector_output, row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)
# capture.output(print.data.frame(Hanford_VGs, row.names = FALSE, print.gap=2), file = selector_output, append = TRUE)
  writeLines(selector_header2, con = selector_output)
# write.table(format(print_times_df, justify = 'right'), file = selector_output, row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)
  capture.output(print.data.frame(print_times_df, row.names = FALSE, print.gap=6), file = selector_output, append = TRUE)
  writeLines("*** END OF INPUT FILE 'SELECTOR.IN' ************************************", con = selector_output)
  close(selector_output)
}
#test function
# write_selector(mat_number = delhi$mat_number, flood_duration = 4, modelDir = modelDir, VGs = delhi$VGs, soil_name = 'Delhi', MaxIT = 30)
# write_selector(mat_number = yolo$mat_number, flood_duration = 4, modelDir = modelDir, VGs = yolo$VGs, soil_name = 'Yolo', MaxIT = 30)
# write_selector(mat_number = yolo_2$mat_number, flood_duration = 4, modelDir = modelDir, VGs = yolo_2$VGs, soil_name = 'Yolo_2', MaxIT = 30)
# write_selector(mat_number = panoche$mat_number, flood_duration = 4, modelDir = modelDir, VGs = panoche$VGs, soil_name = 'Panoche_4', MaxIT = 30)
# write_selector(mat_number = hanford$mat_number, flood_duration = 4, modelDir = modelDir, VGs = hanford$VGs, soil_name = 'Hanford', MaxIT = 30)
# write_selector(mat_number = soil_data[[1]]$mat_number, flood_duration = 4, modelDir = modelDir, VGs = soil_data[[1]]$VGs, soil_name = names(soil_data)[1], MaxIT = 30)

#now modify profile.dat
format_profile <-  function(x) {
  sub('([0-9])$', '0\\1', sprintf("%.6e", x))
}
theta_fc_est <- function(theta_s, theta_r, Ks, n) {
  (n^(-0.6*(2+log10(Ks))))*(theta_s - theta_r) + theta_r
}
vg_theta <- function(theta_r, theta_s, alpha, n, h) {
  theta_r + (theta_s - theta_r) / (1 + (alpha * h)^n)^(1-1/n)
}
find_h_at_theta_fc <- function(h, theta_r, theta_s, alpha, n, theta_fc){
  abs(theta_fc  - theta_r - ((theta_s - theta_r) / (1 + (alpha * h)^n)^(1-1/n)))
}
# vg_theta(theta_r = 0.05872, theta_s = 0.37117, alpha = 0.0335, n=1.34754, h=160)
# theta_fc_est(theta_s = 0.37117, theta_r = 0.05872, Ks=25.381, n=1.34754) #0.2286089 correct according to HYDRUS
# optimize(find_h_at_theta_fc, interval = c(10,10000), theta_r = 0.05872, theta_s = 0.37117, alpha = 0.0335, n=1.34754, theta_fc=0.2286089) #160.129
write_profile <- function(depths, soil, mat_number, modelDir, soil_name) {
  profile.dat <- readLines(file.path(templateDir, 'PROFILE.DAT'))
  profile_header <- profile.dat[1:5]
  profile_tail <- tail(profile.dat, 2)
  z <- c(seq(0,14.5,0.5), c(seq(15,200, 1)))
  z_formatted <- format_profile(-z)
# z_formatted <- formatC(-z, digits=6, format='e')
# mat_choices <- 1:mat_number
  mat_numbers <- sapply(z, function(x) which(x >= depths$hzn_top & x < depths$hzn_bot))
# cbind(z, mat_numbers)
  soil$theta_fc <- theta_fc_est(soil$teta_s, soil$teta_r, soil$Ks..cm.d., soil$n)
  soil$h_fc <- sapply(1:mat_number, function(i) {optimize(find_h_at_theta_fc, interval = c(10,10000), theta_r = soil$teta_r[i], theta_s = soil$teta_s[i], alpha = soil$alpha..1.cm.[i], n=soil$n[i], theta_fc=soil$theta_fc[i])$minimum})
  time_zero_h <- -soil$h_fc[sapply(z, function(x) which(x >= depths$hzn_top & x < depths$hzn_bot))]
# time_zero_h_formatted <- formatC(time_zero_h, digits=6, format='e')'
  time_zero_h_formatted <- format_profile(time_zero_h)
# sub('([1-9])$', '0\\1', sprintf("%.6e", time_zero_h))
  profile_data <- as.data.frame(cbind(1:length(z), z_formatted,  time_zero_h_formatted, mat_numbers, 1, format_profile(0), format_profile(1), format_profile(1), format_profile(1)))
#write PROFILE.DAT to file
  profile_output <- file(file.path(modelDir, soil_name, soil_name, 'PROFILE.DAT'), 'w')
  writeLines(profile_header, con = profile_output)
  write.table(format(profile_data, trim=TRUE, justify = 'right', width = 0), file = profile_output, row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE, sep = ' ')
  writeLines(profile_tail, con = profile_output)
  close(profile_output)
}
# write_profile(depths=delhi$depths, soil = delhi$soil, mat_number = delhi$mat_number, modelDir = modelDir, soil_name = 'Delhi')
# write_profile(depths=yolo$depths, soil = yolo$soil, mat_number = yolo$mat_number, modelDir = modelDir, soil_name = 'Yolo')
# write_profile(depths=yolo_2$depths, soil = yolo_2$soil, mat_number = yolo_2$mat_number, modelDir = modelDir, soil_name = 'Yolo_2')
# write_profile(depths=hanford$depths, soil = hanford$soil, mat_number = hanford$mat_number, modelDir = modelDir, soil_name = 'Hanford')
# write_profile(depths=panoche$depths, soil = panoche$soil, mat_number = panoche$mat_number, modelDir = modelDir, soil_name = 'Panoche_4')
# write_profile(depths=soil_data[[1]]$depths, soil = soil_data[[1]]$soil, mat_number = soil_data[[1]]$mat_number, modelDir = modelDir, soil_name = names(soil_data)[1])

#create directory of input data for all soil names that aren't cemented
#faulty data: soil_2676
for(i in 1:length(soil_data)) {
  write_atmos(modelDir = modelDir, flood_date = "2005-03-15", flood_duration = 4, soil_name=names(soil_data)[i], PET = Panoche_ETo_2005) #PET will be modified in future runs
  write_selector(mat_number = soil_data[[i]]$mat_number, flood_duration = 4, modelDir = modelDir, VGs = soil_data[[i]]$VGs, soil_name = names(soil_data)[i], MaxIT = 30)
  write_profile(depths=soil_data[[i]]$depths, soil = soil_data[[i]]$soil, mat_number = soil_data[[i]]$mat_number, modelDir = modelDir, soil_name = names(soil_data)[i])
}

#write file paths and run.bat file
#5509 soils to model
#copy path1 level_01.dir
#H1D_CALC
paths_shortcut <- paste0('path', 1:5509)
runbat <- file(file.path(modelDir, 'run.bat.txt'), 'w')
for (i in seq_along(paths_shortcut)) {
  writeLines(paste('copy', paths_shortcut[i], 'level_01.dir'), con = runbat)
  writeLines('H1D_CALC', con = runbat)
}
close(runbat)

for (i in seq_along(paths_shortcut)) {
  filepth <- file(file.path(modelDir, 'paths', paths_shortcut[i]), 'w')
  writeLines(file.path('D:\PostDoc\Trafficability\Oct2020test', names(soil_data)[i], names(soil_data)[i], fsep='\\'), con = filepth)
  close(filepth)
}

#look at earlier model run results (before flooding ends)
cokeys_to_investigate <-c('12532', '7183', '15553', '11874', '7322', '16891', '7185', '20199', '20887', '8040', '2747', '6061', '11067', '6035', '2917', '69919', '2750', '21920')
cokeys_to_investigate <- paste0('soil_', cokeys_to_investigate)
indices_to_investigate <- sapply(cokeys_to_investigate, function(x) which(names(soil_data)==x), USE.NAMES = FALSE)
paths_shortcut <- paste0('path', 1:length(cokeys_to_investigate))
runbat <- file(file.path(modelDir, 'run.bat.txt'), 'w')
for (i in seq_along(paths_shortcut)) {
  writeLines(paste('copy', paths_shortcut[i], 'level_01.dir'), con = runbat)
  writeLines('H1D_CALC', con = runbat)
}
close(runbat)
for (i in seq_along(paths_shortcut)) {
  filepth <- file(file.path(modelDir, 'paths', paths_shortcut[i]), 'w')
  writeLines(file.path(modelDir, names(soil_data)[indices_to_investigate[i]], names(soil_data)[indices_to_investigate[i]], fsep='\\'), con = filepth)
  close(filepth)
}
for(i in indices_to_investigate) {
  write_atmos(modelDir = modelDir, flood_date = "2005-03-15", flood_duration = 4, soil_name=names(soil_data)[i], PET = Panoche_ETo_2005) #PET will be modified in future runs
  write_selector(mat_number = soil_data[[i]]$mat_number, flood_duration = 4, modelDir = modelDir, VGs = soil_data[[i]]$VGs, soil_name = names(soil_data)[i], MaxIT = 30)
  write_profile(depths=soil_data[[i]]$depths, soil = soil_data[[i]]$soil, mat_number = soil_data[[i]]$mat_number, modelDir = modelDir, soil_name = names(soil_data)[i])
}

#determine days to trafficability with revised observation node print out
determine_trafficability_v3 <- function(f_path, resultsDir, flood_duration) {
  print(f_path)
  # flood_end <- as.integer(format(as.Date(flood_start), format='%j')) + flood_duration - 1
  flood_end <- flood_duration + 1
  if(!file.exists(file.path(resultsDir, f_path, f_path, 'Obs_Node.out'))) {
    print(paste('No results for', f_path))
    next
  }
  # file.remove(file.path(resultsDir, f_path, f_path, 'Nod_Inf.out'), showWarnings=FALSE)
  # file.remove(file.path(resultsDir, f_path, f_path, 'Nod_Inf_V.out'), showWarnings=FALSE)
  obs <- read.table(file.path(resultsDir, f_path, f_path, 'Obs_Node.out'), col.names = c('time_days', 'h_0cm', 'theta_0cm', 'flux_0cm', 'h_1cm', 'theta_1cm', 'flux_1cm', 'h_2cm', 'theta_2cm', 'flux_2cm', 'h_3cm', 'theta_3cm', 'flux_3cm', 'h_4cm', 'theta_4cm', 'flux_4cm', 'h_5cm', 'theta_5cm', 'flux_5cm', 'h_6cm', 'theta_6cm', 'flux_6cm', 'h_7cm', 'theta_7cm', 'flux_7cm', 'h_8cm', 'theta_8cm', 'flux_8cm', 'h_9cm', 'theta_9cm', 'flux_9cm', 'h_10cm', 'theta_10cm', 'flux_10cm', 'h_15cm', 'theta_15cm', 'flux_15cm', 'h_20cm', 'theta_20cm', 'flux_20cm', 'h_30cm', 'theta_30cm', 'flux_30cm', 'h_50cm', 'theta_50cm', 'flux_50cm'), header = FALSE, nrow=948, skip=11)
  obs_10cm_theta <- obs[,c(1,3,6,9,12,15,18,21,24,27,30,33)]
  obs_0_10cm_theta_avg <- apply(obs_10cm_theta[14:nrow(obs_10cm_theta), 2:12], 1, mean) #row 14 is beginnning of day 5 now
  
  #0-10 cm avg
  days_to_traff_0_10cm <- sapply(seq(0.95, 0.7, -0.01), function(x) {
    obs_10cm_theta$time_days[head(which(obs_0_10cm_theta_avg < apply(obs_10cm_theta[1, 3:12], 1, mean) * x), n=1)+13] - flood_end}) #FC estimate ignores upper layer; add 13 to match correct row
  
  results <- matrix(data=days_to_traff_0_10cm, nrow = length(days_to_traff_0_10cm), ncol = 1, byrow = TRUE, dimnames = list(paste0('FC_', seq(0.95, 0.7, -0.01)), c('days_to_traff_0_10cm')))
  write.csv(results, file = file.path(summaryDir, paste0(f_path, '_results.csv')), row.names=TRUE)
}
sapply(cokeys_to_investigate, determine_trafficability_v3, resultsDir=modelDir, flood_duration=4)
fc_thresholds <- c(0.81, 0.81, 0.82, 0.82, 0.89, 0.89, 0.89, 0.91, 0.85, 0.9, 0.9, 0.9, 0.89, 0.89, 0.8, 0.8, 0.82, 0.82)
fnames_results <- list.files(summaryDir, full.names = FALSE, recursive = FALSE)
extract_0_10cm_trafficability_v2 <- function(x) {
  results <- read.csv(file.path(summaryDir, x), row.names=1, na.strings = 'numeric(0)')
}
trafficability_trial <- lapply(fnames_results, extract_0_10cm_trafficability_v2)
cokeys_modeled <- sub(pattern = 'soil_', replacement = '', x = fnames_results)
cokeys_modeled <- sub(pattern = '_results.csv', replacement = '', x = cokeys_modeled)
names(trafficability_trial) <- cokeys_modeled
df_results <- data.frame(cokey=sub(pattern = 'soil_', replacement = '', x = cokeys_to_investigate), fc_threshold=fc_thresholds)
trafficability_trial <- trafficability_trial[match(df_results$cokey, names(trafficability_trial))]
match_test <- sapply(seq_along(df_results$cokey), function(i) names(trafficability_trial)[i] == df_results$cokey[i])
all(match_test)
theta_seq <- as.character(seq(0.95, 0.7, -0.01))
df_results$result_opt2 <- sapply(seq_along(df_results$fc_threshold), function(i) {trafficability_trial[[i]][as.character(df_results$fc_threshold[i])==theta_seq,]})
df_results
original_results <- read.csv('C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/Oct2020test/tables/trafficability_0_10cm_rosetta5_0_10cm_properties.csv', stringsAsFactors = FALSE)
df_results$result_opt2_orig <- original_results$result_opt2[match(df_results$cokey, original_results$cokey)]
df_results
df_results$delta <- df_results$result_opt2 - df_results$result_opt2_orig
summary(df_results$delta)
sum(df_results$delta > 0, na.rm = TRUE)
#scratch
which(names(soil_data)=='soil_10789743') #589
soil_data[[589]]$soil
write_atmos(modelDir = modelDir, flood_date = "2005-03-15", flood_duration = 4, soil_name=names(soil_data)[589], PET = Panoche_ETo_2005)
write_selector(mat_number = soil_data[[589]]$mat_number, flood_duration = 4, modelDir = modelDir, VGs = soil_data[[589]]$VGs, soil_name = names(soil_data)[589], MaxIT = 30)
write_profile(depths=soil_data[[589]]$depths, soil = soil_data[[589]]$soil, mat_number = soil_data[[589]]$mat_number, modelDir = modelDir, soil_name = names(soil_data)[589])

which(names(soil_data)=='soil_12038657') #3072
soil_data[[3072]]$soil
write_atmos(modelDir = modelDir, flood_date = "2005-03-15", flood_duration = 4, soil_name=names(soil_data)[3072], PET = Panoche_ETo_2005)
write_selector(mat_number = soil_data[[3072]]$mat_number, flood_duration = 4, modelDir = modelDir, VGs = soil_data[[3072]]$VGs, soil_name = names(soil_data)[3072], MaxIT = 30)
write_profile(depths=soil_data[[3072]]$depths, soil = soil_data[[3072]]$soil, mat_number = soil_data[[3072]]$mat_number, modelDir = modelDir, soil_name = names(soil_data)[3072])
