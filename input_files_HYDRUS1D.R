# library(aqp)
options(max.print = 15000)
options(width=130)
options()$max.print
options()$width
laptop <- FALSE
if (laptop) {
  workDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/data from Stathis'
  templateDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/Oct2020test/Template/Template'
  modelDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/climate_runs'
  climateDir <- 'C:/Users/smdevine/Desktop/Allowable_Depletion/model_scaffold/run_model/Mar2018'
  prelimDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/Oct2020test/tables/reduce_soils'
} else {
  workDir <- 'C:/Users/smdevine/Desktop/PostDoc/trafficability/data from Stathis'
  templateDir <- 'D:/PostDoc/Trafficability/Oct2020test/Template/Template'
  modelDir <- 'D:/PostDoc/Trafficability/climate_runs'
  climateDir <- 'D:/PostDoc/Trafficability/climate_runs/CIMIS_cell_selection'
  prelimDir <- 'C:/Users/smdevine/Desktop/PostDoc/trafficability/soils_of_interest'
}

#set flood date and duration
# flood_date <- '2005-03-15'
# flood_duration <- 4

#create PET vector for testing general soil differences in trafficability
# ETo <- read.csv(file.path(climateDir, 'SpatialCIMIS.ETo.QCpass.csv'), stringsAsFactors = FALSE)
# head(ETo$dates)
# tail(ETo$dates)
# unique(ETo$year)
# ETo <- ETo[ETo$year %in% 2004:2017,]
# dim(ETo)
# ETo_annual <- aggregate(ETo[,6:ncol(ETo)], by=list(year=ETo$year), FUN = sum)
# dim(ETo_annual)
# ETo_annual$cell_3091
# ETo_annual_means <- apply(ETo_annual[,2:ncol(ETo_annual)], 1, mean)
# names(ETo_annual_means) <- 2004:2017
# ETo_annual_means[order(ETo_annual_means)]
# median(ETo_annual_means) #either 2012 or 2016
# mean(ETo_annual_means) #1393.3, closest to 2009
# #2005=low ET, 2009=avg. ET, 2015=high ET
# 
# quantiles_find <- quantile(ETo_annual[ETo_annual$year==2009,], probs = seq(from=5, to=95, by=5)/100)
# length((quantiles_find))
# 
# ETo_annual_ranked <- as.data.frame(t(ETo_annual[,2:ncol(ETo_annual)]))
# colnames(ETo_annual_ranked) <- paste0('yr_', 2004:2017)
# ETo_annual_ranked <- ETo_annual_ranked[order(ETo_annual_ranked$yr_2009),]
# dim(ETo_annual_ranked)
# rows_to_select <- round(c(0.05, seq(from=0.1, to=0.9, by=0.1), 0.95)*13060, digits = 0)
# rows_to_select
# ETo_annual_ranked_selection <- ETo_annual_ranked[rows_to_select,]
# CIMIS_cells <- rownames(ETo_annual_ranked_selection)
# which(rownames(ETo_annual_ranked)=='cell_170036')/nrow(ETo_annual_ranked) #10876
# rows_to_select/nrow(ETo_annual_ranked)
# ETo_cells_of_interest <- ETo[,c(1:5, which(colnames(ETo) %in% CIMIS_cells))]
# write.csv(ETo_cells_of_interest[ETo_cells_of_interest$year %in% 2004:2017,], file.path(modelDir, 'CIMIS_cell_selection', 'ETo_cells_of_interest.csv'), row.names=FALSE)
ETo_cells_of_interest <- read.csv(file.path(climateDir, 'ETo_cells_of_interest.csv'), stringsAsFactors = FALSE)

# tapply(ETo$cell_170036, ETo$year, sum)
# Panoche_ETo_170036 <- ETo[, c('dates', 'cell_170036')]
# Panoche_ETo_2005 <- Panoche_ETo_170036[which(Panoche_ETo_170036$dates=='01_01_2005'):which(Panoche_ETo_170036$dates=='12_31_2005'),2] / 10 #to convert from cm to mm

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

#trim database according to profiles of interest from preliminary analysis
prelim_results <- read.csv(file.path(prelimDir, "soils_modeled_revised_QCpass_Oct2020.csv"), stringsAsFactors = FALSE)
dim(prelim_results)
length(unique(prelim_results$cokey)) #2911
length(unique(mod_database$cokey)) #5685
mod_database <- mod_database[mod_database$cokey %in% prelim_results$cokey, ]
length(unique(mod_database$cokey)) #2911

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

#previous work that removed cemented soils and soils without a surface horizon
# cemented_soils <- sapply(unique_cokeys, function(x) {ifelse(sum(mod_database$Roseta.model[mod_database$cokey==x]==1)>0, TRUE, FALSE)}, USE.NAMES = TRUE)
# 
# mod_database_aqp <- mod_database
# depths(mod_database_aqp) <- cokey ~ hzn_top + hzn_bot
# depth_logic <- checkHzDepthLogic(mod_database_aqp)
# head(depth_logic)
# sum(!depth_logic$valid) #23 are not valid
# sum(depth_logic$missingDepth) #0
# depth_logic_T <- depth_logic[!depth_logic$valid,]
# depth_logic_T
# for (i in 1:nrow(depth_logic_T)) {
#   print(mod_database[mod_database$cokey==depth_logic_T$cokey[i], c('hzn_top', 'hzn_bot')])
# }  #all have lower horizon top equal to bottom 
# 
# head(cemented_soils)
# sum(cemented_soils) #173
soil_data <- lapply(unique_cokeys, function(x) {soil_by_cokey(x)})
length(soil_data) #was 5685, now 2911
# soil_data <- soil_data[!cemented_soils]
# length(soil_data) #5512, now has cemented soils removed
# 
# #find soils with missing surface horizon
# soils_w_topsoil <- sapply(unique_cokeys, function(x) {0 %in% mod_database$hzn_top[mod_database$cokey==x]}, USE.NAMES = TRUE)
# sum(!soils_w_topsoil) #4 missing topsoil
# names(soils_w_topsoil[!soils_w_topsoil])
# 
# soil_data <- soil_data[!(names(soil_data) %in% names(soils_w_topsoil[!soils_w_topsoil]))]
# length(soil_data) #5509, now has soils with no topsoil removed

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

#modify ATMOSPH.IN file according to flood date and duration; will eventually need to modify PET vector
flood_inputs <- function(vector_input, days_to_flood, application=30.48) {
  # vector_input[as.integer(format(as.Date(start_date), format='%j')):(as.integer(format(as.Date(start_date), format='%j')) + days_to_flood - 1)] <- application #old version running annual model
  vector_input[2:(days_to_flood+1)] <- application
  vector_input
}

write_atmos <- function(modelDir, flood_date, flood_duration, soil_name, PET, climate_rows=100, CIMIS) {
  subDir <- paste0(CIMIS, '_', flood_date)
  flood_yr <- unlist(strsplit(flood_date, '-'))[1]
  if(!dir.exists(file.path(modelDir, subDir))) {
    dir.create(file.path(modelDir, subDir))
  }
  if(!dir.exists(file.path(modelDir, subDir, soil_name))) {
    dir.create(file.path(modelDir, subDir, soil_name))
    dir.create(file.path(modelDir, subDir, soil_name, soil_name))
  }
  atmos_header <- readLines(file.path(templateDir, 'ATMOSPH.IN'), n=9)
  atmos_data <- read.table(file.path(templateDir, 'ATMOSPH.IN'), header = FALSE, nrow=climate_rows, skip=9)
  atmos_data$V2 <- 0
  atmos_data$V2 <- flood_inputs(atmos_data$V2, days_to_flood=flood_duration)
  PET_df <- PET[, c('dates', CIMIS)]
  PET_df <- PET_df[which(PET_df$dates==paste0('01_01_', flood_yr)):which(PET_df$dates==paste0('12_31_', flood_yr)),2] / 10
  atmos_data$V3 <- PET_df[(as.integer(format(as.Date(flood_date), format='%j'))-1):(as.integer(format(as.Date(flood_date), format='%j'))+climate_rows-2)]
  atmos_colnames <- unlist(strsplit(atmos_header[9], ' '))
  atmos_colnames <- atmos_colnames[sapply(atmos_colnames, function(x) x!='')]
  atmos_data$V9 <- ""
  colnames(atmos_data) <- atmos_colnames
  atmos_header <- atmos_header[-9]
  atmos_output <- file(file.path(modelDir, subDir, soil_name, soil_name, 'ATMOSPH.IN'), 'w')
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
# write_atmos(modelDir = modelDir, flood_date = '2005-03-15', flood_duration = 4, soil_name = names(soil_data)[1], PET = ETo_cells_of_interest, CIMIS='cell_5248')
# write_atmos(modelDir = modelDir, flood_date = '2015-03-15', flood_duration = 4, soil_name = names(soil_data)[1], PET = ETo_cells_of_interest, CIMIS='cell_5248')

#SELECTOR.IN file creation
#print_start was 4
write_selector <- function(mat_number, flood_duration, modelDir, VGs, soil_name, MaxIT, flood_date, CIMIS) { #flood_date,
  # if(!dir.exists(file.path(modelDir, soil_name))) {
  #   dir.create(file.path(modelDir, soil_name))
  #   dir.create(file.path(modelDir, soil_name, soil_name))
  # }
  subDir <- paste0(CIMIS, '_', flood_date)
  modelDir <- file.path(modelDir, subDir, soil_name, soil_name)
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
  print_times <- c(0.001, seq(from = print_start, to=99.5, by=0.1), 100)
  substring(selector_header2[3], 67, 69) <- as.character(length(print_times)) #inserts correct number of print times
  print_times_df <- as.data.frame(matrix(data=print_times, nrow=158, ncol=6, byrow = TRUE))
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
# write_selector(mat_number = soil_data[[1]]$mat_number, flood_duration = 4, modelDir = modelDir, VGs = soil_data[[1]]$VGs, soil_name = names(soil_data)[1], MaxIT = 30, CIMIS = 'cell_5248', flood_date = '2005-03-15')
# write_selector(mat_number = soil_data[[1]]$mat_number, flood_duration = 4, modelDir = modelDir, VGs = soil_data[[1]]$VGs, soil_name = names(soil_data)[1], MaxIT = 30, CIMIS = 'cell_5248', flood_date = '2005-03-15')

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
# surface_horizons_db <- mod_database[mod_database$hzn_top==0,]
# surface_horizons_db$theta_critA <- surface_horizons_db$teta_r + 0.005
# surface_horizons_db$hCritA <- sapply(1:nrow(surface_horizons_db), function(i) { optimize(find_h_at_theta_fc, interval = c(1000,100000), theta_r = surface_horizons_db$teta_r[i], theta_s = surface_horizons_db$teta_s[i], alpha = surface_horizons_db$alpha..1.cm.[i], n=surface_horizons_db$n[i], theta_fc=surface_horizons_db$theta_critA[i])$minimum})
# summary(surface_horizons_db$hCritA)
# surface_horizons_db$hCritA
# sum(surface_horizons_db$hCritA>99999) #2717
# hist(surface_horizons_db$hCritA)
# surface_horizons_db[surface_horizons_db$hCritA < 10000,]
# vg_theta(theta_r = 0.05872, theta_s = 0.37117, alpha = 0.0335, n=1.34754, h=160)
# theta_fc_est(theta_s = 0.37117, theta_r = 0.05872, Ks=25.381, n=1.34754) #0.2286089 correct according to HYDRUS
# optimize(find_h_at_theta_fc, interval = c(10,10000), theta_r = 0.05872, theta_s = 0.37117, alpha = 0.0335, n=1.34754, theta_fc=0.2286089) #160.129
write_profile <- function(depths, soil, mat_number, modelDir, soil_name, CIMIS, flood_date) {
  subDir <- paste0(CIMIS, '_', flood_date)
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
  profile_output <- file(file.path(modelDir, subDir, soil_name, soil_name, 'PROFILE.DAT'), 'w')
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
# write_profile(depths=soil_data[[1]]$depths, soil = soil_data[[1]]$soil, mat_number = soil_data[[1]]$mat_number, modelDir = modelDir, soil_name = names(soil_data)[1], CIMIS = 'cell_5248', flood_date = '2005-03-15')
# write_profile(depths=soil_data[[1]]$depths, soil = soil_data[[1]]$soil, mat_number = soil_data[[1]]$mat_number, modelDir = modelDir, soil_name = names(soil_data)[1], CIMIS = 'cell_5248', flood_date = '2015-03-15')

#create directory of input data for all soil names that aren't cemented
cellnames <- colnames(ETo_cells_of_interest)[6:ncol(ETo_cells_of_interest)]
date_of_interest <- '2005-01-15'
for (j in 1:length(cellnames)) {
  for(i in 1:length(soil_data)) {
    write_atmos(modelDir = modelDir, flood_date = date_of_interest, flood_duration = 4, soil_name = names(soil_data)[i], PET = ETo_cells_of_interest, CIMIS=cellnames[j])
    write_selector(mat_number = soil_data[[i]]$mat_number, flood_duration = 4, modelDir = modelDir, VGs = soil_data[[i]]$VGs, soil_name = names(soil_data)[i], MaxIT = 30, CIMIS = cellnames[j], flood_date = date_of_interest)
    write_profile(depths=soil_data[[i]]$depths, soil = soil_data[[i]]$soil, mat_number = soil_data[[i]]$mat_number, modelDir = modelDir, soil_name = names(soil_data)[i], CIMIS = cellnames[j], flood_date = date_of_interest)
  }
  #write file paths and run.bat file
  paths_shortcut <- paste0('path', 1:length(soil_data))
  sub_Dir <-  paste0(cellnames[j], '_', date_of_interest)
  if(!dir.exists(file.path(modelDir, sub_Dir, 'paths'))) {
    dir.create(file.path(modelDir, sub_Dir, 'paths'))
  }
  file.copy(from='D:/PostDoc/Trafficability/climate_runs/2005-03-15/cell_174576_2005-03-15/paths/H1D_CALC.EXE', to=file.path(modelDir, sub_Dir, 'paths', 'H1D_CALC.EXE'))
  runbat <- file(file.path(modelDir, sub_Dir, 'paths', 'run.bat'), 'w')
  for (i in seq_along(paths_shortcut)) {
    writeLines(paste('copy', paths_shortcut[i], 'level_01.dir'), con = runbat)
    writeLines('H1D_CALC', con = runbat)
  }
  close(runbat)
  for (i in seq_along(paths_shortcut)) {
    filepth <- file(file.path(modelDir, sub_Dir, 'paths', paths_shortcut[i]), 'w')
    writeLines(file.path('D:/PostDoc/Trafficability/climate_runs', sub_Dir, names(soil_data)[i], names(soil_data)[i], fsep='/'), con = filepth)
    close(filepth)
  }
  print(options()$width)
}
