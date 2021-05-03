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
  templateDir <- 'D:/PostDoc/Trafficability/climate_runs/Template/Template'
  modelDir <- 'D:/PostDoc/Trafficability/climate_runs'
  climateDir <- 'D:/PostDoc/Trafficability/climate_runs/CIMIS_cell_selection'
  prelimDir <- 'C:/Users/smdevine/Desktop/PostDoc/trafficability/soils_of_interest'
}

#read in ETo data
ETo_cells_of_interest <- read.csv(file.path(climateDir, 'ETo_cells_of_interest.csv'), stringsAsFactors = FALSE)

#modify soil profile according to VG params
#test run with Hanford entry from modelling_database.R
mod_database <- read.csv(file.path(workDir, 'modelling_database.csv'), stringsAsFactors = FALSE)
mod_database <- mod_database[which(mod_database$hzn_top<200),] #delete horizons that start at or below 200 cm depth

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

soil_data <- lapply(unique_cokeys, function(x) {soil_by_cokey(x)})
length(soil_data) #was 5685, now 2911

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

#create directory of input data for all soil names that aren't cemented
cellnames <- colnames(ETo_cells_of_interest)[6:ncol(ETo_cells_of_interest)]
cellnames
date_of_interest <- '2005-04-15'

#error file to select
# error_soils <- read.csv(file.path(modelDir, 'errors', 'H1D_errors_cell_168006_2005-01-15.csv'), stringsAsFactors = FALSE)
# length(error_soils$soilnames)
# soil_data <- soil_data[names(soil_data) %in% error_soils$soilnames]
# length(soil_data)
# cellnames[7]
# j <- 7
for (j in 1:length(cellnames)) {
  options(max.print = 15000)
  options(width=130)
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