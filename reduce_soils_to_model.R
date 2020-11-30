library(aqp)
workDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/data from Stathis'
WBDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/Oct2020test/water_blnc_check'
summaryDir2 <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/Oct2020test/summary_v2'
tablesDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/Oct2020test/tables'

#get textural class--function taken from ssurgo_allCA_aggregation_awc_r.R
textural.class.calc <- function(sand, silt, clay) {
  ifelse(is.na(sand) | is.na(silt) | is.na(clay), NA,
         ifelse(sand + silt + clay > 101 |
                  sand + silt + clay < 99, 'proportions do not sum to 100+-1',
                ifelse(silt + 1.5 * clay < 15, 'sand',
                       ifelse(silt + 1.5 * clay >= 15 & silt + 2 * clay < 30, 'loamy sand',
                              ifelse((clay >= 7 & clay < 20 & sand > 52 & silt + 2 * clay >= 30) | 
                                       (clay < 7 & silt < 50 & silt + 2 * clay >= 30), 'sandy loam',
                                     ifelse(clay >= 7 & clay < 27 & silt >=28 & silt < 50 & sand <= 52, 'loam',
                                            ifelse((silt >= 50 & clay >= 12 & clay < 27) | 
                                                     (silt >=50 & silt < 80 & clay < 12), 'silt loam',
                                                   ifelse(silt >= 80 & clay < 12, 'silt',
                                                          ifelse(clay >= 20 & clay < 35 & silt < 28 & sand > 45, 'sandy clay loam',
                                                                 ifelse(clay >= 27 & clay < 40 & sand > 20 & sand <= 45, 'clay loam',
                                                                        ifelse(clay >= 27 & clay < 40 & sand <= 20, 'silty clay loam',
                                                                               ifelse(clay >= 35 & sand > 45, 'sandy clay',
                                                                                      ifelse(clay >= 40 & silt >= 40, 'silty clay',
                                                                                             ifelse(clay >= 40 & sand <= 45 & silt < 40, 'clay',
                                                                                                    'undefined textural class'))))))))))))))
}

#read in FC thresholds
FC_thresholds <- read.csv(file.path(tablesDir, 'trafficability_defs_11_3_20.csv'), stringsAsFactors = FALSE)
FC_thresholds
FC_thresholds$Option_4 <- 0.9
FC_thresholds$Option_4[FC_thresholds$Textural.Class=='clay'] <- 0.8
theta_seq <- as.character(seq(0.95, 0.7, -0.01)) #how results are order in trafficability_0_10cm_v2

#read in soil data
mod_database <- read.csv(file.path(workDir, 'modelling_database.csv'), stringsAsFactors = FALSE, na.strings = '-9.9')
mod_database <- mod_database[which(mod_database$hzn_top<200),] #delete horizons that start at or below 200 cm depth
lapply(mod_database, summary)
colnames(mod_database)
mod_database$clay <- round(mod_database$clay, digits = 1)
mod_database$silt <- round(mod_database$silt, digits = 1)
mod_database$sand <- round(mod_database$sand, digits = 1)
soil_by_cokey <- function(cokey) {
  soil <- mod_database[mod_database$cokey==cokey,]
  mat_number <- nrow(soil)
  soil <- soil[order(soil$hzn_top),]
  VGs <- soil[,c(20:24)]
  VGs$l <- 0.5
  colnames(VGs) <- NULL
  depths <- soil[,8:9]
  depths$hzn_bot[length(depths$hzn_bot)] <- 201
  depths_calc <- depths
  depths_calc$hzn_bot[length(depths_calc$hzn_bot)] <- 200
  uniqueness <- c(soil$hzn_top, soil$hzn_bot, soil$clay, soil$silt, soil$sand)
  list(soil=soil, VGs=VGs, depths=depths, mat_number=mat_number, theta_s=sum(soil$teta_s * (depths_calc$hzn_bot - depths_calc$hzn_top)), uniqueness=uniqueness)
}

unique_cokeys <- unique(mod_database$cokey)
length(unique_cokeys) #5685
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
length(soil_data) #5509 now
head(soil_data)
final_cokeys <- as.integer(gsub('soil_', '', names(soil_data)))

horizons_modeled <- mod_database[mod_database$cokey %in% final_cokeys, ]

#check for horizons with illuviated clay
unique(mod_database$hzn_desgn[grepl('t', mod_database$hzn_desgn)])

#how many unique profiles
soil_data[[1]]$uniqueness
length(soil_data[[1]]$uniqueness)
uniqueness_strings <- sapply(soil_data, function(x) as.character(x$uniqueness))
length(unique(uniqueness_strings)) #3333 unique profiles
unique_profiles <- unique(uniqueness_strings)
soil_data_unique <- soil_data[match(unique_profiles, uniqueness_strings)]
length(soil_data_unique)
names(soil_data) %in% names(soil_data_unique)

#add h and theta @ FC info
theta_fc_est <- function(theta_s, theta_r, Ks, n) {
  (n^(-0.6*(2+log10(Ks))))*(theta_s - theta_r) + theta_r
}
vg_theta <- function(theta_r, theta_s, alpha, n, h) {
  theta_r + (theta_s - theta_r) / (1 + (alpha * h)^n)^(1-1/n)
}
find_h_at_theta_fc <- function(h, theta_r, theta_s, alpha, n, theta_fc){
  abs(theta_fc  - theta_r - ((theta_s - theta_r) / (1 + (alpha * h)^n)^(1-1/n)))
}

horizons_modeled$theta_fc <- theta_fc_est(horizons_modeled$teta_s, horizons_modeled$teta_r, horizons_modeled$Ks..cm.d., horizons_modeled$n)

horizons_modeled$h_fc <- -sapply(1:nrow(horizons_modeled), function(i) { optimize(find_h_at_theta_fc, interval = c(10,10000), theta_r = horizons_modeled$teta_r[i], theta_s = horizons_modeled$teta_s[i], alpha = horizons_modeled$alpha..1.cm.[i], n=horizons_modeled$n[i], theta_fc=horizons_modeled$theta_fc[i])$minimum})
summary(horizons_modeled$h_fc)

horizons_modeled$theta_wp_est <- vg_theta(theta_r = horizons_modeled$teta_r, theta_s = horizons_modeled$teta_s, alpha = horizons_modeled$alpha..1.cm., n = horizons_modeled$n, h = (1500*10.1972))
summary(horizons_modeled$theta_wp_est)
hist(horizons_modeled$theta_wp_est)
summary(horizons_modeled$wr_15b)
summary(horizons_modeled$wr_15b - horizons_modeled$theta_wp_est)
horizons_modeled$awc_est <- horizons_modeled$theta_fc - horizons_modeled$theta_wp_est
summary(horizons_modeled$awc_est)
horizons_modeled$awc_est_v2 <- horizons_modeled$wr_13b - horizons_modeled$wr_15b
summary(horizons_modeled$awc_est_v2)
sum(horizons_modeled$awc_est_v2==0, na.rm = TRUE) #1
sum(horizons_modeled$awc_est_v2<0.04, na.rm = TRUE) #470
horizons_modeled[which(horizons_modeled$awc_est_v2==0), ]
horizons_modeled$ps_sums <- horizons_modeled$sand + horizons_modeled$silt + horizons_modeled$clay
summary(horizons_modeled$ps_sums)
horizons_modeled$textural_class <- textural.class.calc(sand = horizons_modeled$sand, silt = horizons_modeled$silt, clay = horizons_modeled$clay)
table(horizons_modeled$textural_class)
tapply(horizons_modeled$awc_est, horizons_modeled$textural_class, summary)
tapply(horizons_modeled$awc_est_v2, horizons_modeled$textural_class, summary)

wtd.mean_v2 <- function(x, y) {
  # use horizon thickness as a weight
  thick <- x$hzn_bot - x$hzn_top
  # compute the weighted mean, accounting for the possibility of missing data
  m <- weighted.mean(horizons(x)[[y]], w=thick, na.rm=FALSE)
  m
}

depth_X_properties <- function(horizon_SPC, upper_depth=0, depth, vars_of_interest = c('clay', 'silt', 'sand', 'estimated_oc', 'db_13b', 'wr_13b', 'wr_15b', 'Ks..cm.d.', 'K0.cm.d.', 'teta_r', 'teta_s', 'alpha..1.cm.', 'n', 'theta_fc', 'h_fc', 'theta_wp_est', 'awc_est', 'awc_est_v2', 'Roseta.model'), varnames = c('clay', 'silt', 'sand', 'oc', 'bd_13b', 'theta_0.33b', 'theta_15b', 'ksat', 'f_ksat', 'theta_r', 'theta_s', 'alpha', 'n', 'theta_fc', 'h_fc', 'theta_wp', 'awc', 'awc_v2', 'Rosetta.model')) { #lep is linear extensibility
  columnames <- paste0(varnames, '_', depth, 'cm')
  print(cbind(vars_of_interest, columnames)) #show that it's all lined up
  assign("depth", depth, envir = .GlobalEnv) #this necessary because slice can't find the variable otherwise
  assign("upper_depth", upper_depth, envir = .GlobalEnv)
  sliced_SPC <- slice(horizon_SPC, upper_depth:(depth-1) ~ .) #depth was '0:depth' in previous version
  horizons(sliced_SPC) <- horizons(sliced_SPC)[order(as.integer(sliced_SPC$cokey), sliced_SPC$hzn_top),] #because integer codes are coerced to character by slice with the sliced data.frame re-ordered contrary to order of site-level data 
  stopifnot(unique(sliced_SPC$cokey)==site(sliced_SPC)$cokey)
  for (i in seq_along(vars_of_interest)) {
    s <- site(sliced_SPC)
    s[[columnames[i]]] <- profileApply(sliced_SPC, FUN = wtd.mean_v2, y=vars_of_interest[i])
    site(sliced_SPC) <- s
  }
  s <- site(sliced_SPC)
  rm(depth, envir = .GlobalEnv) #because we had to put it there earlier
  s
}

#create spc object
depths(horizons_modeled) <- cokey ~ hzn_top + hzn_bot

soils_modeled_10cm <- depth_X_properties(horizon_SPC = horizons_modeled, depth = 10)
dim(soils_modeled_10cm)


#analysis with a 0-10 cm soil property summary
soils_modeled_10cm$texture.sums <- soils_modeled_10cm$clay_10cm + soils_modeled_10cm$silt_10cm + soils_modeled_10cm$sand_10cm
sum(is.na(soils_modeled_10cm$texture.sums)) #0 are NA
sum(soils_modeled_10cm$texture.sums > 100.1 | soils_modeled_10cm$texture.sums < 99.9 ) #no hay problemsas

soils_modeled_10cm$textural_class <- textural.class.calc(sand = soils_modeled_10cm$sand_10cm, silt = soils_modeled_10cm$silt_10cm, clay = soils_modeled_10cm$clay_10cm)
table(soils_modeled_10cm$textural_class)

#condense trafficability results from HYDRUS runs using summaryDir2
#for summary_v2 folder
extract_0_10cm_trafficability_v2 <- function(x) {
  results <- read.csv(file.path(summaryDir2, x), row.names=1, na.strings = 'numeric(0)')
}
fnames_results <- list.files(summaryDir2, full.names = FALSE, recursive = FALSE)
head(fnames_results)
tail(fnames_results)
cokeys_modeled <- sub(pattern = 'soil_', replacement = '', x = fnames_results)
cokeys_modeled <- sub(pattern = '_results.csv', replacement = '', x = cokeys_modeled)
head(cokeys_modeled)
tail(cokeys_modeled) 
length(cokeys_modeled)#5494 because 15 soils produced error messages
trafficability_0_10cm_v2 <- lapply(fnames_results, extract_0_10cm_trafficability_v2)
length(trafficability_0_10cm_v2)
names(trafficability_0_10cm_v2) <- cokeys_modeled

#approach to summarize data on summaryDir2
dim(soils_modeled_10cm) #5509
soils_modeled_10cm <- soils_modeled_10cm[soils_modeled_10cm$cokey %in% names(trafficability_0_10cm_v2),]
dim(soils_modeled_10cm) #5494
soils_modeled_10cm$traff_def_opt2 <- FC_thresholds$Option_2[match(soils_modeled_10cm$textural_class, FC_thresholds$Textural.Class)]

#exclude silt cokeys
silt_cokeys <- soils_modeled_10cm$cokey[soils_modeled_10cm$textural_class=='silt'] #because no plasticity based guide for these
trafficability_0_10cm_v2 <- trafficability_0_10cm_v2[-which(names(trafficability_0_10cm_v2) %in% silt_cokeys)]
length(trafficability_0_10cm_v2)
soils_modeled_10cm <- soils_modeled_10cm[!soils_modeled_10cm$textural_class=='silt',]
dim(soils_modeled_10cm)
match_test <- sapply(seq_along(soils_modeled_10cm$cokey), function(i) names(trafficability_0_10cm_v2)[i] == soils_modeled_10cm$cokey[i])
all(match_test) #everything in order

soils_modeled_10cm$result_opt2 <- sapply(seq_along(soils_modeled_10cm$traff_def_opt2), function(i) {trafficability_0_10cm_v2[[i]][as.character(soils_modeled_10cm$traff_def_opt2[i])==theta_seq,]})
summary(soils_modeled_10cm$result_opt2) #139 NA's

#get unique textural classes
textural_classes <- unique(soils_modeled_10cm$textural_class)
names(textural_classes) <- textural_classes

#add water balance QC info
water_blnc_results <- read.csv(file.path(WBDir, 'water_blnc_results.csv'), stringsAsFactors = FALSE)
sum(water_blnc_results$runoff_PF_cm > 0, na.rm = TRUE)
sum(water_blnc_results$infiltration_PF_cm > 0, na.rm = TRUE) #7
soils_modeled_10cm$wb_rel_pct <- water_blnc_results$water_balance_rel_pct[match(soils_modeled_10cm$cokey, water_blnc_results$cokey)]
soils_modeled_10cm$wb_PF_cm <- water_blnc_results$water_balance_PF_cm[match(soils_modeled_10cm$cokey, water_blnc_results$cokey)]
summary(soils_modeled_10cm$wb_PF_cm)
soils_modeled_10cm$wb_PF_cm <- abs(soils_modeled_10cm$wb_PF_cm)
soils_modeled_10cm$water_input_error <- water_blnc_results$water_input_error[match(soils_modeled_10cm$cokey, water_blnc_results$cokey)]
soils_modeled_10cm$theta_pf <- water_blnc_results$thetaS_ck[match(soils_modeled_10cm$cokey, water_blnc_results$cokey)] #thetaS_max is max from time-series; thetaS_ck is from day 5 at end of flooding event
soils_modeled_10cm$theta_s_profile_cm <- sapply(paste0('soil_', soils_modeled_10cm$cokey), function(x) soil_data[[match(x, names(soil_data))]]$theta_s)
soils_modeled_10cm$pf_theta_sat_check <- soils_modeled_10cm$theta_s_profile_cm - soils_modeled_10cm$theta_pf
summary(soils_modeled_10cm$pf_theta_sat_check)
soils_modeled_10cm$sat_pct_pf <- round(100*(soils_modeled_10cm$theta_s_profile_cm - soils_modeled_10cm$pf_theta_sat_check) / soils_modeled_10cm$theta_s_profile_cm, 2)


#look at unique soil profile instances
sum(soils_modeled_10cm$cokey %in% sub('soil_', '', names(soil_data_unique))) #3319
soils_modeled_10cm_revised <- soils_modeled_10cm[soils_modeled_10cm$cokey %in% sub('soil_', '', names(soil_data_unique)),]

# soils_modeled_10cm$properties <- paste(soils_modeled_10cm$clay_10cm, soils_modeled_10cm$silt_10cm, soils_modeled_10cm$sand_10cm, soils_modeled_10cm$bd_13b_10cm, soils_modeled_10cm$theta_0.33b_10cm, soils_modeled_10cm$theta_15b_10cm, soils_modeled_10cm$result_opt2, sep = ',')
# soils_modeled_10cm$properties <- paste(soils_modeled_10cm$clay_10cm, soils_modeled_10cm$silt_10cm, soils_modeled_10cm$sand_10cm, soils_modeled_10cm$result_opt2, sep = ',')
# length(unique(soils_modeled_10cm$properties)) #3144 compared to 3107 in 0-30 cm dataset for rosetta 5 only; only 2798 if only consider sand, silt, and clay and all rosetta models
# length(unique(soils_modeled_10cm$properties[soils_modeled_10cm$Rosetta.model_10cm==5]))
# unique_properties <- unique(soils_modeled_10cm$properties)
# soils_modeled_10cm[unique_properties[4]==soils_modeled_10cm$properties,]
# 
# soils_modeled_10cm_revised <- soils_modeled_10cm[match(unique_properties, soils_modeled_10cm$properties), ]
dim(soils_modeled_10cm_revised) #3319 were derived from unique profiles
table(soils_modeled_10cm_revised$textural_class[!is.na(soils_modeled_10cm_revised$result_opt2)])
lapply(textural_classes, function(x) summary(soils_modeled_10cm_revised$result_opt2[soils_modeled_10cm_revised$textural_class==x]))
lapply(textural_classes, function(x) summary(soils_modeled_10cm$result_opt2[soils_modeled_10cm$textural_class==x]))
lapply(textural_classes, function(x) summary(soils_modeled_10cm$result_opt2[soils_modeled_10cm$textural_class==x])[3] - summary(soils_modeled_10cm_revised$result_opt2[soils_modeled_10cm_revised$textural_class==x])[3]) #only silty clay stats are changed appreciably
lapply(textural_classes, function(x) summary(soils_modeled_10cm_revised$result_opt2[soils_modeled_10cm_revised$textural_class==x & soils_modeled_10cm_revised$Rosetta.model_10cm==5])[3] - summary(soils_modeled_10cm_revised$result_opt2[soils_modeled_10cm_revised$textural_class==x])[3])
lapply(textural_classes, function(x) summary(soils_modeled_10cm_revised$result_opt2[soils_modeled_10cm_revised$textural_class==x & soils_modeled_10cm_revised$Rosetta.model_10cm==5])[3] - summary(soils_modeled_10cm_revised$result_opt2[soils_modeled_10cm_revised$textural_class==x & soils_modeled_10cm_revised$Rosetta.model_10cm==2])[3])

#look at water balance issues by textural class
summary(soils_modeled_10cm_revised$sat_pct_pf)
sum(soils_modeled_10cm_revised$sat_pct_pf < 99, na.rm = TRUE) #672 less than 99% saturated
sum(soils_modeled_10cm_revised$sat_pct_pf < 98, na.rm = TRUE) #550 less than 98% saturated
sum(soils_modeled_10cm_revised$sat_pct_pf < 95, na.rm = TRUE) #399 less than 98% saturated
sum(soils_modeled_10cm_revised$sat_pct_pf < 90, na.rm = TRUE) #235 less than 90% saturated
tapply(soils_modeled_10cm_revised$sat_pct_pf[!is.na(soils_modeled_10cm_revised$result_opt2)], soils_modeled_10cm_revised$textural_class[!is.na(soils_modeled_10cm_revised$result_opt2)], summary)
tapply(soils_modeled_10cm_revised$sat_pct_pf[!is.na(soils_modeled_10cm_revised$result_opt2)], soils_modeled_10cm_revised$textural_class[!is.na(soils_modeled_10cm_revised$result_opt2)], function(x) sum(x < 95, na.rm = TRUE))
tapply(soils_modeled_10cm_revised$sat_pct_pf[!is.na(soils_modeled_10cm_revised$result_opt2)], soils_modeled_10cm_revised$textural_class[!is.na(soils_modeled_10cm_revised$result_opt2)], function(x) sum(x < 90, na.rm = TRUE))
tapply(soils_modeled_10cm_revised$result_opt2[!is.na(soils_modeled_10cm_revised$result_opt2)], soils_modeled_10cm_revised$textural_class[!is.na(soils_modeled_10cm_revised$result_opt2)], function(x) median(x))
tapply(soils_modeled_10cm_revised$result_opt2[!is.na(soils_modeled_10cm_revised$result_opt2)], soils_modeled_10cm_revised$textural_class[!is.na(soils_modeled_10cm_revised$result_opt2)], function(x) median(x))
tapply(soils_modeled_10cm_revised$result_opt2[!is.na(soils_modeled_10cm_revised$result_opt2)], soils_modeled_10cm_revised$textural_class[!is.na(soils_modeled_10cm_revised$result_opt2)], function(x) sd(x))
tapply(soils_modeled_10cm_revised$result_opt2[!is.na(soils_modeled_10cm_revised$result_opt2)], soils_modeled_10cm_revised$textural_class[!is.na(soils_modeled_10cm_revised$result_opt2)], function(x) median(x))
tapply(soils_modeled_10cm_revised$result_opt2[!is.na(soils_modeled_10cm_revised$result_opt2)], soils_modeled_10cm_revised$textural_class[!is.na(soils_modeled_10cm_revised$result_opt2)], function(x) sd(x))

soils_modeled_10cm_revised[soils_modeled_10cm_revised$cokey==20199,]
tapply(soils_modeled_10cm_revised$wb_rel_pct[!is.na(soils_modeled_10cm_revised$result_opt2)], soils_modeled_10cm_revised$textural_class[!is.na(soils_modeled_10cm_revised$result_opt2)], summary)
tapply(soils_modeled_10cm_revised$wb_rel_pct[!is.na(soils_modeled_10cm_revised$result_opt2)], soils_modeled_10cm_revised$textural_class[!is.na(soils_modeled_10cm_revised$result_opt2)], function(x) sum(x > 5, na.rm=TRUE))
tapply(soils_modeled_10cm_revised$wb_rel_pct[!is.na(soils_modeled_10cm_revised$result_opt2)], soils_modeled_10cm_revised$textural_class[!is.na(soils_modeled_10cm_revised$result_opt2)], function(x) sum(x > 3, na.rm=TRUE))
tapply(soils_modeled_10cm_revised$wb_PF_cm[!is.na(soils_modeled_10cm_revised$result_opt2)], soils_modeled_10cm_revised$textural_class[!is.na(soils_modeled_10cm_revised$result_opt2)], function(x) sum(x > 0.5, na.rm=TRUE))
tapply(soils_modeled_10cm_revised$wb_PF_cm[!is.na(soils_modeled_10cm_revised$result_opt2)], soils_modeled_10cm_revised$textural_class[!is.na(soils_modeled_10cm_revised$result_opt2)], function(x) sum(x > 0.25, na.rm=TRUE)) #10

tapply(soils_modeled_10cm_revised$result_opt2[which(soils_modeled_10cm_revised$sat_pct_pf >=90 & soils_modeled_10cm_revised$wb_PF_cm < 0.25)], soils_modeled_10cm_revised$textural_class[which(soils_modeled_10cm_revised$sat_pct_pf >=90 & soils_modeled_10cm_revised$wb_PF_cm < 0.25)], function(x) median(x, na.rm = TRUE))
tapply(soils_modeled_10cm_revised$result_opt2[which(soils_modeled_10cm_revised$sat_pct_pf >=90 & soils_modeled_10cm_revised$wb_PF_cm < 0.25)], soils_modeled_10cm_revised$textural_class[which(soils_modeled_10cm_revised$sat_pct_pf >=90 & soils_modeled_10cm_revised$wb_PF_cm < 0.25)], function(x) sd(x, na.rm = TRUE))
table(soils_modeled_10cm_revised$textural_class[which(soils_modeled_10cm_revised$sat_pct_pf >=90 & soils_modeled_10cm_revised$wb_PF_cm < 0.25 & !is.na(soils_modeled_10cm_revised$result_opt2))])

tapply(soils_modeled_10cm_revised$result_opt2[which(soils_modeled_10cm_revised$sat_pct_pf >=80 & soils_modeled_10cm_revised$wb_PF_cm < 0.25 & soils_modeled_10cm_revised$wb_rel_pct <= 3)], soils_modeled_10cm_revised$textural_class[which(soils_modeled_10cm_revised$sat_pct_pf >=80 & soils_modeled_10cm_revised$wb_PF_cm < 0.25 & soils_modeled_10cm_revised$wb_rel_pct <= 3)], function(x) median(x, na.rm = TRUE))
tapply(soils_modeled_10cm_revised$result_opt2[which(soils_modeled_10cm_revised$sat_pct_pf >=80 & soils_modeled_10cm_revised$wb_PF_cm < 0.25 & soils_modeled_10cm_revised$wb_rel_pct <= 3)], soils_modeled_10cm_revised$textural_class[which(soils_modeled_10cm_revised$sat_pct_pf >=80 & soils_modeled_10cm_revised$wb_PF_cm < 0.25 & soils_modeled_10cm_revised$wb_rel_pct <= 3)], function(x) sd(x, na.rm = TRUE))
tapply(soils_modeled_10cm_revised$result_opt2[which(soils_modeled_10cm_revised$sat_pct_pf >=80 & soils_modeled_10cm_revised$wb_PF_cm < 0.25 & soils_modeled_10cm_revised$wb_rel_pct <= 3)], soils_modeled_10cm_revised$textural_class[which(soils_modeled_10cm_revised$sat_pct_pf >=80 & soils_modeled_10cm_revised$wb_PF_cm < 0.25 & soils_modeled_10cm_revised$wb_rel_pct <= 3)], function(x) max(x, na.rm = TRUE))
tapply(soils_modeled_10cm_revised$result_opt2[which(soils_modeled_10cm_revised$sat_pct_pf >=80 & soils_modeled_10cm_revised$wb_PF_cm < 0.25 & soils_modeled_10cm_revised$wb_rel_pct <= 3)], soils_modeled_10cm_revised$textural_class[which(soils_modeled_10cm_revised$sat_pct_pf >=80 & soils_modeled_10cm_revised$wb_PF_cm < 0.25 & soils_modeled_10cm_revised$wb_rel_pct <= 3)], function(x) min(x, na.rm = TRUE))
table(soils_modeled_10cm_revised$textural_class[which(soils_modeled_10cm_revised$sat_pct_pf >=80 & soils_modeled_10cm_revised$wb_PF_cm < 0.25 & soils_modeled_10cm_revised$wb_rel_pct <= 3 & !is.na(soils_modeled_10cm_revised$result_opt2))])
lapply(textural_classes, function(x) hist(soils_modeled_10cm_revised$result_opt2[which(soils_modeled_10cm_revised$sat_pct_pf >=80 & soils_modeled_10cm_revised$wb_PF_cm < 0.25 & soils_modeled_10cm_revised$wb_rel_pct <= 3 & soils_modeled_10cm_revised$textural_class==x)], main = x))


tapply(soils_modeled_10cm_revised$result_opt2[which(soils_modeled_10cm_revised$sat_pct_pf >=80 & soils_modeled_10cm_revised$wb_PF_cm < 0.25 & soils_modeled_10cm_revised$wb_rel_pct <= 5)], soils_modeled_10cm_revised$textural_class[which(soils_modeled_10cm_revised$sat_pct_pf >=80 & soils_modeled_10cm_revised$wb_PF_cm < 0.25 & soils_modeled_10cm_revised$wb_rel_pct <= 5)], function(x) median(x, na.rm = TRUE))
tapply(soils_modeled_10cm_revised$result_opt2[which(soils_modeled_10cm_revised$sat_pct_pf >=80 & soils_modeled_10cm_revised$wb_PF_cm < 0.25 & soils_modeled_10cm_revised$wb_rel_pct <= 5)], soils_modeled_10cm_revised$textural_class[which(soils_modeled_10cm_revised$sat_pct_pf >=80 & soils_modeled_10cm_revised$wb_PF_cm < 0.25 & soils_modeled_10cm_revised$wb_rel_pct <= 5)], function(x) sd(x, na.rm = TRUE))
tapply(soils_modeled_10cm_revised$result_opt2[which(soils_modeled_10cm_revised$sat_pct_pf >=80 & soils_modeled_10cm_revised$wb_PF_cm < 0.25 & soils_modeled_10cm_revised$wb_rel_pct <= 5)], soils_modeled_10cm_revised$textural_class[which(soils_modeled_10cm_revised$sat_pct_pf >=80 & soils_modeled_10cm_revised$wb_PF_cm < 0.25 & soils_modeled_10cm_revised$wb_rel_pct <= 5)], function(x) max(x, na.rm = TRUE))
tapply(soils_modeled_10cm_revised$result_opt2[which(soils_modeled_10cm_revised$sat_pct_pf >=80 & soils_modeled_10cm_revised$wb_PF_cm < 0.25 & soils_modeled_10cm_revised$wb_rel_pct <= 5)], soils_modeled_10cm_revised$textural_class[which(soils_modeled_10cm_revised$sat_pct_pf >=80 & soils_modeled_10cm_revised$wb_PF_cm < 0.25 & soils_modeled_10cm_revised$wb_rel_pct <= 5)], function(x) min(x, na.rm = TRUE))
table(soils_modeled_10cm_revised$textural_class[which(soils_modeled_10cm_revised$sat_pct_pf >=80 & soils_modeled_10cm_revised$wb_PF_cm < 0.25 & soils_modeled_10cm_revised$wb_rel_pct <= 5 & !is.na(soils_modeled_10cm_revised$result_opt2))])
lapply(textural_classes, function(x) hist(soils_modeled_10cm_revised$result_opt2[which(soils_modeled_10cm_revised$sat_pct_pf >=80 & soils_modeled_10cm_revised$wb_PF_cm < 0.25 & soils_modeled_10cm_revised$wb_rel_pct <= 5 & soils_modeled_10cm_revised$textural_class==x)], main = x))

table(soils_modeled_10cm_revised$textural_class[which(soils_modeled_10cm_revised$sat_pct_pf >=80 & soils_modeled_10cm_revised$wb_PF_cm < 0.25 & !is.na(soils_modeled_10cm_revised$result_opt2))])


#look at extremes by texture
x <- 'silty clay'
crit_val <- 80
soils_modeled_10cm_revised[which(soils_modeled_10cm_revised$sat_pct_pf >=80 & soils_modeled_10cm_revised$wb_PF_cm < 0.25 & soils_modeled_10cm_revised$wb_rel_pct <= 3 & soils_modeled_10cm_revised$textural_class==x & soils_modeled_10cm_revised$result_opt2 > crit_val),]
soils_modeled_10cm_revised[which(soils_modeled_10cm_revised$textural_class==x & soils_modeled_10cm_revised$result_opt2 > crit_val),]

#create a QC limited soil database
soils_modeled_10cm_revised_QC <- soils_modeled_10cm_revised[which(soils_modeled_10cm_revised$wb_PF_cm < 0.25 & soils_modeled_10cm_revised$wb_rel_pct <= 3), ]
dim(soils_modeled_10cm_revised_QC) #2920

#create summary based on soil names
soils_modeled_10cm$soil_name <- mod_database$taxonname[match(soils_modeled_10cm$cokey, mod_database$cokey)]
soil_name_counts <- table(soils_modeled_10cm$soil_name)[order(table(soils_modeled_10cm$soil_name), decreasing = TRUE)]
soil_name_mas30 <- soil_name_counts[soil_name_counts >= 30]

summarize_by_names <- function(traff_def) {
  soilname_summary <- do.call(cbind, lapply(names(soil_name_mas30), function(x) as.matrix(summary(soils_modeled_10cm[[traff_def]][soils_modeled_10cm$soil_name==x]))[1:6,]))
  colnames(soilname_summary) <- names(soil_name_mas30)
  soilname_summary <- rbind(soilname_summary, n=as.integer(sapply(names(soil_name_mas30), function(x) sum(!is.na(soils_modeled_10cm[[traff_def]][soils_modeled_10cm$soil_name==x])))))
  soilname_summary <- rbind(soilname_summary, mean_clay=sapply(names(soil_name_mas30), function(x) mean(soils_modeled_10cm$clay_10cm[soils_modeled_10cm$soil_name==x], na.rm=TRUE)))
  soilname_summary <- rbind(soilname_summary, texture=sapply(names(soil_name_mas30), function(x) paste(unique(soils_modeled_10cm$textural_class[soils_modeled_10cm$soil_name==x]), collapse = ', ')))
  print(soilname_summary)
  write.csv(soilname_summary, file.path(tablesDir, 'reduce_soils', paste0('summary_by_soilname_0_10cm_properties_', traff_def, '.csv')), row.names = TRUE)
}
summarize_by_names('result_opt2')


summarize_by_texture <- function(df, traff_def, fname) {
  textural_class_summary <- do.call(cbind, lapply(unique(df$textural_class), function(x) as.matrix(summary(df[[traff_def]][df$textural_class==x]))[1:6,]))
  colnames(textural_class_summary) <- unique(df$textural_class)
  textural_class_summary <- rbind(textural_class_summary, n=as.integer(sapply(unique(df$textural_class), function(x) sum(!is.na(df[[traff_def]][df$textural_class==x])))))
  textural_class_summary <- rbind(textural_class_summary, mean_clay=sapply(unique(df$textural_class), function(x) mean(df$clay_10cm[df$textural_class==x], na.rm=TRUE)))
  textural_class_summary <- rbind(textural_class_summary, Not_Determined=as.integer(sapply(unique(df$textural_class), function(x) sum(is.na(df[[traff_def]][df$textural_class==x])))))
  colnames(textural_class_summary)[order(textural_class_summary[8,])]
textural_class_summary <- textural_class_summary[,order(textural_class_summary[8,])]
  print(textural_class_summary)
  write.csv(textural_class_summary, file.path(tablesDir, 'reduce_soils', paste0('summary_by_textural_class_0_10cm_', traff_def, '_', fname, '.csv')), row.names = TRUE)
}
summarize_by_texture(df=soils_modeled_10cm_revised, traff_def = 'result_opt2', fname = 'unique')
summarize_by_texture(df=soils_modeled_10cm, traff_def = 'result_opt2', fname = 'all')
summarize_by_texture(df=soils_modeled_10cm_revised_QC, traff_def = 'result_opt2', fname = 'unique_QC')

#plot results vs clay
plot(soils_modeled_10cm_revised_QC$clay_10cm, soils_modeled_10cm_revised_QC$result_opt2)

#make a vioplot




lapply(1:length(textural_classes), function(i) {
  hist(soils_modeled_10cm_revised_QC$result_opt2[soils_modeled_10cm_revised_QC$textural_class==textural_classes[i]], main=textural_classes[i])
})

quantile_breaks_result2 <- lapply(textural_classes, function(x) {
  quantile(soils_modeled_10cm_revised_QC$result_opt2[soils_modeled_10cm_revised_QC$textural_class==x], probs = c(0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99), na.rm = TRUE)
})
quantile_breaks_result2 <- rbind(as.data.frame(quantile_breaks_result2), sd=sapply(textural_classes, function(x) {sd(soils_modeled_10cm_revised_QC$result_opt2[soils_modeled_10cm_revised_QC$textural_class==x], na.rm=TRUE)}), mean=sapply(textural_classes, function(x) {mean(soils_modeled_10cm_revised_QC$result_opt2[soils_modeled_10cm_revised_QC$textural_class==x], na.rm=TRUE)}))
quantile_breaks_result2 <- quantile_breaks_result2[,order(quantile_breaks_result2[6,])]
write.csv(quantile_breaks_result2, file.path(tablesDir, 'reduce_soils', 'quantiles_result2_unique.csv'), row.names = TRUE)

#write reduced database results to file
dim(soils_modeled_10cm_revised_QC) #2920
soils_modeled_10cm_revised_QC <- soils_modeled_10cm_revised_QC[!is.na(soils_modeled_10cm_revised_QC$result_opt2),]
dim(soils_modeled_10cm_revised_QC) #2911
write.csv(soils_modeled_10cm_revised_QC, file.path(tablesDir, 'reduce_soils', 'soils_modeled_revised_QCpass_Oct2020.csv'), row.names = FALSE)
