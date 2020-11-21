library(aqp)
workDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/data from Stathis'
summaryDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/Oct2020test/summary'
summaryDir2 <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/Oct2020test/summary_v2'
tablesDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/Oct2020test/tables'

#read in soil data
mod_database <- read.csv(file.path(workDir, 'modelling_database.csv'), stringsAsFactors = FALSE, na.strings = '-9.9')
mod_database <- mod_database[which(mod_database$hzn_top<200),] #delete horizons that start at or below 200 cm depth
lapply(mod_database, summary)

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

final_cokeys <- as.integer(gsub('soil_', '', names(soil_data)))

horizons_modeled <- mod_database[mod_database$cokey %in% final_cokeys, ]


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

wtd.mean_v2 <- function(x, y) {
  # use horizon thickness as a weight
  thick <- x$hzn_bot - x$hzn_top
  # compute the weighted mean, accounting for the possibility of missing data
  m <- weighted.mean(horizons(x)[[y]], w=thick, na.rm=FALSE)
  m
}

depth_X_properties <- function(horizon_SPC, upper_depth=0, depth, vars_of_interest = c('clay', 'silt', 'sand', 'estimated_oc', 'db_13b', 'wr_13b', 'wr_15b', 'Ks..cm.d.', 'K0.cm.d.', 'teta_r', 'teta_s', 'alpha..1.cm.', 'n', 'theta_fc', 'h_fc', 'Roseta.model'), varnames = c('clay', 'silt', 'sand', 'oc', 'bd_13b', 'theta_0.33b', 'theta_15b', 'ksat', 'f_ksat', 'theta_r', 'theta_s', 'alpha', 'n', 'theta_fc', 'h_fc', 'Rosetta.model')) { #lep is linear extensibility
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

soils_modeled_30cm <- depth_X_properties(horizon_SPC = horizons_modeled, depth = 30)
# dim(soils_modeled_30cm)

soils_modeled_10cm <- depth_X_properties(horizon_SPC = horizons_modeled, depth = 10)
dim(soils_modeled_10cm)

# soils_modeled_10_50cm <- depth_X_properties(horizon_SPC = horizons_modeled, upper_depth = 10, depth = 50)

soils_modeled_50cm <- depth_X_properties(horizon_SPC = horizons_modeled, depth = 50)

#check textural class--function taken from ssurgo_allCA_aggregation_awc_r.R
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

#condense trafficability results from HYDRUS runs
fnames_results <- list.files(summaryDir2, full.names = FALSE, recursive = FALSE)
head(fnames_results)
tail(fnames_results)
cokeys_modeled <- sub(pattern = 'soil_', replacement = '', x = fnames_results)
cokeys_modeled <- sub(pattern = '_results.csv', replacement = '', x = cokeys_modeled)
head(cokeys_modeled)
tail(cokeys_modeled) #5494 because 15 soils produced error messages
length(cokeys_modeled)

#add textural class to 10-50 cm summary
soils_modeled_10_50cm$texture.sums <- soils_modeled_10_50cm$clay_50cm + soils_modeled_10_50cm$silt_50cm + soils_modeled_10_50cm$sand_50cm
sum(is.na(soils_modeled_10_50cm$texture.sums)) #196 are NA
sum(soils_modeled_10_50cm$texture.sums > 100.1 | soils_modeled_10_50cm$texture.sums < 99.9, na.rm = TRUE) #1 problem
summary(soils_modeled_10_50cm$texture.sums)

soils_modeled_10_50cm$textural_class <- textural.class.calc(sand = soils_modeled_10_50cm$sand_50cm, silt = soils_modeled_10_50cm$silt_50cm, clay = soils_modeled_10_50cm$clay_50cm)
table(soils_modeled_10_50cm$textural_class)

#analysis with a 0-10 cm soil property summary
soils_modeled_10cm$texture.sums <- soils_modeled_10cm$clay_10cm + soils_modeled_10cm$silt_10cm + soils_modeled_10cm$sand_10cm
sum(is.na(soils_modeled_10cm$texture.sums)) #0 are NA
sum(soils_modeled_10cm$texture.sums > 100.1 | soils_modeled_10cm$texture.sums < 99.9 ) #no hay problemsas

soils_modeled_10cm$textural_class <- textural.class.calc(sand = soils_modeled_10cm$sand_10cm, silt = soils_modeled_10cm$silt_10cm, clay = soils_modeled_10cm$clay_10cm)
table(soils_modeled_10cm$textural_class)
soils_modeled_10cm[soils_modeled_10cm$cokey==2688,]
# write.csv(soils_modeled_10cm, file.path(tablesDir, 'soils_modeled_0_10cm_properties.csv'))

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
tail(cokeys_modeled) #5494 because 15 soils produced error messages
length(cokeys_modeled)
trafficability_0_10cm_v2 <- lapply(fnames_results, extract_0_10cm_trafficability_v2)
length(trafficability_0_10cm_v2)
names(trafficability_0_10cm_v2) <- cokeys_modeled

#approach to summarize data on summaryDir2
list.files(tablesDir)
FC_thresholds <- read.csv(file.path(tablesDir, 'trafficability_defs_11_3_20.csv'), stringsAsFactors = FALSE)
FC_thresholds
FC_thresholds$Option_4 <- 0.9
FC_thresholds$Option_4[FC_thresholds$Textural.Class=='clay'] <- 0.8
dim(soils_modeled_10cm) #5509
soils_modeled_10cm <- soils_modeled_10cm[soils_modeled_10cm$cokey %in% names(trafficability_0_10cm_v2),]
dim(soils_modeled_10cm) #5494
soils_modeled_10cm$traff_def_opt1 <- FC_thresholds$Option_1[match(soils_modeled_10cm$textural_class, FC_thresholds$Textural.Class)]
soils_modeled_10cm$traff_def_opt2 <- FC_thresholds$Option_2[match(soils_modeled_10cm$textural_class, FC_thresholds$Textural.Class)]
soils_modeled_10cm$traff_def_opt3 <- FC_thresholds$Option_3[match(soils_modeled_10cm$textural_class, FC_thresholds$Textural.Class)]
soils_modeled_10cm$traff_def_opt4 <- FC_thresholds$Option_4[match(soils_modeled_10cm$textural_class, FC_thresholds$Textural.Class)]
head(soils_modeled_10cm)
# names(trafficability_0_10cm_v2)[947] == soils_modeled_10cm$cokey[947]
theta_seq <- as.character(seq(0.95, 0.7, -0.01)) #how results are order in trafficability_0_10cm_v2
#exclude silt cokeys
silt_cokeys <- soils_modeled_10cm$cokey[soils_modeled_10cm$textural_class=='silt'] #because no plasticity based guide for these
trafficability_0_10cm_v2 <- trafficability_0_10cm_v2[-which(names(trafficability_0_10cm_v2) %in% silt_cokeys)]
length(trafficability_0_10cm_v2)
soils_modeled_10cm <- soils_modeled_10cm[!soils_modeled_10cm$textural_class=='silt',]
dim(soils_modeled_10cm)
match_test <- sapply(seq_along(soils_modeled_10cm$cokey), function(i) names(trafficability_0_10cm_v2)[i] == soils_modeled_10cm$cokey[i])
all(match_test) #everything in order
# test <- sapply(1:1200, function(i) {trafficability_0_10cm_v2[[i]][as.character(soils_modeled_10cm$traff_def_opt1[i])==theta_seq,]})
# test
# rm(test)

soils_modeled_10cm$result_opt1 <- sapply(seq_along(soils_modeled_10cm$traff_def_opt1), function(i) {trafficability_0_10cm_v2[[i]][as.character(soils_modeled_10cm$traff_def_opt1[i])==theta_seq,]})
summary(soils_modeled_10cm$result_opt1) #626 NA's
soils_modeled_10cm$result_opt2 <- sapply(seq_along(soils_modeled_10cm$traff_def_opt2), function(i) {trafficability_0_10cm_v2[[i]][as.character(soils_modeled_10cm$traff_def_opt2[i])==theta_seq,]})
summary(soils_modeled_10cm$result_opt2) #139 NA's
soils_modeled_10cm$result_opt3 <- sapply(seq_along(soils_modeled_10cm$traff_def_opt3), function(i) {trafficability_0_10cm_v2[[i]][as.character(soils_modeled_10cm$traff_def_opt3[i])==theta_seq,]})
summary(soils_modeled_10cm$result_opt3) #147 NA's
soils_modeled_10cm$result_opt4 <- sapply(seq_along(soils_modeled_10cm$traff_def_opt4), function(i) {trafficability_0_10cm_v2[[i]][as.character(soils_modeled_10cm$traff_def_opt4[i])==theta_seq,]})
summary(soils_modeled_10cm$result_opt4)
traffic3_by_texture_lm <- lm(result_opt3 ~ textural_class, data = soils_modeled_10cm)
summary(traffic3_by_texture_lm) #r^2=0.64

#get some names by textural class
soils_modeled_10cm[soils_modeled_10cm$textural_class=='clay' & is.na(soils_modeled_10cm$result_opt1),]
clay_cokeys_noresults <- soils_modeled_10cm$cokey[soils_modeled_10cm$textural_class=='clay' & is.na(soils_modeled_10cm$result_opt3)]
write.csv(clay_cokeys_noresults, file.path(tablesDir, 'clay_cokeys_no_results.csv')) #hay 50

#add 10-50 cm info
soils_modeled_10cm$textural_class_10_50cm <- soils_modeled_10_50cm$textural_class[match(soils_modeled_10cm$cokey, soils_modeled_10_50cm$cokey)]
table(soils_modeled_10cm$textural_class_10_50cm)
summary(lm(result_opt3 ~ textural_class + textural_class_10_50cm, data = soils_modeled_10cm))
soils_modeled_10cm$ksat_10_50cm <- soils_modeled_10_50cm$ksat_50cm[match(soils_modeled_10cm$cokey, soils_modeled_10_50cm$cokey)]

#look at oc and bd correlation
correlation_results <- soils_modeled_10cm
correlation_results <- correlation_results[!is.na(correlation_results$result_opt3), ]
dim(correlation_results) #5344
correlation_results <- correlation_results[!is.na(correlation_results$oc_10cm), ]
correlation_results <- correlation_results[!is.na(correlation_results$bd_13b_10cm), ]
dim(correlation_results)
sum(correlation_results$oc_10cm > 10) #16
sum(correlation_results$oc_10cm > 8) #26
sum(correlation_results$oc_10cm > 7) #44
sum(correlation_results$oc_10cm > 5) #191
correlation_results <- correlation_results[correlation_results$oc_10cm > 0 & correlation_results$oc_10cm <=8,]
hist(correlation_results$oc_10cm)
dim(correlation_results) #4917
summary(correlation_results$bd_13b_10cm)
hist(correlation_results$bd_13b_10cm)
sum(correlation_results$bd_13b_10cm < 0.8) #18
sum(correlation_results$bd_13b_10cm < 0.9) #49
sum(correlation_results$bd_13b_10cm < 1) #64
correlation_results <- correlation_results[!correlation_results$bd_13b_10cm < 1, ]
textural_classes <- unique(correlation_results$textural_class)
for (i in seq_along(textural_classes)) {
  print(textural_classes[i])
  print(summary(lm(correlation_results$result_opt3[correlation_results$textural_class==textural_classes[i]] ~ correlation_results$oc_10cm[correlation_results$textural_class==textural_classes[i]])))
}
for (i in seq_along(textural_classes)) {
  print(textural_classes[i])
  print(summary(lm(correlation_results$result_opt3[correlation_results$textural_class==textural_classes[i]] ~ correlation_results$bd_13b_10cm[correlation_results$textural_class==textural_classes[i]])))
}

plot(correlation_results$oc_10cm[correlation_results$textural_class=='loam'], correlation_results$result_opt3[correlation_results$textural_class=='loam'])
plot(correlation_results$oc_10cm[correlation_results$textural_class=='silt loam'], correlation_results$result_opt3[correlation_results$textural_class=='silt loam'])
plot(correlation_results$oc_10cm[correlation_results$textural_class=='silty clay loam'], correlation_results$result_opt3[correlation_results$textural_class=='silty clay loam'])
plot(correlation_results$oc_10cm[correlation_results$textural_class=='silty clay loam'], correlation_results$result_opt3[correlation_results$textural_class=='silty clay loam'])
traffic_by_texture_lm <- lm(result_opt3 ~ textural_class, data = correlation_results)
summary(traffic_by_texture_lm)
data.frame(texture=textural_classes, predict(traffic_by_texture_lm, data.frame(textural_class=textural_classes), se.fit = TRUE, interval = 'prediction', level = 0.5))
predict(traffic_by_texture_lm, data.frame(textural_class='clay'), se.fit = TRUE, interval = 'prediction', level = 0.9)
quantile(correlation_results$result_opt3[correlation_results$textural_class=='clay'], probs=c(0.975, 0.95, 0.9, 0.75, 0.5, 0.25, 0.1, 0.05, 0.025))
hist(correlation_results$result_opt3[correlation_results$textural_class=='clay'])
predict(traffic_by_texture_lm, data.frame(textural_class='sand'), se.fit = TRUE, interval = 'prediction', level = 0.85)
traffic_by_texture_om_bd_lm <- lm(result_opt3 ~ textural_class + oc_10cm + bd_13b_10cm, data = correlation_results)
summary(traffic_by_texture_om_bd_lm)

traffic_by_clay_lm <- lm(result_opt3 ~ clay_10cm, data = correlation_results)
summary(traffic_by_clay_lm) #0.60
data.frame(clay=100*seq(0, 0.95, by=0.05), predict(traffic_by_clay_lm, data.frame(clay_10cm=100*seq(0, 0.95, by=0.05)), se.fit = TRUE, interval = 'prediction', level = 0.9)) #negative predictions for clay < 15%

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
  write.csv(soilname_summary, file.path(tablesDir, paste0('summary_by_soilname_0_10cm_properties_', traff_def, '.csv')), row.names = TRUE)
}
summarize_by_names('result_opt1')
summarize_by_names('result_opt2')
summarize_by_names('result_opt3')
summarize_by_names('result_opt4')

summarize_by_texture <- function(traff_def) {
  textural_class_summary <- do.call(cbind, lapply(unique(soils_modeled_10cm$textural_class), function(x) as.matrix(summary(soils_modeled_10cm[[traff_def]][soils_modeled_10cm$textural_class==x]))[1:6,]))
  colnames(textural_class_summary) <- unique(soils_modeled_10cm$textural_class)
  textural_class_summary <- rbind(textural_class_summary, n=as.integer(sapply(unique(soils_modeled_10cm$textural_class), function(x) sum(!is.na(soils_modeled_10cm[[traff_def]][soils_modeled_10cm$textural_class==x])))))
  textural_class_summary <- rbind(textural_class_summary, mean_clay=sapply(unique(soils_modeled_10cm$textural_class), function(x) mean(soils_modeled_10cm$clay_10cm[soils_modeled_10cm$textural_class==x], na.rm=TRUE)))
  textural_class_summary <- rbind(textural_class_summary, Not_Determined=as.integer(sapply(unique(soils_modeled_10cm$textural_class), function(x) sum(is.na(soils_modeled_10cm[[traff_def]][soils_modeled_10cm$textural_class==x])))))
  colnames(textural_class_summary)[order(textural_class_summary[8,])]
textural_class_summary <- textural_class_summary[,order(textural_class_summary[8,])]
  print(textural_class_summary)
  write.csv(textural_class_summary, file.path(tablesDir, paste0('summary_by_textural_class_0_10cm_', traff_def, '.csv')), row.names = TRUE)
}
summarize_by_texture('result_opt1')
summarize_by_texture('result_opt2')
summarize_by_texture('result_opt3')
summarize_by_texture('result_opt4')

#look at rosetta model=5 only
table(soils_modeled_10cm$Rosetta.model_10cm)
sum(soils_modeled_10cm$Rosetta.model_10cm==5) #5130
soils_modeled_10cm_rosetta5 <- soils_modeled_10cm[soils_modeled_10cm$Rosetta.model_10cm==5,]
dim(soils_modeled_10cm_rosetta5) #5130
summary(lm(result_opt3 ~ textural_class, data = soils_modeled_10cm_rosetta5)) #r2=0.6981, rse=5.77, 134 are NA
summary(lm(result_opt3 ~ textural_class + bd_13b_10cm, data = soils_modeled_10cm_rosetta5)) #r2=0.6998, rse=5.758
summary(lm(result_opt3 ~ textural_class + bd_13b_10cm + oc_10cm, data = soils_modeled_10cm_rosetta5)) #r2=0.7095, rse=0.71
summary(lm(result_opt3 ~ textural_class + textural_class_10_50cm, data = soils_modeled_10cm_rosetta5)) #r2=0.7059, rse=5.752, 310 are NA
summary(lm(result_opt3 ~ textural_class + textural_class_10_50cm + bd_13b_10cm + oc_10cm, data = soils_modeled_10cm_rosetta5)) #r2=0.7191, rse=5.634, 362 are NA


summary(lm(result_opt2 ~ textural_class, data = soils_modeled_10cm_rosetta5)) #r2=0.6675, rse=5.289, 126 are NA
summary(lm(result_opt2 ~ textural_class + bd_13b_10cm, data = soils_modeled_10cm_rosetta5)) #r2=0.6682, rse=5.284
summary(lm(result_opt2 ~ textural_class + bd_13b_10cm + oc_10cm, data = soils_modeled_10cm_rosetta5)) #r2=0.682, rse=5.164
summary(lm(result_opt2 ~ textural_class + textural_class_10_50cm, data = soils_modeled_10cm_rosetta5)) #r2=0.6759, rse=5.288, 302 are NA
summary(lm(result_opt2 ~ clay_10cm + sand_10cm + bd_13b_10cm, data = soils_modeled_10cm_rosetta5)) #r2=0.5271, rse=6.292
summary(lm(result_opt2 ~ ksat_10cm, data = soils_modeled_10cm_rosetta5)) #r2=0.05252, rse=8.921
summary(lm(result_opt2 ~ ksat_10cm + ksat_10_50cm, data = soils_modeled_10cm_rosetta5)) #r2=0.06094, rse=8.985
summary(lm(result_opt2 ~ textural_class + textural_class_10_50cm + bd_13b_10cm + oc_10cm, data = soils_modeled_10cm_rosetta5)) #r2=0.6922, rse=5.145, 355 are NA
summary(lm(result_opt2 ~ ksat_10cm + textural_class, data=soils_modeled_10cm_rosetta5)) #r2=0.668, rse=5.286

lapply(1:length(textural_classes), function(i) {
  hist(soils_modeled_10cm_rosetta5$result_opt3[soils_modeled_10cm_rosetta5$textural_class==textural_classes[i]], main=textural_classes[i])
})
lapply(1:length(textural_classes), function(i) {
  hist(soils_modeled_10cm_rosetta5$result_opt1[soils_modeled_10cm_rosetta5$textural_class==textural_classes[i]], main=textural_classes[i])
})
names(textural_classes) <- textural_classes
quantile_breaks_result3 <- lapply(textural_classes, function(x) {
  quantile(soils_modeled_10cm_rosetta5$result_opt3[soils_modeled_10cm_rosetta5$textural_class==x], probs = c(0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99), na.rm = TRUE)
})
quantile_breaks_result3 <- rbind(as.data.frame(quantile_breaks_result3), sd=sapply(textural_classes, function(x) {sd(soils_modeled_10cm_rosetta5$result_opt3[soils_modeled_10cm_rosetta5$textural_class==x], na.rm=TRUE)}), mean=sapply(textural_classes, function(x) {mean(soils_modeled_10cm_rosetta5$result_opt3[soils_modeled_10cm_rosetta5$textural_class==x], na.rm=TRUE)}))
write.csv(quantile_breaks_result3, file.path(tablesDir, 'quantiles_result3_rosetta5.csv'), row.names = TRUE)

quantile_breaks_result2 <- lapply(textural_classes, function(x) {
  quantile(soils_modeled_10cm_rosetta5$result_opt2[soils_modeled_10cm_rosetta5$textural_class==x], probs = c(0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99), na.rm = TRUE)
})
quantile_breaks_result2 <- rbind(as.data.frame(quantile_breaks_result2), sd=sapply(textural_classes, function(x) {sd(soils_modeled_10cm_rosetta5$result_opt2[soils_modeled_10cm_rosetta5$textural_class==x], na.rm=TRUE)}), mean=sapply(textural_classes, function(x) {mean(soils_modeled_10cm_rosetta5$result_opt2[soils_modeled_10cm_rosetta5$textural_class==x], na.rm=TRUE)}))
write.csv(quantile_breaks_result2, file.path(tablesDir, 'quantiles_result2_rosetta5.csv'), row.names = TRUE)
for (i in seq_along(textural_classes)) {
  if (length(unique(soils_modeled_10cm_rosetta5$textural_class_10_50cm[soils_modeled_10cm_rosetta5$textural_class==textural_classes[i]])) >= 2) {
    print(textural_classes[i])
    print(summary(lm(soils_modeled_10cm_rosetta5$result_opt2[soils_modeled_10cm_rosetta5$textural_class==textural_classes[i]] ~ soils_modeled_10cm_rosetta5$textural_class_10_50cm[soils_modeled_10cm_rosetta5$textural_class==textural_classes[i]])))
  } else {next}
}
#now look at BD as factor
for (i in seq_along(textural_classes)) {
  print(textural_classes[i])
  print(summary(lm(soils_modeled_10cm_rosetta5$result_opt2[soils_modeled_10cm_rosetta5$textural_class==textural_classes[i]] ~ soils_modeled_10cm_rosetta5$bd_13b_10cm[soils_modeled_10cm_rosetta5$textural_class==textural_classes[i]])))
}
#now look at Ks as factor
for (i in seq_along(textural_classes)) {
  print(textural_classes[i])
  print(summary(lm(soils_modeled_10cm_rosetta5$result_opt2[soils_modeled_10cm_rosetta5$textural_class==textural_classes[i]] ~ soils_modeled_10cm_rosetta5$ksat_10cm[soils_modeled_10cm_rosetta5$textural_class==textural_classes[i]])))
}
#now look at h_fc assumption as factor
for (i in seq_along(textural_classes)) {
  print(textural_classes[i])
  print(summary(lm(soils_modeled_10cm_rosetta5$result_opt2[soils_modeled_10cm_rosetta5$textural_class==textural_classes[i]] ~ soils_modeled_10cm_rosetta5$h_fc_10cm[soils_modeled_10cm_rosetta5$textural_class==textural_classes[i]])))
}
#now look at h_fc, bd, and ksat as explaining variation within each 0-10 cm texture class
for (i in seq_along(textural_classes)) {
  print(textural_classes[i])
  print(summary(lm(soils_modeled_10cm_rosetta5$result_opt2[soils_modeled_10cm_rosetta5$textural_class==textural_classes[i]] ~ soils_modeled_10cm_rosetta5$h_fc_10cm[soils_modeled_10cm_rosetta5$textural_class==textural_classes[i]] + soils_modeled_10cm_rosetta5$ksat_10cm[soils_modeled_10cm_rosetta5$textural_class==textural_classes[i]] + soils_modeled_10cm_rosetta5$bd_13b_10cm[soils_modeled_10cm_rosetta5$textural_class==textural_classes[i]])))
}

#now look at h_fc, bd, and ksat as explaining variation within each 0-10 cm texture class
for (i in seq_along(textural_classes)) {
  if (length(unique(soils_modeled_10cm_rosetta5$textural_class_10_50cm[soils_modeled_10cm_rosetta5$textural_class==textural_classes[i]])) >= 2) {
    print(textural_classes[i])
    print(summary(lm(soils_modeled_10cm_rosetta5$result_opt2[soils_modeled_10cm_rosetta5$textural_class==textural_classes[i]] ~ soils_modeled_10cm_rosetta5$h_fc_10cm[soils_modeled_10cm_rosetta5$textural_class==textural_classes[i]] + soils_modeled_10cm_rosetta5$ksat_10cm[soils_modeled_10cm_rosetta5$textural_class==textural_classes[i]] + soils_modeled_10cm_rosetta5$bd_13b_10cm[soils_modeled_10cm_rosetta5$textural_class==textural_classes[i]] +  soils_modeled_10cm_rosetta5$textural_class_10_50cm[soils_modeled_10cm_rosetta5$textural_class==textural_classes[i]])))} else {next}
}

#results as function of 0-10 cm ksat within each textural class 
plot(soils_modeled_10cm_rosetta5$ksat_10cm[soils_modeled_10cm_rosetta5$textural_class=='clay'], soils_modeled_10cm_rosetta5$result_opt2[soils_modeled_10cm_rosetta5$textural_class=='clay'])
plot(soils_modeled_10cm_rosetta5$ksat_10cm[soils_modeled_10cm_rosetta5$textural_class=='sand'], soils_modeled_10cm_rosetta5$result_opt2[soils_modeled_10cm_rosetta5$textural_class=='sand'])
plot(soils_modeled_10cm_rosetta5$ksat_10cm[soils_modeled_10cm_rosetta5$textural_class=='sandy loam'], soils_modeled_10cm_rosetta5$result_opt2[soils_modeled_10cm_rosetta5$textural_class=='sandy loam'])
plot(soils_modeled_10cm_rosetta5$ksat_10cm[soils_modeled_10cm_rosetta5$textural_class=='loam'], soils_modeled_10cm_rosetta5$result_opt2[soils_modeled_10cm_rosetta5$textural_class=='loam'])
plot(soils_modeled_10cm_rosetta5$ksat_10cm[soils_modeled_10cm_rosetta5$textural_class=='loamy sand'], soils_modeled_10cm_rosetta5$result_opt2[soils_modeled_10cm_rosetta5$textural_class=='loamy sand'])
plot(soils_modeled_10cm_rosetta5$ksat_10cm[soils_modeled_10cm_rosetta5$textural_class=='loamy sand'], soils_modeled_10cm_rosetta5$result_opt2[soils_modeled_10cm_rosetta5$textural_class=='loamy sand'])

#boxplots
for (i in seq_along(textural_classes)) {
  test <- boxplot(soils_modeled_10cm_rosetta5$result_opt2[soils_modeled_10cm_rosetta5$textural_class==textural_classes[i]], main=textural_classes[i], plot = TRUE, range = 3)
  print(textural_classes[i])
  print(test)
}

#analysis to do a summary by 0-30 cm textural classes
#analysis with a 0-10 cm soil property summary
soils_modeled_30cm$textural_class <- textural.class.calc(sand = soils_modeled_30cm$sand_30cm, silt = soils_modeled_30cm$silt_30cm, clay = soils_modeled_30cm$clay_30cm)
table(soils_modeled_30cm$textural_class)
soils_modeled_30cm$texture.sums <- soils_modeled_30cm$clay_30cm + soils_modeled_30cm$silt_30cm + soils_modeled_30cm$sand_30cm
sum(is.na(soils_modeled_30cm$texture.sums)) #45 are NA
sum(soils_modeled_30cm$texture.sums > 100.1 | soils_modeled_30cm$texture.sums < 99.9, na.rm = TRUE) #no hay problemsas

# write.csv(soils_modeled_30cm, file.path(tablesDir, 'soils_modeled_0_30cm_properties.csv'))

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
length(cokeys_modeled) #5494 because 15 soils produced error messages
trafficability_0_10cm_v2 <- lapply(fnames_results, extract_0_10cm_trafficability_v2)
length(trafficability_0_10cm_v2)
names(trafficability_0_10cm_v2) <- cokeys_modeled

#approach to summarize data on summaryDir2
dim(soils_modeled_30cm) #5509
soils_modeled_30cm <- soils_modeled_30cm[soils_modeled_30cm$cokey %in% names(trafficability_0_10cm_v2),]
dim(soils_modeled_30cm) #5447
soils_modeled_30cm$traff_def_opt1 <- FC_thresholds$Option_1[match(soils_modeled_30cm$textural_class, FC_thresholds$Textural.Class)]
soils_modeled_30cm$traff_def_opt2 <- FC_thresholds$Option_2[match(soils_modeled_30cm$textural_class, FC_thresholds$Textural.Class)]
soils_modeled_30cm$traff_def_opt3 <- FC_thresholds$Option_3[match(soils_modeled_30cm$textural_class, FC_thresholds$Textural.Class)]
soils_modeled_30cm$traff_def_opt4 <- FC_thresholds$Option_4[match(soils_modeled_30cm$textural_class, FC_thresholds$Textural.Class)]
head(soils_modeled_30cm)
theta_seq <- as.character(seq(0.95, 0.7, -0.01)) #how results are order in trafficability_0_10cm_v2

#exclude 
no_texture_cokeys <- soils_modeled_30cm$cokey[is.na(soils_modeled_30cm$textural_class)]
no_texture_cokeys
trafficability_0_10cm_v2 <- trafficability_0_10cm_v2[-which(names(trafficability_0_10cm_v2) %in% no_texture_cokeys)]
soils_modeled_30cm <- soils_modeled_30cm[!is.na(soils_modeled_30cm$textural_class),]

#exclude silt cokeys
silt_cokeys <- soils_modeled_30cm$cokey[soils_modeled_30cm$textural_class=='silt'] #because no plasticity based guide for these
trafficability_0_10cm_v2 <- trafficability_0_10cm_v2[-which(names(trafficability_0_10cm_v2) %in% silt_cokeys)]
length(trafficability_0_10cm_v2) #5447
soils_modeled_30cm <- soils_modeled_30cm[-which(soils_modeled_30cm$textural_class=='silt'),] #use which because there are NA textural class
dim(soils_modeled_30cm) #5447
match_test <- sapply(seq_along(soils_modeled_30cm$cokey), function(i) names(trafficability_0_10cm_v2)[i] == soils_modeled_30cm$cokey[i])
all(match_test) #everything in order
# test <- sapply(1:1200, function(i) {trafficability_0_10cm_v2[[i]][as.character(soils_modeled_30cm$traff_def_opt1[i])==theta_seq,]})
# test
# rm(test)

soils_modeled_30cm$result_opt1 <- sapply(seq_along(soils_modeled_30cm$traff_def_opt1), function(i) {trafficability_0_10cm_v2[[i]][as.character(soils_modeled_30cm$traff_def_opt1[i])==theta_seq,]})
summary(soils_modeled_30cm$result_opt1) #640 NA's
soils_modeled_30cm$result_opt2 <- sapply(seq_along(soils_modeled_30cm$traff_def_opt2), function(i) {trafficability_0_10cm_v2[[i]][as.character(soils_modeled_30cm$traff_def_opt2[i])==theta_seq,]})
summary(soils_modeled_30cm$result_opt2) #140 NA's
soils_modeled_30cm$result_opt3 <- sapply(seq_along(soils_modeled_30cm$traff_def_opt3), function(i) {trafficability_0_10cm_v2[[i]][as.character(soils_modeled_30cm$traff_def_opt3[i])==theta_seq,]})
summary(soils_modeled_30cm$result_opt3) #148 NA's
soils_modeled_30cm$result_opt4 <- sapply(seq_along(soils_modeled_30cm$traff_def_opt4), function(i) {trafficability_0_10cm_v2[[i]][as.character(soils_modeled_30cm$traff_def_opt4[i])==theta_seq,]})
summary(soils_modeled_30cm$result_opt4) #141 NAs

#limit to rosetta 5 profiles
soils_modeled_30cm <- soils_modeled_30cm[soils_modeled_30cm$Rosetta.model_30cm==5,]
dim(soils_modeled_30cm) #4979
5447-4979 #another 468 dropped

traffic2_by_texture_lm <- lm(result_opt2 ~ textural_class, data = soils_modeled_30cm)
summary(traffic2_by_texture_lm) #r2=0.6454, RSE=5.567
traffic3_by_texture_lm <- lm(result_opt3 ~ textural_class, data = soils_modeled_30cm)
summary(traffic3_by_texture_lm) #r^2=0.6719, RSE=6.07

summary(lm(result_opt2 ~ textural_class + bd_13b_30cm, data = soils_modeled_30cm)) #r2=0.6458, rse=5.564
summary(lm(result_opt2 ~ textural_class + bd_13b_30cm + oc_30cm, data = soils_modeled_30cm)) #r2=0.6614, rse=5.434
# summary(lm(result_opt2 ~ textural_class + textural_class_10_50cm, data = soils_modeled_30cm)) #r2=0.6759, rse=5.288, 302 are NA
summary(lm(result_opt2 ~ clay_30cm + sand_30cm + bd_13b_30cm, data = soils_modeled_30cm)) #r2=0.4944, rse=6.642
summary(lm(result_opt2 ~ clay_30cm + sand_30cm + bd_13b_30cm + theta_0.33b_30cm + theta_15b_30cm, data = soils_modeled_30cm))
summary(lm(result_opt2 ~ ksat_30cm, data = soils_modeled_30cm)) #r2=0.06, rse=9.051
# summary(lm(result_opt2 ~ ksat_30cm + ksat_10_50cm, data = soils_modeled_30cm)) #r2=0.06094, rse=8.985
# summary(lm(result_opt2 ~ textural_class + textural_class_10_50cm + bd_13b_30cm + oc_30cm, data = soils_modeled_30cm)) #r2=0.6922, rse=5.145, 355 are NA
summary(lm(result_opt2 ~ ksat_30cm + textural_class, data=soils_modeled_30cm)) #r2=0.6463, rse=5.56


summarize_by_texture <- function(traff_def) {
  textural_class_summary <- do.call(cbind, lapply(unique(soils_modeled_30cm$textural_class), function(x) as.matrix(summary(soils_modeled_30cm[[traff_def]][soils_modeled_30cm$textural_class==x]))[1:6,]))
  colnames(textural_class_summary) <- unique(soils_modeled_30cm$textural_class)
  textural_class_summary <- rbind(textural_class_summary, n=as.integer(sapply(unique(soils_modeled_30cm$textural_class), function(x) sum(!is.na(soils_modeled_30cm[[traff_def]][soils_modeled_30cm$textural_class==x])))))
  textural_class_summary <- rbind(textural_class_summary, mean_clay=sapply(unique(soils_modeled_30cm$textural_class), function(x) mean(soils_modeled_30cm$clay_30cm[soils_modeled_30cm$textural_class==x], na.rm=TRUE)))
  textural_class_summary <- rbind(textural_class_summary, Not_Determined=as.integer(sapply(unique(soils_modeled_30cm$textural_class), function(x) sum(is.na(soils_modeled_30cm[[traff_def]][soils_modeled_30cm$textural_class==x])))))
  colnames(textural_class_summary)[order(textural_class_summary[8,])]
  textural_class_summary <- textural_class_summary[,order(textural_class_summary[8,])]
  print(textural_class_summary)
  write.csv(textural_class_summary, file.path(tablesDir, paste0('summary_by_textural_class_0_30cm_', traff_def, '.csv')), row.names = TRUE)
}
summarize_by_texture('result_opt1')
summarize_by_texture('result_opt2')
summarize_by_texture('result_opt3')
summarize_by_texture('result_opt4')

textural_classes <- unique(soils_modeled_30cm$textural_class)
names(textural_classes) <- textural_classes
quantile_breaks_result2 <- lapply(textural_classes, function(x) {
  quantile(soils_modeled_30cm$result_opt2[soils_modeled_30cm$textural_class==x], probs = c(0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99), na.rm = TRUE)
})
quantile_breaks_result2 <- rbind(as.data.frame(quantile_breaks_result2), sd=sapply(textural_classes, function(x) {sd(soils_modeled_30cm$result_opt2[soils_modeled_30cm$textural_class==x], na.rm=TRUE)}), mean=sapply(textural_classes, function(x) {mean(soils_modeled_30cm$result_opt2[soils_modeled_30cm$textural_class==x], na.rm=TRUE)}))
quantile_breaks_result2[6,]
quantile_breaks_result2 <- quantile_breaks_result2[,order(quantile_breaks_result2[6,])]
write.csv(quantile_breaks_result2, file.path(tablesDir, 'quantiles_result2_0_30cm_properties_rosetta5.csv'), row.names = TRUE)

#now look at BD as factor
for (i in seq_along(textural_classes)) {
  plot(soils_modeled_30cm$bd_13b_30cm[soils_modeled_30cm$textural_class==textural_classes[i]], soils_modeled_30cm$result_opt2[soils_modeled_30cm$textural_class==textural_classes[i]], main=textural_classes[i])
  print(textural_classes[i])
  print(summary(lm(soils_modeled_30cm$result_opt2[soils_modeled_30cm$textural_class==textural_classes[i]] ~ soils_modeled_30cm$bd_13b_30cm[soils_modeled_30cm$textural_class==textural_classes[i]])))
}
#now look at Ks as factor
for (i in seq_along(textural_classes)) {
  plot(soils_modeled_30cm$ksat_30cm[soils_modeled_30cm$textural_class==textural_classes[i]], soils_modeled_30cm$result_opt2[soils_modeled_30cm$textural_class==textural_classes[i]], main=textural_classes[i])
  print(textural_classes[i])
  print(summary(lm(soils_modeled_30cm$result_opt2[soils_modeled_30cm$textural_class==textural_classes[i]] ~ soils_modeled_30cm$ksat_30cm[soils_modeled_30cm$textural_class==textural_classes[i]])))
}
#now look at h_fc assumption as factor
for (i in seq_along(textural_classes)) {
  plot(soils_modeled_30cm$h_fc_30cm[soils_modeled_30cm$textural_class==textural_classes[i]], soils_modeled_30cm$result_opt2[soils_modeled_30cm$textural_class==textural_classes[i]], main=textural_classes[i], xlim=c(-500,-50))
  print(textural_classes[i])
  print(summary(lm(soils_modeled_30cm$result_opt2[soils_modeled_30cm$textural_class==textural_classes[i]] ~ soils_modeled_30cm$h_fc_30cm[soils_modeled_30cm$textural_class==textural_classes[i]])))
}
#now look at %clay within each textural class
for (i in seq_along(textural_classes)) {
  plot(soils_modeled_30cm$clay_30cm[soils_modeled_30cm$textural_class==textural_classes[i]], soils_modeled_30cm$result_opt2[soils_modeled_30cm$textural_class==textural_classes[i]], main=textural_classes[i])
  print(textural_classes[i])
  print(summary(lm(soils_modeled_30cm$result_opt2[soils_modeled_30cm$textural_class==textural_classes[i]] ~ soils_modeled_30cm$clay_30cm[soils_modeled_30cm$textural_class==textural_classes[i]])))
}
#now look at wilting point within each textural class
for (i in seq_along(textural_classes)) {
  plot(soils_modeled_30cm$theta_15b_30cm[soils_modeled_30cm$textural_class==textural_classes[i]], soils_modeled_30cm$result_opt2[soils_modeled_30cm$textural_class==textural_classes[i]], main=textural_classes[i])
  print(textural_classes[i])
  print(summary(lm(soils_modeled_30cm$result_opt2[soils_modeled_30cm$textural_class==textural_classes[i]] ~ soils_modeled_30cm$theta_15b_30cm[soils_modeled_30cm$textural_class==textural_classes[i]])))
}
#now look at theta at 0.33b within each textural class
for (i in seq_along(textural_classes)) {
  plot(soils_modeled_30cm$theta_0.33b_30cm[soils_modeled_30cm$textural_class==textural_classes[i]], soils_modeled_30cm$result_opt2[soils_modeled_30cm$textural_class==textural_classes[i]], main=textural_classes[i])
  print(textural_classes[i])
  print(summary(lm(soils_modeled_30cm$result_opt2[soils_modeled_30cm$textural_class==textural_classes[i]] ~ soils_modeled_30cm$theta_0.33b_30cm[soils_modeled_30cm$textural_class==textural_classes[i]])))
}
#now look at theta at 0.33b within each textural class
for (i in seq_along(textural_classes)) {
  plot(soils_modeled_30cm$theta_fc_30cm[soils_modeled_30cm$textural_class==textural_classes[i]], soils_modeled_30cm$result_opt2[soils_modeled_30cm$textural_class==textural_classes[i]], main=textural_classes[i])
  print(textural_classes[i])
  print(summary(lm(soils_modeled_30cm$result_opt2[soils_modeled_30cm$textural_class==textural_classes[i]] ~ soils_modeled_30cm$theta_fc_30cm[soils_modeled_30cm$textural_class==textural_classes[i]])))
}
#now look at silt within each textural clas
for (i in seq_along(textural_classes)) {
  plot(soils_modeled_30cm$silt_30cm[soils_modeled_30cm$textural_class==textural_classes[i]], soils_modeled_30cm$result_opt2[soils_modeled_30cm$textural_class==textural_classes[i]], main=textural_classes[i])
  print(textural_classes[i])
  print(summary(lm(soils_modeled_30cm$result_opt2[soils_modeled_30cm$textural_class==textural_classes[i]] ~ soils_modeled_30cm$silt_30cm[soils_modeled_30cm$textural_class==textural_classes[i]])))
}
#now look at sand within each class
for (i in seq_along(textural_classes)) {
  plot(soils_modeled_30cm$sand_30cm[soils_modeled_30cm$textural_class==textural_classes[i]], soils_modeled_30cm$result_opt2[soils_modeled_30cm$textural_class==textural_classes[i]], main=textural_classes[i])
  print(textural_classes[i])
  print(summary(lm(soils_modeled_30cm$result_opt2[soils_modeled_30cm$textural_class==textural_classes[i]] ~ soils_modeled_30cm$sand_30cm[soils_modeled_30cm$textural_class==textural_classes[i]])))
}
plot(soils_modeled_30cm$clay_30cm, soils_modeled_30cm$result_opt2)
