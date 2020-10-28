library(aqp)
workDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/data from Stathis'
summaryDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/Oct2020test/summary'
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
depths(horizons_modeled) <- cokey ~ hzn_top + hzn_bot

wtd.mean_v2 <- function(x, y) {
  # use horizon thickness as a weight
  thick <- x$hzn_bot - x$hzn_top
  # compute the weighted mean, accounting for the possibility of missing data
  m <- weighted.mean(horizons(x)[[y]], w=thick, na.rm=FALSE)
  m
}

depth_X_properties <- function(horizon_SPC, depth, vars_of_interest = c('clay', 'silt', 'sand', 'estimated_oc', 'db_13b', 'wr_13b', 'wr_15b', 'Ks..cm.d.', 'K0.cm.d.', 'teta_r', 'teta_s', 'alpha..1.cm.', 'n'), varnames = c('clay', 'silt', 'sand', 'oc', 'bd_13b', 'theta_0.33b', 'theta_15b', 'ksat', 'f_ksat', 'theta_r', 'theta_s', 'alpha', 'n')) { #lep is linear extensibility
  columnames <- paste0(varnames, '_', depth, 'cm')
  print(cbind(vars_of_interest, columnames)) #show that it's all lined up
  assign("depth", depth, envir = .GlobalEnv) #this necessary because slice can't find the variable otherwise
  sliced_SPC <- slice(horizon_SPC, 0:(depth-1) ~ .) #depth was '0:depth' in previous version
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


soils_modeled_30cm <- depth_X_properties(horizon_SPC = horizons_modeled, depth = 30)

dim(soils_modeled_30cm)

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
soils_modeled_30cm$texture.sums <- soils_modeled_30cm$clay_30cm + soils_modeled_30cm$silt_30cm + soils_modeled_30cm$sand_30cm
sum(is.na(soils_modeled_30cm$texture.sums)) #45 are NA
incomplete_texture <- soils_modeled_30cm[is.na(soils_modeled_30cm$texture.sums),]
soils_modeled_30cm <- soils_modeled_30cm[!is.na(soils_modeled_30cm$texture.sums),]
dim(soils_modeled_30cm) #5464 rows

sum(soils_modeled_30cm$texture.sums > 100.1 | soils_modeled_30cm$texture.sums < 99.9 ) #no hay problemsas

soils_modeled_30cm$textural_class <- textural.class.calc(sand = soils_modeled_30cm$sand_30cm, silt = soils_modeled_30cm$silt_30cm, clay = soils_modeled_30cm$clay_30cm)
table(soils_modeled_30cm$textural_class)
soils_modeled_30cm[soils_modeled_30cm$cokey==2688,]
write.csv(soils_modeled_30cm, file.path(tablesDir, 'soils_modeled_0_30cm_properties.csv'))

#condense trafficability results from HYDRUS runs
fnames_results <- list.files(summaryDir, full.names = FALSE, recursive = FALSE)
head(fnames_results)
tail(fnames_results)
cokeys_modeled <- sub(pattern = 'soil_', replacement = '', x = fnames_results)
cokeys_modeled <- sub(pattern = '_results.csv', replacement = '', x = cokeys_modeled)
head(cokeys_modeled)
tail(cokeys_modeled) #5494 because 15 soils produced error messages
length(cokeys_modeled)
extract_0_10cm_trafficability <- function(x) {
  results <- read.csv(file.path(summaryDir, x), row.names=1)
  results[3,]
}
trafficability_0_10cm <- do.call(rbind, lapply(fnames_results, extract_0_10cm_trafficability))
dim(trafficability_0_10cm)
head(trafficability_0_10cm)
row.names(trafficability_0_10cm) <- cokeys_modeled
lapply(trafficability_0_10cm, function(x) sum(is.na(x)))
write.csv(trafficability_0_10cm, file.path(tablesDir, 'trafficability_0_10cm_allSoils.csv'), row.names = TRUE)
trafficability_0_10cm$textural_class_30cm <- soils_modeled_30cm$textural_class[match(row.names(trafficability_0_10cm), soils_modeled_30cm$cokey)]
trafficability_0_10cm$ksat_30cm <- soils_modeled_30cm$ksat_30cm[match(row.names(trafficability_0_10cm), soils_modeled_30cm$cokey)]
trafficability_0_10cm$alpha_30cm <- soils_modeled_30cm$alpha_30cm[match(row.names(trafficability_0_10cm), soils_modeled_30cm$cokey)]
trafficability_0_10cm$theta_r_30cm <- soils_modeled_30cm$theta_r_30cm[match(row.names(trafficability_0_10cm), soils_modeled_30cm$cokey)]
trafficability_0_10cm$theta_s_30cm <- soils_modeled_30cm$theta_s_30cm[match(row.names(trafficability_0_10cm), soils_modeled_30cm$cokey)]
trafficability_0_10cm$clay_30cm <- soils_modeled_30cm$clay_30cm[match(row.names(trafficability_0_10cm), soils_modeled_30cm$cokey)]
trafficability_0_10cm$silt_30cm <- soils_modeled_30cm$silt_30cm[match(row.names(trafficability_0_10cm), soils_modeled_30cm$cokey)]
trafficability_0_10cm$sand_30cm <- soils_modeled_30cm$sand_30cm[match(row.names(trafficability_0_10cm), soils_modeled_30cm$cokey)]
trafficability_0_10cm$oc_30cm <- soils_modeled_30cm$oc_30cm[match(row.names(trafficability_0_10cm), soils_modeled_30cm$cokey)]
trafficability_0_10cm$bd_30cm <- soils_modeled_30cm$bd_13b_30cm[match(row.names(trafficability_0_10cm), soils_modeled_30cm$cokey)]
trafficability_0_10cm$theta_0.33b_30cm <- soils_modeled_30cm$theta_0.33b_30cm[match(row.names(trafficability_0_10cm), soils_modeled_30cm$cokey)]
trafficability_0_10cm$theta_15b_30cm <- soils_modeled_30cm$theta_15b_30cm[match(row.names(trafficability_0_10cm), soils_modeled_30cm$cokey)]
trafficability_0_10cm$rosetta.model <- mod_database$Roseta.model[match(row.names(trafficability_0_10cm), mod_database$cokey)]
trafficability_0_10cm$soil_name <- mod_database$taxonname[match(row.names(trafficability_0_10cm), mod_database$cokey)]


trafficability_0_10cm[,1:4] <- lapply(trafficability_0_10cm[,1:4], function(x) as.numeric(x))
tapply(trafficability_0_10cm$FC_0.9, trafficability_0_10cm$textural_class_30cm, mean, na.rm=TRUE)
tapply(trafficability_0_10cm$FC_0.9, trafficability_0_10cm$textural_class_30cm, summary)
tapply(trafficability_0_10cm$FC_0.95, trafficability_0_10cm$textural_class_30cm, mean, na.rm=TRUE)
tapply(trafficability_0_10cm$FC_0.9, trafficability_0_10cm$rosetta.model, mean, na.rm=TRUE) #rosetta model 3 is biased
table(trafficability_0_10cm$rosetta.model)
#2    3    5 
#344   10 5140 
sum(!is.na(trafficability_0_10cm$FC_0.9) & trafficability_0_10cm$rosetta.model==5 & !is.na(trafficability_0_10cm$textural_class_30cm)) #4938

#only consider rosetta model 5
trafficability_0_10cm <- trafficability_0_10cm[trafficability_0_10cm$rosetta.model==5,]
dim(trafficability_0_10cm) #5140
soil_name_counts <- table(trafficability_0_10cm$soil_name)[order(table(trafficability_0_10cm$soil_name), decreasing=TRUE)]
soil_name_mas30 <- soil_name_counts[soil_name_counts >= 30]

soilname_summary <- do.call(cbind, lapply(names(soil_name_mas30), function(x) as.matrix(summary(trafficability_0_10cm$FC_0.9[trafficability_0_10cm$soil_name==x]))[1:6,]))
colnames(soilname_summary) <- names(soil_name_mas30) 
soilname_summary <- rbind(soilname_summary, n=as.integer(sapply(names(soil_name_mas30), function(x) sum(!is.na(trafficability_0_10cm$FC_0.9[trafficability_0_10cm$soil_name==x])))))
soilname_summary
soilname_summary <- rbind(soilname_summary, mean_clay=sapply(names(soil_name_mas30), function(x) mean(trafficability_0_10cm$clay_30cm[trafficability_0_10cm$soil_name==x], na.rm=TRUE)))
write.csv(soilname_summary, file.path(tablesDir, 'summary_by_soilname.csv'), row.names = TRUE)
mod_database[mod_database$cokey==10610083,]

textural_class_summary <- do.call(cbind, lapply(unique(trafficability_0_10cm$textural_class_30cm), function(x) as.matrix(summary(trafficability_0_10cm$FC_0.9[trafficability_0_10cm$textural_class_30cm==x]))[1:6,]))
colnames(textural_class_summary) <- unique(trafficability_0_10cm$textural_class_30cm)
textural_class_summary
textural_class_summary <- rbind(textural_class_summary, n=as.integer(sapply(unique(trafficability_0_10cm$textural_class_30cm), function(x) sum(!is.na(trafficability_0_10cm$FC_0.9[trafficability_0_10cm$textural_class_30cm==x])))))
textural_class_summary <- rbind(textural_class_summary, mean_clay=sapply(unique(trafficability_0_10cm$textural_class_30cm), function(x) mean(trafficability_0_10cm$clay_30cm[trafficability_0_10cm$textural_class_30cm==x], na.rm=TRUE)))
textural_class_summary <- rbind(textural_class_summary, FC_0.95=sapply(unique(trafficability_0_10cm$textural_class_30cm), function(x) mean(trafficability_0_10cm$FC_0.95[trafficability_0_10cm$textural_class_30cm==x], na.rm=TRUE)))
textural_class_summary <- rbind(textural_class_summary, FC_0.85=sapply(unique(trafficability_0_10cm$textural_class_30cm), function(x) mean(trafficability_0_10cm$FC_0.85[trafficability_0_10cm$textural_class_30cm==x], na.rm=TRUE)))
textural_class_summary <- rbind(textural_class_summary, FC_0.8=sapply(unique(trafficability_0_10cm$textural_class_30cm), function(x) mean(trafficability_0_10cm$FC_0.8[trafficability_0_10cm$textural_class_30cm==x], na.rm=TRUE)))
colnames(textural_class_summary)[order(textural_class_summary[8,])]
textural_class_summary <- textural_class_summary[,order(textural_class_summary[8,])]
row.names(textural_class_summary)
write.csv(textural_class_summary, file.path(tablesDir, 'summary_by_textural_class.csv'), row.names = TRUE)

plot(log(trafficability_0_10cm$ksat_30cm), trafficability_0_10cm$FC_0.9)
plot(trafficability_0_10cm$alpha_30cm, trafficability_0_10cm$FC_0.9)
plot(trafficability_0_10cm$clay_30cm, trafficability_0_10cm$FC_0.9)
plot(trafficability_0_10cm$sand_30cm, trafficability_0_10cm$FC_0.9)
plot(trafficability_0_10cm$silt_30cm, trafficability_0_10cm$FC_0.9)
plot(trafficability_0_10cm$theta_0.33b_30cm, trafficability_0_10cm$FC_0.9)
plot(trafficability_0_10cm$theta_15b_30cm, trafficability_0_10cm$FC_0.9)

summary(lm(FC_0.9 ~ clay_30cm, data = trafficability_0_10cm)) #Multiple R-squared:  0.4762
summary(lm(FC_0.9 ~ clay_30cm+sand_30cm, data = trafficability_0_10cm)) #Multiple R-squared:  0.4813
summary(lm(FC_0.9 ~ clay_30cm+sand_30cm+bd_30cm, data = trafficability_0_10cm)) #Multiple R-squared:  0.5059
summary(lm(FC_0.9 ~ clay_30cm+sand_30cm+bd_30cm+theta_0.33b_30cm, data = trafficability_0_10cm)) #Multiple R-squared:  0.5844
summary(lm(FC_0.9 ~ clay_30cm+sand_30cm+bd_30cm+theta_0.33b_30cm+theta_15b_30cm, data = trafficability_0_10cm)) #Multiple R-squared:  0.5988

write.csv(trafficability_0_10cm, file.path(tablesDir, 'trafficability_0_10cm_rosetta5_0_30cm_properties.csv'), row.names = TRUE)
