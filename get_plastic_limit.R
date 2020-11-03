workDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/data from Stathis'
soilsDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/soils data'
tablesDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/Oct2020test/tables'
library(aqp)
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
mod_database <- read.csv(file.path(workDir, 'modelling_database.csv'), stringsAsFactors = FALSE)
mod_database <- mod_database[which(mod_database$hzn_top<200),] #delete horizons that start at or below 200 cm depth
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

final_cokeys <- as.integer(gsub('soil_', '', names(soil_data)))
SSURGO_final_cokeys <- final_cokeys[final_cokeys > 99999]
length(SSURGO_final_cokeys)
SSURGO_horizons_modeled <- mod_database[mod_database$cokey %in% SSURGO_final_cokeys, ]

soil_names <- unique(mod_database$taxonname[match(final_cokeys, mod_database$cokey)])
length(soil_names) #1322 unique names

#test
soilname_to_cokey <- function(x) {paste0("SELECT co.mukey, cokey, compname, comppct_r
FROM legend
INNER JOIN mapunit mu ON mu.lkey = legend.lkey
INNER JOIN component co ON mu.mukey = co.mukey 
WHERE
legend.areasymbol != 'US'
AND compname ='", x,  "';")}
soilname_queries <- lapply(soil_names, soilname_to_cokey)
comps_of_interest <- do.call(rbind, lapply(soilname_queries, SDA_query))
head(comps_of_interest)
comps_of_interest <- comps_of_interest[which(comps_of_interest$comppct_r >= 50), ]
length(comps_of_interest$cokey) #15308
write.csv(comps_of_interest, file.path(soilsDir, 'comps_of_interest.csv'), row.names = FALSE)

library(soilDB)
query_horizon <- function(x) { #this will not return cokey NAs
  print(x)
  SDA_query(paste0("SELECT comp.compname, comp.mukey, comp.cokey, ch.chkey, hzname, hzdept_r, hzdepb_r, awc_r, ec_r, claytotal_r, silttotal_r, sandtotal_r, dbthirdbar_r, wsatiated_r, wthirdbar_r, wfifteenbar_r, ksat_r, om_r, ll_r, pi_r
    FROM component comp
      LEFT OUTER JOIN chorizon ch on ch.cokey = comp.cokey
    WHERE comp.cokey = '", x, "'"))
}

SSURGO_horizons <- do.call(rbind, lapply(comps_of_interest$cokey, query_horizon))
write.csv(SSURGO_horizons, file.path(soilsDir, 'SSURGO_horizons_data.csv'), row.names=FALSE)
head(SSURGO_horizons)

#perform horizon-level with matching with VG parameter database
dim(SSURGO_horizons)
SSURGO_horizons_modeled[1,]
# i <- 1
# SSURGO_horizons[which(SSURGO_horizons$compname==SSURGO_horizons_modeled$taxonname[i]),]
# SSURGO_horizons[which(SSURGO_horizons$compname==SSURGO_horizons_modeled$taxonname[i] & SSURGO_horizons$hzname==SSURGO_horizons_modeled$hzn_desgn[i]), ]
SSURGO_horizons_modeled$ll_r <- NA
SSURGO_horizons_modeled$pi_r <- NA
for (i in seq_along(SSURGO_horizons_modeled$cokey)) {
  df <- SSURGO_horizons[which(SSURGO_horizons$compname==SSURGO_horizons_modeled$taxonname[i] & SSURGO_horizons$hzname==SSURGO_horizons_modeled$hzn_desgn[i] & SSURGO_horizons$hzdept_r==SSURGO_horizons_modeled$hzn_top[i] & SSURGO_horizons$hzdepb_r==SSURGO_horizons_modeled$hzn_bot[i] & SSURGO_horizons$claytotal_r==SSURGO_horizons_modeled$clay[i]), c('ll_r', 'pi_r')]
  SSURGO_horizons_modeled$ll_r[i] <- if (length(unique(df$ll_r))==1) {unique(df$ll_r)} else {NA}
  SSURGO_horizons_modeled$pi_r[i] <- if (length(unique(df$pi_r))==1) {unique(df$pi_r)} else {NA}
}
head(SSURGO_horizons_modeled, 10)
SSURGO_horizons_modeled$lpl_r <- SSURGO_horizons_modeled$ll_r - SSURGO_horizons_modeled$pi_r
summary(SSURGO_horizons_modeled$lpl_r) #5052 NAs

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

SSURGO_horizons_modeled$theta_fc <- theta_fc_est(SSURGO_horizons_modeled$teta_s, SSURGO_horizons_modeled$teta_r, SSURGO_horizons_modeled$Ks..cm.d., SSURGO_horizons_modeled$n)
summary(SSURGO_horizons_modeled$theta_fc)

SSURGO_horizons_modeled$h_fc <- -sapply(1:nrow(SSURGO_horizons_modeled), function(i) { optimize(find_h_at_theta_fc, interval = c(10,10000), theta_r = SSURGO_horizons_modeled$teta_r[i], theta_s = SSURGO_horizons_modeled$teta_s[i], alpha = SSURGO_horizons_modeled$alpha..1.cm.[i], n=SSURGO_horizons_modeled$n[i], theta_fc=SSURGO_horizons_modeled$theta_fc[i])$minimum})
summary(SSURGO_horizons_modeled$h_fc)

SSURGO_horizons_modeled$lpl_to_theta_fc <- (SSURGO_horizons_modeled$lpl_r/100) / SSURGO_horizons_modeled$theta_fc
summary(SSURGO_horizons_modeled$lpl_to_theta_fc)
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
SSURGO_horizons_modeled$textural_class <- textural.class.calc(sand = SSURGO_horizons_modeled$sand, silt = SSURGO_horizons_modeled$silt, clay = SSURGO_horizons_modeled$clay)
sum(is.na(SSURGO_horizons_modeled$textural_class))
sum(is.na(SSURGO_horizons_modeled$lpl_r)) #5052
write.csv(SSURGO_horizons_modeled, file.path(soilsDir, 'SSURGO_horizons_modeled_plasticity_data.csv'), row.names = FALSE)
SSURGO_horizons_modeled <- SSURGO_horizons_modeled[!is.na(SSURGO_horizons_modeled$lpl_to_theta_fc),]
dim(SSURGO_horizons_modeled) #11466 rows
sum(SSURGO_horizons_modeled$lpl_to_theta_fc <= 0) #786 less than 0
SSURGO_horizons_modeled <- SSURGO_horizons_modeled[!SSURGO_horizons_modeled$lpl_to_theta_fc <= 0,]
sum(SSURGO_horizons_modeled$lpl_to_theta_fc > 1) #2581 greater than 1
SSURGO_horizons_modeled <- SSURGO_horizons_modeled[!SSURGO_horizons_modeled$lpl_to_theta_fc > 1,]
sum(SSURGO_horizons_modeled$lpl_to_theta_fc <=0.5) #another 256 less than 0.5
SSURGO_horizons_modeled <- SSURGO_horizons_modeled[!SSURGO_horizons_modeled$lpl_to_theta_fc <= 0.5,]
tapply(SSURGO_horizons_modeled$lpl_to_theta_fc, SSURGO_horizons_modeled$textural_class, function(x) sum(!is.na(x)))
lpl_to_theata_fc_by_texture <- aggregate(SSURGO_horizons_modeled$lpl_to_theta_fc, list(SSURGO_horizons_modeled$textural_class), summary)
lpl_to_theata_fc_by_texture
lpl_to_theata_fc_by_texture$n <- tapply(SSURGO_horizons_modeled$lpl_to_theta_fc, SSURGO_horizons_modeled$textural_class, function(x) sum(!is.na(x)))
write.csv(lpl_to_theata_fc_by_texture, file.path(tablesDir, 'lpl_to_theta_fc_by_texture.csv'), row.names = TRUE)
