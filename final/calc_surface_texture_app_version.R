#moved code to "make_trafficability_maps.R" to add trafficability estimates to soil map unit shapfile  
library(raster)
library(rgeos)
library(aqp)
CIMISrawDir <- 'D:/Dissertation/Allowable_Depletion/SpatialCIMIS'
SAGBIdir <- 'C:/Users/smdevine/Desktop/SpatialData/SAGBI'
mainDir <- 'C:/Users/smdevine/Desktop/PostDoc'
trafficDir <- file.path(mainDir, 'trafficability')
SSURGOdir <- ssurgoDir <- file.path(trafficDir, 'soilweb app data')

#general functions
min_modified <- function(x) {
  if(all(is.na(x))) {
    return(NA)
  }
  else {min(x, na.rm = TRUE)}
}

#ssurgo functions
concat_names <- function(x, decat=FALSE) {
  if (decat) {x <- unlist(strsplit(x, split = '-'))}
  if(all(is.na(x))) {
    NA
  } else if (length(unique(x[!is.na(x)]))==1) {
    unique(x[!is.na(x)])
  } else {
    paste(unique(x[!is.na(x)])[order(unique(x[!is.na(x)]))], collapse = '-')
  }
}

#functions to work with ProfileApply
textural.class.calc <- function(sand, silt, clay, QC_param=1) {
  ifelse(is.na(sand) | is.na(silt) | is.na(clay), NA,
         ifelse(sand + silt + clay > (100 + QC_param) | sand + silt + clay < (100 - QC_param), paste('proportions do not sum to 100+-', QC_param), 
                ifelse(silt + 1.5 * clay < 15, 'sand',
                       ifelse(silt + 1.5 * clay >= 15 & silt + 2 * clay < 30, 'loamy sand',
                              ifelse((clay >= 7 & clay < 20 & sand > 52 & silt + 2 * clay >= 30) | (clay < 7 & silt < 50 & silt + 2 * clay >= 30), 'sandy loam',
                                     ifelse(clay >= 7 & clay < 27 & silt >=28 & silt < 50 & sand <= 52, 'loam',
                                            ifelse((silt >= 50 & clay >= 12 & clay < 27) | (silt >=50 & silt < 80 & clay < 12), 'silt loam',
                                                   ifelse(silt >= 80 & clay < 12, 'silt',
                                                          ifelse(clay >= 20 & clay < 35 & silt < 28 & sand > 45, 'sandy clay loam',
                                                                 ifelse(clay >= 27 & clay < 40 & sand > 20 & sand <= 45, 'clay loam',
                                                                        ifelse(clay >= 27 & clay < 40 & sand <= 20, 'silty clay loam',
                                                                               ifelse(clay >= 35 & sand > 45, 'sandy clay',
                                                                                      ifelse(clay >= 40 & silt >= 40, 'silty clay',
                                                                                             ifelse(clay >= 40 & sand <= 45 & silt < 40, 'clay','undefined textural class'))))))))))))))
}
wtd.mean <- function(x, y) {
  # use horizon thickness as a weight
  thick <- x$hzdepb_r - x$hzdept_r
  # compute the weighted mean, accounting for the possibility of missing data
  m <- weighted.mean(horizons(x)[[y]], w=thick, na.rm=TRUE)
  m
}
horizon_to_comp <- function(horizon_SPC, depth, comp_df, vars_of_interest = c('claytotal_r', 'silttotal_r', 'sandtotal_r'), varnames = c('clay', 'silt', 'sand'), SOC_content=FALSE, sum_AWC=FALSE) {
  columnames <- paste0(varnames, '_', depth, 'cm')
  print(cbind(vars_of_interest, columnames)) #show that it's all lined up
  assign("depth", depth, envir = .GlobalEnv) #this necessary because slice can't find the variable otherwise
  sliced_SPC <- slice(horizon_SPC, 0:(depth-1) ~ .) #depth was '0:depth' in previous version
  stopifnot(unique(sliced_SPC$pedon_key)==site(sliced_SPC)$pedon_key)
  for (i in seq_along(vars_of_interest)) {
    s <- site(sliced_SPC)
    s[[columnames[i]]] <- profileApply(sliced_SPC, FUN = wtd.mean, y=vars_of_interest[i])
    site(sliced_SPC) <- s
  }
  s <- site(sliced_SPC)
  if (SOC_content) {
    s[[paste0('kgOrg.m2_', depth, 'cm')]] <- profileApply(sliced_SPC, FUN = kgOrgC_sum)
    columnames <- c(columnames, paste0('kgOrg.m2_', depth, 'cm'))
  }
  
  if (sum_AWC) {
    s[[paste0('awc_', depth, 'cm')]] <- profileApply(sliced_SPC, FUN = awc_sum)
    columnames <- c(columnames, paste0('awc_', depth, 'cm'))
  }
  rm(depth, envir = .GlobalEnv) #because we had to put it there earlier
  s$compname <- comp_df$compname[match(s$cokey, comp_df$cokey)]
  s$mukey <- comp_df$mukey[match(s$cokey, comp_df$cokey)]
  s$comppct <- comp_df$comppct_r[match(s$cokey, comp_df$cokey)]
  s <- s[,c('mukey', 'cokey', 'compname', 'comppct', columnames)]
  s
}
MUaggregate <- function(df1, varname) {
  sapply(split(x=df1, f=df1$mukey), FUN=function(x) {if(sum(!is.na(x[[varname]]))==0) {NA} 
    else{sum(x$comppct[!is.na(x[[varname]])] * x[[varname]][!is.na(x[[varname]])] / sum(x$comppct[!is.na(x[[varname]])]))}
  })
}
MUAggregate_wrapper <- function(df1, varnames) {
  x <- sapply(varnames, FUN=MUaggregate, df1=df1)
  as.data.frame(cbind(mukey=as.integer(row.names(x)), x))
}

#read in soil data, provided by Mike W. on 4/19/22
list.files(ssurgoDir)
ca_mapunits <- shapefile(file.path(SSURGOdir, 'farmland_mapunits_2022.shp'))
crs(ca_mapunits) #+proj=aea +lat_0=0 +lon_0=-120 +lat_1=34 +lat_2=40.5 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=m +no_defs
ca_mapunits$area_ha <- area(ca_mapunits) / 10000
sum(ca_mapunits$area_ha) #7,334,505
sum(ca_mapunits$area_ha < 0.001) #57,426
names(ca_mapunits)
# ca_mapunits$mu_id <- 1:length(ca_mapunits)
# plot(ca_mapunits)

#read-in soil horizon texture data
farmland_texture <- read.csv(file.path(ssurgoDir, 'farmland_texture_data.csv'), stringsAsFactors = FALSE, na.strings = c('', ' '))

#add cimis cell numbers to map units
#export unique raster cell numbers of interest
#used spatialCIMIS_traffic.R to extract daily ET from these cells and export to consolidated csv
centroids_mu <- gCentroid(ca_mapunits, byid = TRUE)
centroids_mu_df <- SpatialPointsDataFrame(centroids_mu, data.frame(ca_mapunits))
sum(centroids_mu_df$area_ha)
length(unique(ca_mapunits$mukey)) #8847
spatialCIMIS <- raster(file.path(CIMISrawDir, 'ETo', '2009', 'ETo20090115.tif'))
as.character(crs(spatialCIMIS)) == as.character(crs(ca_mapunits))
CIMIScellnumbers <- as.integer(cellFromXY(object=spatialCIMIS, xy=centroids_mu_df))
length(CIMIScellnumbers)
length(unique(CIMIScellnumbers)) #25,478 CIMIS cells need sampling
CIMIScellunique_df <- data.frame(CIMIS_cells=unique(CIMIScellnumbers))
write.csv(CIMIScellunique_df, file.path(trafficDir, 'CIMIS', 'CIMIS_cells_unique_5.20.22.csv'), row.names = FALSE)

#add CIMIS_cells_unique
centroids_mu_df$CIMIScell <- cellFromXY(object=spatialCIMIS, xy=centroids_mu_df)
length(unique(centroids_mu_df$CIMIScell)) #25478
shapefile(centroids_mu_df, file.path(ssurgoDir, 'intermediate files', 'mu_centroids_cimis.shp'), overwrite=TRUE)
ca_mapunits$CIMIScell <- centroids_mu_df$CIMIScell[match(ca_mapunits$id, centroids_mu_df$id)]
# shapefile(ca_mapunits, file.path(ssurgoDir, 'intermediate files', 'mu_cimis_ref.shp')) #not necessary

#extract comp data from horizon data file provided by Mike W. on 4/19/22
colnames(farmland_texture)
comp_data <- data.frame(cokey=unique(farmland_texture$cokey))
comp_data$mukey <- farmland_texture$mukey[match(comp_data$cokey, farmland_texture$cokey)]
comp_data$compname <- farmland_texture$compname[match(comp_data$cokey, farmland_texture$cokey)]
comp_data$comppct_r <- farmland_texture$comppct_r[match(comp_data$cokey, farmland_texture$cokey)]
dim(comp_data)
head(comp_data)
length(unique(comp_data$mukey)) #8846
length(unique(ca_mapunits$mukey)) #8847
comp_data_AOI <- comp_data
# rm(comp_data)
# if(sum(is.na(comp_data_AOI$majcompflag)) > 0) {stop(print('there are NAs in majcomp column!'))}
if(sum(is.na(comp_data_AOI$comppct_r)) > 0) {stop(print('there are NAs in the comppct column!'))}
sum(is.na(comp_data_AOI$comppct_r)) #0
# comp_data_AOI <- comp_data_AOI[!is.na(comp_data_AOI$comppct_r),]
length(unique(comp_data_AOI$compname)) #2014 unique component names
# summary(comp_data_AOI$comppct_r[comp_data_AOI$majcompflag=='Yes'])
# sum(comp_data_AOI$comppct_r[comp_data_AOI$majcompflag=='Yes'] < 15) #18 instance of <15% comppct_r flagged as majcomps
# comp_data_AOI[comp_data_AOI$majcompflag=='Yes' & comp_data_AOI$comppct_r < 15,]
# sum(comp_data_AOI$comppct_r[comp_data_AOI$majcompflag=='No '] >= 15) #110 not flagged as majcomp in these valleys with > 15%
# comp_data_AOI[comp_data_AOI$majcompflag=='No ' & comp_data_AOI$comppct_r>=15, ]
# sum(comp_data_AOI$majcompflag=='No ' & comp_data_AOI$comppct_r>=15 & !is.na(comp_data_AOI$castorieindex)) #12 instances

#check length of each cokey
cokey_lengths <- unlist(lapply(strsplit(as.character(comp_data_AOI$cokey[1:2]), ""), length))
summary(cokey_lengths) #all 8 so slice won't create problems later
rm(cokey_lengths)

#define majcompflag
comp_data_AOI$majcompflag[comp_data_AOI$comppct_r>=15] <- 'Yes'
comp_data_AOI$majcompflag[comp_data_AOI$comppct_r < 15] <- 'No '

#aggregate horizon data to a 0-10 cm summary to define map unit and major component textural classes
colnames(farmland_texture)
horizon_data_AOI <- farmland_texture[,c(1:2,5:10)]

#horizon fixes identified as necessary below are applied now above
horizon_data_AOI$majcompflag <- comp_data_AOI$majcompflag[match(horizon_data_AOI$cokey, comp_data_AOI$cokey)]
table(horizon_data_AOI$majcompflag)
#convert to Soil Profile collection class with no minor comps per filtering below
horizon_AOI_majcomps <- horizon_data_AOI[horizon_data_AOI$majcompflag=='Yes', ]
# horizon_AOI_majcomps <- horizon_AOI_majcomps[!is.na(horizon_AOI_majcomps$hzname),] #removing these empty rows associated with non-soil major components eliminates error associated with "missing depths"; but these are kept to maintain full mukey representation in comp level table
depths(horizon_AOI_majcomps) <- cokey ~ hzdept_r + hzdepb_r
class(horizon_AOI_majcomps)
print(horizon_AOI_majcomps)
depth_ck <- checkHzDepthLogic(horizon_AOI_majcomps)
lapply(depth_ck, summary) #760 missing depth
lapply(depth_ck, table)
all(depth_ck$valid) #not all good
horizon_data_AOI[horizon_data_AOI$cokey==depth_ck$cokey[depth_ck$overlapOrGap==TRUE],] #this one error is below 10 cm, so can ignore
rm(depth_ck)

comp_AOI_10cm <- horizon_to_comp(horizon_SPC = horizon_AOI_majcomps, depth = 10, comp_df = comp_data_AOI)

comp_AOI_10cm$textural_class <- textural.class.calc(sand=comp_AOI_10cm$sand_10cm, silt = comp_AOI_10cm$silt_10cm, clay = comp_AOI_10cm$clay_10cm)
unique(comp_AOI_10cm$textural_class)
table(comp_AOI_10cm$textural_class)
sum(is.na(comp_AOI_10cm$textural_class)) #1344 are NA
unique(comp_AOI_10cm$compname[is.na(comp_AOI_10cm$textural_class)])
write.csv(comp_AOI_10cm, file.path(SSURGOdir, 'intermediate files', 'comp_AOI_10cm_5.20.22.csv'), row.names = FALSE) #look at removing NA to create table

colnames(comp_AOI_10cm)
length(unique(comp_AOI_10cm$mukey[!is.na(comp_AOI_10cm$clay_10cm)]))#8430 unique mukeys with components having at least texture data
length(unique(comp_AOI_10cm$cokey[!is.na(comp_AOI_10cm$clay_10cm)])) #10321 components
length(unique(comp_AOI_10cm$mukey[!is.na(comp_AOI_10cm$textural_class)])) #8406 mukeys with components having a textural class determination


#add metadata and restrictive layer info to shapefile
#clay
compnames_by_mukey <- data.frame(mukey=row.names(tapply(comp_data_AOI$compname, comp_data_AOI$mukey, function(x) unique(x))), compnames = as.character(tapply(comp_data_AOI$compname, comp_data_AOI$mukey, concat_names)), stringsAsFactors = FALSE)
length(unique(compnames_by_mukey$compnames)) #5284

majcompnames_by_mukey <- data.frame(mukey=row.names(tapply(comp_data_AOI$compname[comp_data_AOI$majcompflag=='Yes'], comp_data_AOI$mukey[comp_data_AOI$majcompflag=='Yes'], concat_names)), majcompnames=as.character(tapply(comp_data_AOI$compname[comp_data_AOI$majcompflag=='Yes'], comp_data_AOI$mukey[comp_data_AOI$majcompflag=='Yes'], concat_names)), stringsAsFactors = FALSE)
length(unique(majcompnames_by_mukey$majcompnames)) #2328 unique major component name combos

comppct_by_mukey_clay_data <- data.frame(mukey=row.names(tapply(comp_AOI_10cm$comppct[!is.na(comp_AOI_10cm$clay_10cm)], comp_AOI_10cm$mukey[!is.na(comp_AOI_10cm$clay_10cm)], sum)), comppct_tot=as.numeric(tapply(comp_AOI_10cm$comppct[!is.na(comp_AOI_10cm$clay_10cm)], comp_AOI_10cm$mukey[!is.na(comp_AOI_10cm$clay_10cm)], sum)), stringsAsFactors = FALSE)
summary(comppct_by_mukey_clay_data$comppct_tot)
sum(comppct_by_mukey_clay_data$comppct_tot < 70) #564 mukeys have less than 70% area with clay data

comppct_by_mukey_textural_data <- data.frame(mukey=row.names(tapply(comp_AOI_10cm$comppct[!is.na(comp_AOI_10cm$textural_class)], comp_AOI_10cm$mukey[!is.na(comp_AOI_10cm$textural_class)], sum)), comppct_tot=as.numeric(tapply(comp_AOI_10cm$comppct[!is.na(comp_AOI_10cm$textural_class)], comp_AOI_10cm$mukey[!is.na(comp_AOI_10cm$textural_class)], sum)), stringsAsFactors = FALSE)
summary(comppct_by_mukey_textural_data$comppct_tot)

domcomp_pct_by_mukey <- data.frame(mukey=row.names(tapply(comp_data_AOI$comppct_r, comp_data_AOI$mukey, function(x) max(x, na.rm=TRUE))), docomppct=as.numeric(tapply(comp_data_AOI$comppct_r, comp_data_AOI$mukey, function(x) max(x, na.rm=TRUE))), stringsAsFactors = FALSE)
summary(domcomp_pct_by_mukey$docomppct)

majcomp_pct_by_mukey <- data.frame(mukey=row.names(tapply(comp_data_AOI$comppct_r[comp_data_AOI$majcompflag=='Yes'], comp_data_AOI$mukey[comp_data_AOI$majcompflag=='Yes'], function(x) sum(x))), majcomppct=as.numeric(tapply(comp_data_AOI$comppct_r[comp_data_AOI$majcompflag=='Yes'], comp_data_AOI$mukey[comp_data_AOI$majcompflag=='Yes'], function(x) sum(x))), stringsAsFactors = FALSE)
summary(majcomp_pct_by_mukey$majcomppct)

majcomps_no_by_mukey <- data.frame(mukey=row.names(tapply(comp_data_AOI$majcompflag, comp_data_AOI$mukey, function(x) sum(x=='Yes'))), majcomp_no = as.numeric(tapply(comp_data_AOI$majcompflag, comp_data_AOI$mukey, function(x) sum(x=='Yes'))), stringsAsFactors = FALSE)
unique(majcomps_no_by_mukey$majcomp_no)
sum(majcomps_no_by_mukey$majcomp_no==0) #2
table(majcomps_no_by_mukey$majcomp_no)

#now add this summary info
ca_mapunits$mjcps_no <- majcomps_no_by_mukey$majcomp_no[match(ca_mapunits$mukey, majcomps_no_by_mukey$mukey)]
ca_mapunits$mjcmpnms <- majcompnames_by_mukey$majcompnames[match(ca_mapunits$mukey, majcompnames_by_mukey$mukey)]
ca_mapunits$dmcmp_pct <- domcomp_pct_by_mukey$docomppct[match(ca_mapunits$mukey, domcomp_pct_by_mukey$mukey)]
summary(ca_mapunits$dmcmp_pct)
ca_mapunits$mjcmp_pct <- majcomp_pct_by_mukey$majcomppct[match(ca_mapunits$mukey, majcomp_pct_by_mukey$mukey)]
summary(ca_mapunits$mjcmp_pct) #2 NAs
ca_mapunits$compct_clay <- comppct_by_mukey_clay_data$comppct_tot[match(ca_mapunits$mukey, comppct_by_mukey_clay_data$mukey)]
summary(ca_mapunits$compct_clay) #25696 are NA
ca_mapunits$compct_tex <- comppct_by_mukey_textural_data$comppct_tot[match(ca_mapunits$mukey, comppct_by_mukey_textural_data$mukey)]
summary(ca_mapunits$compct_tex) #26468

colnames(comp_AOI_10cm)
colnames(comp_AOI_10cm)[5:7]
AOI_10cm_muagg <- MUAggregate_wrapper(df1=comp_AOI_10cm, varnames = colnames(comp_AOI_10cm)[5:7])
AOI_10cm_muagg$textural_class <- textural.class.calc(sand = AOI_10cm_muagg$sand_10cm, silt = AOI_10cm_muagg$silt_10cm, clay = AOI_10cm_muagg$clay_10cm)
table(AOI_10cm_muagg$textural_class)
write.csv(AOI_10cm_muagg, file.path(SSURGOdir, 'intermediate files', 'AOI_10cm_muagg_5.20.22.csv'), row.names = FALSE)
colnames(AOI_10cm_muagg)
length(unique(AOI_10cm_muagg$mukey)) #8846 unique mukeys

ca_mapunits_10cm <- merge(ca_mapunits, AOI_10cm_muagg, by = 'mukey')
names(ca_mapunits_10cm)
sum(ca_mapunits_10cm$area_ha[!is.na(ca_mapunits_10cm$textural_class)]) #7116062
sum(ca_mapunits_10cm$area_ha) #7334505
7116062/7334505 #97% coverage

#acres by textural class
textural_class_area <- tapply(ca_mapunits_10cm$area_ha, ca_mapunits_10cm$textural_class, sum)
write.csv(textural_class_area, file.path(SSURGOdir, 'tables', 'texture_area.5.20.22.csv'), row.names = TRUE)

#export unique cimis_cell x textural class combos
ca_mapunits_10cm$climsoil_code <- paste0(ca_mapunits_10cm$CIMIScell, '_', ca_mapunits_10cm$textural_class)
shapefile(ca_mapunits_10cm, file.path(SSURGOdir, 'intermediate files', 'ca_mapunits_10cm.shp'), overwrite=TRUE)

climsoil_codes <- unique(ca_mapunits_10cm$climsoil_code)
length(climsoil_codes) #65630
climsoil_df <- data.frame(code=climsoil_codes, CIMIS=NA, textural_class=NA, stringsAsFactors = FALSE)
climsoil_df$CIMIS <- as.integer(sapply(climsoil_df$code, function(x) unlist(strsplit(x, '_'))[1]))
climsoil_df$textural_class <- sapply(climsoil_df$code, function(x) unlist(strsplit(x, '_'))[2])
unique(climsoil_df$textural_class)
table(climsoil_df$textural_class)
table(climsoil_df$textural_class[!climsoil_df$textural_class %in% c('proportions do not sum to 100+- 1', 'sandy clay', 'silt', 'NA')])
climsoil_df <- climsoil_df[!climsoil_df$textural_class %in% c('proportions do not sum to 100+- 1', 'sandy clay', 'silt', 'NA'), ]
dim(climsoil_df) #57804 to calc
write.csv(climsoil_df, file.path(SSURGOdir, 'tables', 'climsoil_unique_5.20.22.csv'), row.names = FALSE)
