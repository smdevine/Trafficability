#moved code to "make_trafficability_maps.R" to add trafficability estimates to soil map unit shapfile  
library(raster)
library(rgeos)
library(aqp)
CIMISrawDir <- 'D:/Dissertation/Allowable_Depletion/SpatialCIMIS'
SAGBIdir <- 'C:/Users/smdevine/Desktop/SpatialData/SAGBI'
mainDir <- 'C:/Users/smdevine/Desktop/PostDoc'
SSURGOdir <- ssurgoDir <- file.path(mainDir, 'soil health/ssurgo_data')
trafficDir <- file.path(mainDir, 'trafficability')
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
horizon_to_comp <- function(horizon_SPC, depth, comp_df, vars_of_interest = c('claytotal_r', 'silttotal_r', 'sandtotal_r', 'ksat_r'), varnames = c('clay', 'silt', 'sand', 'ksat'), SOC_content=FALSE, sum_AWC=FALSE) {
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

#read in soil data
ca_mapunits <- shapefile(file.path(ssurgoDir, 'ca_mapunits', 'ca_mapunits.shp'))
ca_mapunits$area_ha <- area(ca_mapunits) / 10000
sum(ca_mapunits$area_ha < 0.001) #only 23 less than 0.01 ha; 15 less than 0.001 ha
# crs(ca_mapunits)



#read in shapefile created in ArcGIS from sagbi x ca mapunit intersection
list.files(file.path(trafficDir, 'shapefiles'))
mu_sagbi <- shapefile(file.path(trafficDir, 'shapefiles', 'mu_sagbi.shp'))
2.47105 * sum(area(mu_sagbi) / 10000)
mu_sagbi #188319 features
mu_sagbi$area_ha <- area(mu_sagbi) / 10000
sum(mu_sagbi$area_ha) * 2.47105 #18101851 acres
sum(mu_sagbi$area_ha < 0.001) #4106 features < 0.01 ha; 1103 features < 0.01 ha

# shapefile(mu_sagbi, file.path(trafficDir, 'shapefiles', 'mu_sagbi.shp'), overwrite=TRUE)

#crop ca_mapunits with sagbiMod
# ca_mapunits_ca_ta_cropped <- crop(ca_mapunits_ca_ta, sagbiMod_dissolved) #had problems running
# ca_mapunits_ca_ta_SAGBI <- gIntersection(ca_mapunits_ca_ta, sagbiMod_dissolved, drop_lower_td=TRUE) #also failed

#read in mu data
ca_mu_data <- read.csv(file.path(ssurgoDir, 'ca_mapunit_data.csv'), stringsAsFactors = FALSE)
AOI_mu_data <- ca_mu_data[ca_mu_data$mukey %in% mu_sagbi$mukey,]

#read in comp and horizon data and limit to AOI
comp_data <- read.csv(file.path(ssurgoDir, 'component_data', 'ca_component_data.csv'), stringsAsFactors = FALSE, na.strings = c('', ' '))
dim(comp_data)
comp_data_AOI <- comp_data[comp_data$mukey %in% mu_sagbi$mukey,]
dim(comp_data_AOI)
length(unique(comp_data_AOI$mukey)) #8781 unique mukeys
if(sum(is.na(comp_data_AOI$majcompflag)) > 0) {stop(print('there are NAs in majcomp column!'))}
if(sum(is.na(comp_data_AOI$comppct_r)) > 0) {stop(print('there are NAs in the comppct column!'))}
sum(is.na(comp_data_AOI$comppct_r)) #3
comp_data_AOI <- comp_data_AOI[!is.na(comp_data_AOI$comppct_r),]
dim(comp_data_AOI) #41173 rows
length(unique(comp_data_AOI$mukey)) #8781 map units match above
length(unique(comp_data_AOI$cokey)) #41173 unique components
length(unique(comp_data_AOI$compname)) #2026 unique component names
summary(comp_data_AOI$comppct_r[comp_data_AOI$majcompflag=='Yes'])
sum(comp_data_AOI$comppct_r[comp_data_AOI$majcompflag=='Yes'] < 15) #18 instance of <15% comppct_r flagged as majcomps
comp_data_AOI[comp_data_AOI$majcompflag=='Yes' & comp_data_AOI$comppct_r < 15,]
sum(comp_data_AOI$comppct_r[comp_data_AOI$majcompflag=='No '] >= 15) #110 not flagged as majcomp in these valleys with > 15%
comp_data_AOI[comp_data_AOI$majcompflag=='No ' & comp_data_AOI$comppct_r>=15, ]
sum(comp_data_AOI$majcompflag=='No ' & comp_data_AOI$comppct_r>=15 & !is.na(comp_data_AOI$castorieindex)) #12 instances

#check length of each cokey
cokey_lengths <- unlist(lapply(strsplit(as.character(comp_data_AOI$cokey[1:2]), ""), length))
summary(cokey_lengths) #all 8 so slice won't create problems later
rm(cokey_lengths)

#fix a few majcompflag errors
comp_data_AOI$majcompflag[comp_data_AOI$majcompflag=='No ' & comp_data_AOI$comppct_r>=15] <- 'Yes'
comp_data_AOI$majcompflag[comp_data_AOI$majcompflag=='Yes' & comp_data_AOI$comppct_r < 15] <- 'No '



#aggregate horizon data to a 0-10 cm summary to map textural classes
horizon_data <- read.csv(file.path(ssurgoDir, 'ca_horizon_data.csv'), stringsAsFactors = FALSE)
horizon_data_AOI <- horizon_data[horizon_data$cokey %in% comp_data_AOI$cokey, ]
#horizon fixes identified as necessary
horizon_data_AOI[horizon_data_AOI$cokey==16387801,]
horizon_data_AOI[horizon_data_AOI$cokey==16387801 & horizon_data_AOI$hzname=='Ck','hzdept_r'] <- 114 #to fix horizonation to something acceptable
horizon_data_AOI[horizon_data_AOI$cokey==16597708, ] #same issue as above
horizon_data_AOI[horizon_data_AOI$cokey==16597708 & horizon_data_AOI$hzname=='Ck','hzdept_r'] <- 114
horizon_data_AOI$mukey <- comp_data_AOI$mukey[match(horizon_data_AOI$cokey, comp_data_AOI$cokey)]
horizon_data_AOI$majcompflag <- comp_data_AOI$majcompflag[match(horizon_data_AOI$cokey, comp_data_AOI$cokey)]
#convert to Soil Profile collection class with no minor comps per filtering above
horizon_AOI_majcomps <- horizon_data_AOI[horizon_data_AOI$majcompflag=='Yes', ]
depths(horizon_AOI_majcomps) <- cokey ~ hzdept_r + hzdepb_r
class(horizon_AOI_majcomps)
print(horizon_AOI_majcomps)
comp_AOI_10cm <- horizon_to_comp(horizon_SPC = horizon_AOI_majcomps, depth = 10, comp_df = comp_data_AOI)

comp_AOI_10cm$textural_class <- textural.class.calc(sand=comp_AOI_10cm$sand_10cm, silt = comp_AOI_10cm$silt_10cm, clay = comp_AOI_10cm$clay_10cm)
unique(comp_AOI_10cm$textural_class)
table(comp_AOI_10cm$textural_class)
write.csv(comp_AOI_10cm, file.path(trafficDir, 'ssurgo_intermediates', 'comp_AOI_10cm_2.9.21.csv'), row.names = FALSE)

comp_AOI_10cm <- read.csv(file.path(trafficDir, 'ssurgo_intermediates', 'comp_AOI_10cm_2.9.21.csv'))
colnames(comp_AOI_10cm)
length(unique(comp_AOI_10cm$mukey[!is.na(comp_AOI_10cm$clay_10cm)]))#8363 unique mukeys with components having at least texture data
length(unique(comp_AOI_10cm$cokey[!is.na(comp_AOI_10cm$clay_10cm)])) #10224


#add metadata and restrictive layer info to shapefile
#clay
compnames_by_mukey <- data.frame(mukey=row.names(tapply(comp_data_AOI$compname, comp_data_AOI$mukey, function(x) unique(x))), compnames = as.character(tapply(comp_data_AOI$compname, comp_data_AOI$mukey, concat_names)), stringsAsFactors = FALSE)
length(unique(compnames_by_mukey$compnames))

majcompnames_by_mukey <- data.frame(mukey=row.names(tapply(comp_data_AOI$compname[comp_data_AOI$majcompflag=='Yes'], comp_data_AOI$mukey[comp_data_AOI$majcompflag=='Yes'], concat_names)), majcompnames=as.character(tapply(comp_data_AOI$compname[comp_data_AOI$majcompflag=='Yes'], comp_data_AOI$mukey[comp_data_AOI$majcompflag=='Yes'], concat_names)), stringsAsFactors = FALSE)
length(unique(majcompnames_by_mukey$majcompnames)) #2330 unique major component name combos

comppct_by_mukey_clay_data <- data.frame(mukey=row.names(tapply(comp_AOI_10cm$comppct[!is.na(comp_AOI_10cm$clay_10cm)], comp_AOI_10cm$mukey[!is.na(comp_AOI_10cm$clay_10cm)], sum)), comppct_tot=as.numeric(tapply(comp_AOI_10cm$comppct[!is.na(comp_AOI_10cm$clay_10cm)], comp_AOI_10cm$mukey[!is.na(comp_AOI_10cm$clay_10cm)], sum)), stringsAsFactors = FALSE)
head(comppct_by_mukey_clay_data)
summary(comppct_by_mukey_clay_data$comppct_tot)
sum(comppct_by_mukey_clay_data$comppct_tot < 70) #574

mu_sagbi$muname <- AOI_mu_data$muname[match(mu_sagbi$mukey, AOI_mu_data$mukey)]
mu_sagbi$mjcmpnms <- majcompnames_by_mukey$majcompnames[match(mu_sagbi$mukey, majcompnames_by_mukey$mukey)]
mu_sagbi$complex <- ifelse(grepl('complex', mu_sagbi$muname), 'Yes', 'No')
table(mu_sagbi$complex)
mu_sagbi$associan <- ifelse(grepl('association', mu_sagbi$muname), 'Yes', 'No') #272 yes
table(mu_sagbi$associan)
mu_sagbi$dmcmp_pct <- domcomp_pct_by_mukey$docomppct[match(mu_sagbi$mukey, domcomp_pct_by_mukey$mukey)]
summary(mu_sagbi$dmcmp_pct)
mu_sagbi$mjcmp_pct <- majcomp_pct_by_mukey$majcomppct[match(mu_sagbi$mukey, majcomp_pct_by_mukey$mukey)]
summary(mu_sagbi$mjcmp_pct) #114 NAs
mu_sagbi$compct_clay <- comppct_by_mukey_clay_data$comppct_tot[match(mu_sagbi$mukey, comppct_by_mukey_clay_data$mukey)]

colnames(comp_AOI_10cm)
colnames(comp_AOI_10cm)[5:21]
AOI_10cm_muagg <- MUAggregate_wrapper(df1=comp_AOI_10cm, varnames = colnames(comp_AOI_10cm)[5:21])
AOI_10cm_muagg$textural_class <- textural.class.calc(sand = AOI_10cm_muagg$sand_10cm, silt = AOI_10cm_muagg$silt_10cm, clay = AOI_10cm_muagg$clay_10cm)
table(AOI_10cm_muagg$textural_class)
write.csv(AOI_10cm_muagg, file.path(trafficDir, 'ssurgo_intermediates', 'AOI_10cm_muagg.csv'), row.names = FALSE)
colnames(AOI_10cm_muagg)
AOI_10cm_muagg <- read.csv(file.path(trafficDir, 'ssurgo_intermediates', 'AOI_10cm_muagg.csv'))
length(unique(AOI_10cm_muagg$mukey)) #8560 unique mukeys

mu_sagbi_10cm <- merge(mu_sagbi, AOI_10cm_muagg, by = 'mukey')
names(mu_sagbi_10cm)
shapefile(mu_sagbi_10cm, file.path(trafficDir, 'shapefiles', 'mu_sagbi_10cm.shp'), overwrite=TRUE)

#acres by textural class
textural_class_area <- tapply(mu_sagbi_10cm$area_ha, mu_sagbi_10cm$textural_class, sum)
write.csv(textural_class_area, file.path(trafficDir, 'tables', 'texture_area_sagbi.csv'), row.names = TRUE)

#export unique cimis_cell x textural class combos
mu_sagbi_10cm$climsoil_code <- paste0(mu_sagbi_10cm$CIMIScell, '_', mu_sagbi_10cm$textural_class)
climsoil_codes <- unique(mu_sagbi_10cm$climsoil_code)
length(climsoil_codes)
climsoil_df <- data.frame(code=climsoil_codes, CIMIS=NA, textural_class=NA, stringsAsFactors = FALSE)
climsoil_df$CIMIS <- as.integer(sapply(climsoil_df$code, function(x) unlist(strsplit(x, '_'))[1]))
climsoil_df$textural_class <- sapply(climsoil_df$code, function(x) unlist(strsplit(x, '_'))[2])
unique(climsoil_df$textural_class)
table(climsoil_df$textural_class)
table(climsoil_df$textural_class[!climsoil_df$textural_class %in% c('proportions do not sum to 100+- 1', 'sandy clay', 'silt', 'NA')])
climsoil_df <- climsoil_df[!climsoil_df$textural_class %in% c('proportions do not sum to 100+- 1', 'sandy clay', 'silt', 'NA'), ]
dim(climsoil_df) #53709 to calc
write.csv(climsoil_df, file.path(trafficDir, 'climate_soil', 'climsoil_unique.csv'), row.names = FALSE)