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
#this used to count up component percentages by reskind
reskind_comppct <- function(reskind, comp_df, reskind_df) {
  data.frame(mukey=row.names(tapply(comp_df$comppct_r[comp_df$cokey %in% reskind_df$cokey[grepl(reskind, reskind_df$reskinds)]], comp_df$mukey[comp_df$cokey %in% reskind_df$cokey[grepl(reskind, reskind_df$reskinds)]], sum)), compct_sum = as.numeric(tapply(comp_df$comppct_r[comp_df$cokey %in% reskind_df$cokey[grepl(reskind, reskind_df$reskinds)]], comp_df$mukey[comp_df$cokey %in% reskind_df$cokey[grepl(reskind, reskind_df$reskinds)]], sum)), stringsAsFactors = FALSE)
}

#functions to work with ProfileApply
kgOrgC_sum <- function(x, slice_it=FALSE, depth, rm.NAs=FALSE, om_to_c=1.72) {
  if (slice_it) {
    x <- horizons(x)[1:depth, ]
    depths(x) <- cokey ~ hzdept_r + hzdepb_r
  }
  thick <- x$hzdepb_r - x$hzdept_r
  sum((thick / 10) * (x$om_r / om_to_c) * x$dbthirdbar_r * (1 - x$fragvol_r_sum / 100), na.rm = rm.NAs)
}

awc_sum <- function(x, rm.NAs=TRUE) {
  thick <- x$hzdepb_r - x$hzdept_r
  sum(thick * x$awc_r, na.rm = rm.NAs)
}
wtd.mean <- function(x, y) {
  # use horizon thickness as a weight
  thick <- x$hzdepb_r - x$hzdept_r
  # compute the weighted mean, accounting for the possibility of missing data
  m <- weighted.mean(horizons(x)[[y]], w=thick, na.rm=TRUE)
  m
}
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
horizon_to_comp <- function(horizon_SPC, depth, comp_df, vars_of_interest = c('claytotal_r', 'silttotal_r', 'sandtotal_r', 'om_r', 'cec7_r', 'dbthirdbar_r', 'fragvol_r_sum', 'kwfact', 'ec_r', 'ph1to1h2o_r', 'sar_r', 'caco3_r', 'gypsum_r', 'lep_r', 'ksat_r'), varnames = c('clay', 'silt', 'sand', 'om', 'cec', 'bd', 'frags', 'kwf', 'ec', 'pH', 'sar', 'caco3', 'gyp', 'lep', 'ksat')) { #lep is linear extensibility
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
  s[[paste0('kgOrg.m2_', depth, 'cm')]] <- profileApply(sliced_SPC, FUN = kgOrgC_sum)
  s[[paste0('awc_', depth, 'cm')]] <- profileApply(sliced_SPC, FUN = awc_sum)
  columnames <- c(columnames, paste0('kgOrg.m2_', depth, 'cm'), paste0('awc_', depth, 'cm')) 
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

list.files(SAGBIdir)
sagbiMod <- shapefile(file.path(SAGBIdir, 'sagbi_mod.shp'))
names(sagbiMod)
sagbiUnMod <- shapefile(file.path(SAGBIdir, 'sagbi_unmod.shp'))
names(sagbiUnMod)
(sum(area(sagbiMod)) / 10000) * 2.47105 #18136950 acres
(sum(area(sagbiUnMod)) / 10000) * 2.47105 #same
summary(sagbiMod$sagbi)
crs(sagbiMod)

#dissolve sagbi mod into single polygons
sagbiMod_dissolved <- gUnaryUnion(sagbiMod)
sagbiMod_dissolved #1 feature
2.47105 * area(sagbiMod_dissolved) / 10000 #18123928 acres
shapefile(sagbiMod_dissolved, file.path(SAGBIdir, 'sagbi_mod_diss.shp'))
sagbiMod_dissolved <- shapefile(file.path(SAGBIdir, 'sagbi_mod_diss.shp'))
sagbiMod_dissolved <- shapefile(file.path(SAGBIdir, 'sagbi_mod_diss_v2.shp')) #arcgis produced same change in geometry


#read in soil data
ca_mapunits <- shapefile(file.path(ssurgoDir, 'ca_mapunits', 'ca_mapunits.shp'))
ca_mapunits$area_ha <- area(ca_mapunits) / 10000
sum(ca_mapunits$area_ha < 0.001) #only 23 less than 0.01 ha; 15 less than 0.001 ha
# crs(ca_mapunits)

# ca_mapunits_ca_ta <- spTransform(ca_mapunits, crs(sagbiMod))
# shapefile(ca_mapunits_ca_ta, file.path(ssurgoDir, 'ca_mapunits', 'ca_mapunits_CA_TA.shp'))
ca_mapunits_ca_ta <- shapefile(file.path(ssurgoDir, 'ca_mapunits', 'ca_mapunits_CA_TA.shp'))
ca_mapunits_ca_ta #464430 features

#read in shapefile created in ArcGIS from sagbi x ca mapunit intersection
list.files(file.path(trafficDir, 'shapefiles'))
mu_sagbi <- shapefile(file.path(trafficDir, 'shapefiles', 'mu_sagbi.shp'))
2.47105 * sum(area(mu_sagbi) / 10000)
mu_sagbi #188319 features
mu_sagbi$area_ha <- area(mu_sagbi) / 10000
sum(mu_sagbi$area_ha) * 2.47105 #18101851 acres
sum(mu_sagbi$area_ha < 0.001) #4106 features < 0.01 ha; 1103 features < 0.01 ha
centroids_mu <- gCentroid(mu_sagbi, byid = TRUE)
centroids_mu_df <- SpatialPointsDataFrame(centroids_mu, data.frame(mu_sagbi))
sum(centroids_mu_df$area_ha) * 2.47105
mukeys_unique <- unique(mu_sagbi$mukey)
length(mukeys_unique) #8782

#export unique raster cell numbers of interest
#used spatialCIMIS_traffic.R to extract daily ET from these cells and export to consolidated csv
spatialCIMIS <- raster(file.path(CIMISrawDir, 'ETo', '2009', 'ETo20090115.tif'))
CIMIScellnumbers <- as.integer(cellFromXY(object=spatialCIMIS, xy=centroids_mu_df))
length(CIMIScellnumbers)
length(unique(CIMIScellnumbers)) #24,561 CIMIS cells need sampling
CIMIScellunique_df <- data.frame(CIMIS_cells=unique(CIMIScellnumbers))
write.csv(CIMIScellunique_df, file.path(trafficDir, 'CIMIS', 'CIMIS_cells_unique.csv'), row.names = FALSE)

#add CIMIS_cells_unique
length(unique(mu_sagbi$FID_ca_map)) #188309
length(unique(centroids_mu_df$FID_ca_map)) #188309
centroids_mu_df$CIMIScell <- cellFromXY(object=spatialCIMIS, xy=centroids_mu_df)
length(unique(centroids_mu_df$CIMIScell)) #24561
shapefile(centroids_mu_df, file.path(trafficDir, 'shapefiles', 'mu_centroids_sagbi.shp'), overwrite=TRUE)
mu_sagbi$CIMIScell <- centroids_mu_df$CIMIScell[match(mu_sagbi$FID_ca_map, centroids_mu_df$FID_ca_map)]
shapefile(mu_sagbi, file.path(trafficDir, 'shapefiles', 'mu_sagbi.shp'), overwrite=TRUE)

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

restrictions <- read.csv(file.path(ssurgoDir, 'component_data', 'ca_restrictions.csv'), stringsAsFactors = FALSE, na.strings = c("", " "))

restrictions_AOI <- restrictions[restrictions$cokey %in% comp_data_AOI$cokey, ]
restrictions_AOI[is.na(restrictions_AOI$reskind),]
restrictions_AOI <- restrictions_AOI[!is.na(restrictions_AOI$reskind),]
restrictions_AOI$majcompflag <- comp_data_AOI$majcompflag[match(restrictions_AOI$cokey, comp_data_AOI$cokey)]
unique(restrictions_AOI$reskind)
restrictions_AOI$mukey <- comp_data_AOI$mukey[match(restrictions_AOI$cokey, comp_data_AOI$cokey)]
sum(duplicated(restrictions_AOI$cokey)) #551
restrictions_AOI$comppct <- comp_data_AOI$comppct_r[match(restrictions_AOI$cokey, comp_data_AOI$cokey)]

#depths to majcomp restrictions by reskind type
#which is necessary given that reskind has a NA
table(restrictions_AOI$reskind[restrictions_AOI$majcompflag=='Yes'])
restrictions_AOI[restrictions_AOI$reskind=='Densic bedrock',]
comp_data_AOI[comp_data_AOI$cokey %in% c(16560135, 16560147),]
AOI_mu_data[AOI_mu_data$mukey %in% c(1406599, 1406601),]

lithic_by_cokey <- restrictions_AOI[which(restrictions_AOI$reskind=='Lithic bedrock' & restrictions_AOI$majcompflag=='Yes'), ]
head(lithic_by_cokey)
# dim(lithic_by_cokey)
# sum(duplicated(lithic_by_cokey$cokey))
lithic_by_mukey <- MUAggregate_wrapper(df1 = lithic_by_cokey, varnames = c('resdept_r', 'resdepb_r'))
head(lithic_by_mukey)

paralithic_by_cokey <- restrictions_AOI[which(restrictions_AOI$reskind=='Paralithic bedrock' & restrictions_AOI$majcompflag=='Yes'), ]
# dim(paralithic_by_cokey)
# sum(duplicated(paralithic_by_cokey$cokey))
paralithic_by_mukey <- MUAggregate_wrapper(df1 = paralithic_by_cokey, varnames = c('resdept_r', 'resdepb_r'))
head(paralithic_by_mukey)

duripan_by_cokey <- restrictions_AOI[which(restrictions_AOI$reskind=='Duripan' & restrictions_AOI$majcompflag=='Yes'), ]
# dim(duripan_by_cokey)
# sum(duplicated(duripan_by_cokey$cokey))
duripan_by_mukey <- MUAggregate_wrapper(df1 = duripan_by_cokey, varnames = c('resdept_r', 'resdepb_r'))
head(duripan_by_mukey)

ATC_by_cokey <- restrictions_AOI[which(restrictions_AOI$reskind=='Abrupt textural change' & restrictions_AOI$majcompflag=='Yes'), ]
# dim(ATC_by_cokey)
# sum(duplicated(ATC_by_cokey$cokey))
ATC_by_mukey <- MUAggregate_wrapper(df1 = ATC_by_cokey, varnames = c('resdept_r', 'resdepb_r'))
head(ATC_by_mukey)

natric_by_cokey <- restrictions_AOI[which(restrictions_AOI$reskind=='Natric' & restrictions_AOI$majcompflag=='Yes'), ]
# dim(natric_by_cokey)
# sum(duplicated(natric_by_cokey$cokey))
natric_by_mukey <- MUAggregate_wrapper(df1 = natric_by_cokey, varnames = c('resdept_r', 'resdepb_r'))
natric_by_mukey

salic_by_cokey <- restrictions_AOI[which(restrictions_AOI$reskind=='Salic' & restrictions_AOI$majcompflag=='Yes'), ]
salic_by_mukey <- MUAggregate_wrapper(df1 = salic_by_cokey, varnames = c('resdept_r', 'resdepb_r'))
salic_by_mukey

SCTS_by_cokey <- restrictions_AOI[which(restrictions_AOI$reskind=='Strongly contrasting textural stratification' & restrictions_AOI$majcompflag=='Yes'), ]
SCTS_by_mukey <- MUAggregate_wrapper(df1 = SCTS_by_cokey, varnames = c('resdept_r', 'resdepb_r'))
SCTS_by_mukey

#now needs to be fixed so that if there's more than one restriction for a given cokey, the minimum is chosen before averaging with other cokeys
restrictions_AOI[restrictions_AOI$reskind=='Densic bedrock',]
misc_res_by_cokey <- restrictions_AOI[which(restrictions_AOI$reskind %in% c('Densic material', 'Cemented horizon', 'Petrocalcic', 'Densic bedrock') & restrictions_AOI$majcompflag=='Yes'), ]
# dim(misc_res_by_cokey) #
sum(duplicated(misc_res_by_cokey$cokey))
sum(duplicated(misc_res_by_cokey$mukey)) #5 will be averaged
misc_res_by_cokey[duplicated(misc_res_by_cokey$mukey),]
restrictions_AOI[restrictions_AOI$mukey %in% misc_res_by_cokey$mukey[duplicated(misc_res_by_cokey$mukey)],] #467123 has different depths for cokey 16620622
misc_res_by_mukey <- MUAggregate_wrapper(df1 = misc_res_by_cokey, varnames = c('resdept_r', 'resdepb_r'))
misc_res_by_mukey[misc_res_by_mukey$mukey==461133, ]

#consolidate reskinds by cokey and mukey
reskinds_by_cokey <- data.frame(cokey = row.names(tapply(restrictions_AOI$reskind[restrictions_AOI$majcompflag=='Yes'], restrictions_AOI$cokey[restrictions_AOI$majcompflag=='Yes'],  function(x) unique(x))), reskinds = as.character(tapply(restrictions_AOI$reskind[restrictions_AOI$majcompflag=='Yes'], restrictions_AOI$cokey[restrictions_AOI$majcompflag=='Yes'], concat_names)), stringsAsFactors = FALSE)

reskinds_by_cokey$mukey <- comp_data_AOI$mukey[match(reskinds_by_cokey$cokey, comp_data_AOI$cokey)]
length(unique(reskinds_by_cokey$mukey)) #4547 mukeys
reskinds_by_cokey$comp_pct <- comp_data_AOI$comppct_r[match(reskinds_by_cokey$cokey, comp_data_AOI$cokey)]
reskinds_by_cokey$majcompflag <- comp_data_AOI$majcompflag[match(reskinds_by_cokey$cokey, comp_data_AOI$cokey)]
#get rid of rock outcrop cokeys from this
reskinds_by_cokey$compname <- comp_data_AOI$compname[match(reskinds_by_cokey$cokey, comp_data_AOI$cokey)]
sum(reskinds_by_cokey$compname=='Rock outcrop') #hay 211
reskinds_by_cokey <- reskinds_by_cokey[reskinds_by_cokey$compname != 'Rock outcrop', ]
summary(as.factor(reskinds_by_cokey$majcompflag)) #minor components were left out in the original creation of the table above
summary(as.factor(tapply(reskinds_by_cokey$cokey, reskinds_by_cokey$mukey, function(x) length(x)))) #up to 4 per mukey, makes sense given that there are at most 3 major components in a map-unit for this AOI
#reskinds needs to be unconcatenated first before finding unique

reskinds_by_mukey <- data.frame(mukey = row.names(tapply(reskinds_by_cokey$cokey, reskinds_by_cokey$mukey, function(x) unique(x))), reskinds = as.character(tapply(reskinds_by_cokey$reskind, reskinds_by_cokey$mukey, concat_names, decat=TRUE)), stringsAsFactors = FALSE)
unique(reskinds_by_mukey$reskinds) #now this appears to be ok, using the decat=TRUE
table(reskinds_by_mukey$reskinds)

domcomp_pct_by_mukey <- data.frame(mukey=row.names(tapply(comp_data_AOI$comppct_r, comp_data_AOI$mukey, function(x) max(x, na.rm=TRUE))), docomppct=as.numeric(tapply(comp_data_AOI$comppct_r, comp_data_AOI$mukey, function(x) max(x, na.rm=TRUE))), stringsAsFactors = FALSE)
summary(domcomp_pct_by_mukey$docomppct)

majcomp_pct_by_mukey <- data.frame(mukey=row.names(tapply(comp_data_AOI$comppct_r[comp_data_AOI$majcompflag=='Yes'], comp_data_AOI$mukey[comp_data_AOI$majcompflag=='Yes'], function(x) sum(x))), majcomppct=as.numeric(tapply(comp_data_AOI$comppct_r[comp_data_AOI$majcompflag=='Yes'], comp_data_AOI$mukey[comp_data_AOI$majcompflag=='Yes'], function(x) sum(x))), stringsAsFactors = FALSE)
summary(majcomp_pct_by_mukey$majcomppct)

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
mu_sagbi$restrict <- reskinds_by_mukey$reskinds[match(mu_sagbi$mukey, reskinds_by_mukey$mukey)]
table(mu_sagbi$restrict)
sum(is.na(mu_sagbi$restrict))
mu_sagbi$restrict[is.na(mu_sagbi$restrict)] <- 'None'
mu_sagbi$Rock_OC <- ifelse(grepl('Rock outcrop', compnames_by_mukey$compnames[match(mu_sagbi$mukey, compnames_by_mukey$mukey)]), 'Yes', 'No')

mu_sagbi$Lithic <- ifelse(grepl('Lithic bedrock', reskinds_by_mukey$reskinds[match(mu_sagbi$mukey, reskinds_by_mukey$mukey)]), 'Yes', 'No') #632 were 'yes' after accounting for mukeys with more than one cokey with restrictions | mu_sagbi$Rock_OC=='Yes' add to conditional
table(mu_sagbi$Lithic)

mu_sagbi$Paralith <- ifelse(grepl('Paralithic bedrock', reskinds_by_mukey$reskinds[match(mu_sagbi$mukey, reskinds_by_mukey$mukey)]), 'Yes', 'No')
table(mu_sagbi$Paralith)

mu_sagbi$Duripan <- ifelse(grepl('Duripan', reskinds_by_mukey$reskinds[match(mu_sagbi$mukey, reskinds_by_mukey$mukey)]), 'Yes', 'No')
table(mu_sagbi$Duripan)

mu_sagbi$ATC <- ifelse(grepl('Abrupt textural change', reskinds_by_mukey$reskinds[match(mu_sagbi$mukey, reskinds_by_mukey$mukey)]), 'Yes', 'No') #ATC=abrupt textural change
table(mu_sagbi$ATC)

mu_sagbi$Natric <- ifelse(grepl('Natric', reskinds_by_mukey$reskinds[match(mu_sagbi$mukey, reskinds_by_mukey$mukey)]), 'Yes', 'No')
table(mu_sagbi$Natric)

mu_sagbi$Salic <- ifelse(grepl('Salic', reskinds_by_mukey$reskinds[match(mu_sagbi$mukey, reskinds_by_mukey$mukey)]), 'Yes', 'No')
table(mu_sagbi$Salic)

mu_sagbi$SCTS <- ifelse(grepl('Strongly contrasting textural stratification', reskinds_by_mukey$reskinds[match(mu_sagbi$mukey, reskinds_by_mukey$mukey)]), 'Yes', 'No')
table(mu_sagbi$SCTS)

mu_sagbi$Misc_Res <- ifelse(grepl('Densic material|Cemented horizon|Petrocalcic', reskinds_by_mukey$reskinds[match(mu_sagbi$mukey, reskinds_by_mukey$mukey)]), 'Yes', 'No')
table(mu_sagbi$Misc_Res)

AOI_lithic_comppct <- data.frame(mukey=row.names(tapply(comp_data_AOI$comppct_r[comp_data_AOI$cokey %in% reskinds_by_cokey$cokey[grepl('Lithic', reskinds_by_cokey$reskinds)] & comp_data_AOI$compname != 'Rock outcrop'], comp_data_AOI$mukey[comp_data_AOI$cokey %in% reskinds_by_cokey$cokey[grepl('Lithic', reskinds_by_cokey$reskinds)] & comp_data_AOI$compname != 'Rock outcrop'], sum)), compct_sum = as.numeric(tapply(comp_data_AOI$comppct_r[comp_data_AOI$cokey %in% reskinds_by_cokey$cokey[grepl('Lithic', reskinds_by_cokey$reskinds)] & comp_data_AOI$compname != 'Rock outcrop'], comp_data_AOI$mukey[comp_data_AOI$cokey %in% reskinds_by_cokey$cokey[grepl('Lithic', reskinds_by_cokey$reskinds)] & comp_data_AOI$compname != 'Rock outcrop'], sum)), stringsAsFactors = FALSE)
head(AOI_lithic_comppct)

AOI_rockOC_comppct <- data.frame(mukey=row.names(tapply(comp_data_AOI$comppct_r[comp_data_AOI$compname=='Rock outcrop'], comp_data_AOI$mukey[comp_data_AOI$compname=='Rock outcrop'], function(x) sum(x, na.rm = TRUE))), compct_sum = as.numeric(tapply(comp_data_AOI$comppct_r[comp_data_AOI$compname=='Rock outcrop'], comp_data_AOI$mukey[comp_data_AOI$compname=='Rock outcrop'], function(x) sum(x, na.rm = TRUE))), stringsAsFactors = FALSE)
head(AOI_rockOC_comppct)

#paralithic component pct
AOI_paralithic_comppct <- reskind_comppct(reskind = 'Paralithic', comp_df = comp_data_AOI, reskind_df = reskinds_by_cokey)

#duripan component pct
AOI_duripan_comppct <- reskind_comppct(reskind = 'Duripan', comp_df = comp_data_AOI, reskind_df = reskinds_by_cokey)

#abrupt textural contrast (ATC) comppct
AOI_ATC_comppct <- reskind_comppct(reskind = 'Abrupt textural change', comp_df = comp_data_AOI, reskind_df = reskinds_by_cokey)

#natric comppct
AOI_Natric_comppct <- reskind_comppct(reskind = 'Natric', comp_df = comp_data_AOI, reskind_df = reskinds_by_cokey)
dim(AOI_Natric_comppct)
summary(AOI_Natric_comppct$compct_sum)

##Salic comppct
AOI_Salic_comppct <- reskind_comppct(reskind = 'Salic', comp_df = comp_data_AOI, reskind_df = reskinds_by_cokey)
dim(AOI_Salic_comppct)
AOI_Salic_comppct

#SCTS comppct
AOI_SCTS_comppct <- reskind_comppct(reskind = 'Strongly contrasting textural stratification', comp_df = comp_data_AOI, reskind_df = reskinds_by_cokey)
AOI_SCTS_comppct

AOI_Misc_comppct <- reskind_comppct(reskind = 'Densic material|Cemented horizon|Petrocalcic', comp_df = comp_data_AOI, reskind_df = reskinds_by_cokey)
AOI_Misc_comppct

#add reskind comppct to mapunit
mu_sagbi$Lthc_pct <- 0
mu_sagbi$Lthc_pct[mu_sagbi$Lithic=='Yes'] <- AOI_lithic_comppct$compct_sum[match(mu_sagbi$mukey[mu_sagbi$Lithic=='Yes'], AOI_lithic_comppct$mukey)]
# summary(mu_sagbi$Lthc_pct)
mu_sagbi$RckOC_pct <- 0
mu_sagbi$RckOC_pct[mu_sagbi$Rock_OC=='Yes'] <- AOI_rockOC_comppct$compct_sum[match(mu_sagbi$mukey[mu_sagbi$Rock_OC=='Yes'], AOI_rockOC_comppct$mukey)]
# summary(mu_sagbi$RckOC_pct)
# summary(rowSums(as.data.frame(mu_sagbi[c('Lthc_pct', 'RckOC_pct')])))
mu_sagbi$Plth_pct <- 0
mu_sagbi$Plth_pct[mu_sagbi$Paralith=='Yes'] <- AOI_paralithic_comppct$compct_sum[match(mu_sagbi$mukey[mu_sagbi$Paralith=='Yes'], AOI_paralithic_comppct$mukey)]
# summary(mu_sagbi$Plth_pct)
mu_sagbi$Drpn_pct <- 0
mu_sagbi$Drpn_pct[mu_sagbi$Duripan=='Yes'] <- AOI_duripan_comppct$compct_sum[match(mu_sagbi$mukey[mu_sagbi$Duripan=='Yes'], AOI_duripan_comppct$mukey)]
# summary(mu_sagbi$Drpn_pct)
mu_sagbi$ATC_pct <- 0
mu_sagbi$ATC_pct[mu_sagbi$ATC=='Yes'] <- AOI_ATC_comppct$compct_sum[match(mu_sagbi$mukey[mu_sagbi$ATC=='Yes'], AOI_ATC_comppct$mukey)]
# summary(mu_sagbi$ATC_pct)
mu_sagbi$Natr_pct <- 0
mu_sagbi$Natr_pct[mu_sagbi$Natric=='Yes'] <- AOI_Natric_comppct$compct_sum[match(mu_sagbi$mukey[mu_sagbi$Natric=='Yes'], AOI_Natric_comppct$mukey)]
# summary(mu_sagbi$Natr_pct)
mu_sagbi$Salc_pct <- 0
mu_sagbi$Salc_pct[mu_sagbi$Salic=='Yes'] <- AOI_Salic_comppct$compct_sum[match(mu_sagbi$mukey[mu_sagbi$Salic=='Yes'], AOI_Salic_comppct$mukey)]
mu_sagbi$SCTS_pct <- 0
mu_sagbi$SCTS_pct[mu_sagbi$SCTS=='Yes'] <- AOI_SCTS_comppct$compct_sum[match(mu_sagbi$mukey[mu_sagbi$SCTS=='Yes'], AOI_SCTS_comppct$mukey)]
mu_sagbi$MRes_pct <- 0
mu_sagbi$MRes_pct[mu_sagbi$Misc_Res=='Yes'] <- AOI_Misc_comppct$compct_sum[match(mu_sagbi$mukey[mu_sagbi$Misc_Res=='Yes'], AOI_Misc_comppct$mukey)]

#add restriction depth info
mu_sagbi$Lthc_dep <- NA
mu_sagbi$Lthc_dep[mu_sagbi$Lithic=='Yes'] <- lithic_by_mukey$resdept_r[match(mu_sagbi$mukey[mu_sagbi$Lithic=='Yes'], lithic_by_mukey$mukey)]
summary(mu_sagbi$Lthc_dep)
mu_sagbi$Plth_dep <- NA
mu_sagbi$Plth_dep[mu_sagbi$Paralith=='Yes'] <- paralithic_by_mukey$resdept_r[match(mu_sagbi$mukey[mu_sagbi$Paralith=='Yes'], paralithic_by_mukey$mukey)]
summary(mu_sagbi$Plth_dep)
mu_sagbi$Drpn_dep <- NA
mu_sagbi$Drpn_dep[mu_sagbi$Duripan=='Yes'] <- duripan_by_mukey$resdept_r[match(mu_sagbi$mukey[mu_sagbi$Duripan=='Yes'], duripan_by_mukey$mukey)]
summary(mu_sagbi$Drpn_dep)
mu_sagbi$ATC_dep <- NA
mu_sagbi$ATC_dep[mu_sagbi$ATC=='Yes'] <- ATC_by_mukey$resdept_r[match(mu_sagbi$mukey[mu_sagbi$ATC=='Yes'], ATC_by_mukey$mukey)]
summary(mu_sagbi$ATC_dep)
mu_sagbi$Natr_dep <- NA
mu_sagbi$Natr_dep[mu_sagbi$Natric=='Yes'] <- natric_by_mukey$resdept_r[match(mu_sagbi$mukey[mu_sagbi$Natric=='Yes'], natric_by_mukey$mukey)]
summary(mu_sagbi$Natr_dep)
mu_sagbi$Salc_dep <- NA
mu_sagbi$Salc_dep[mu_sagbi$Salic=='Yes'] <- salic_by_mukey$resdept_r[match(mu_sagbi$mukey[mu_sagbi$Salic=='Yes'], salic_by_mukey$mukey)]
summary(mu_sagbi$Salc_dep)
mu_sagbi$SCTS_dep <- NA
mu_sagbi$SCTS_dep[mu_sagbi$SCTS=='Yes'] <- SCTS_by_mukey$resdept_r[match(mu_sagbi$mukey[mu_sagbi$SCTS=='Yes'], SCTS_by_mukey$mukey)]
summary(mu_sagbi$SCTS_dep)
mu_sagbi$MRes_dep <- NA
mu_sagbi$MRes_dep[mu_sagbi$Misc_Res=='Yes'] <- misc_res_by_mukey$resdept_r[match(mu_sagbi$mukey[mu_sagbi$Misc_Res=='Yes'], misc_res_by_mukey$mukey)]
summary(mu_sagbi$MRes_dep)

#minimum of all reskind depths
mu_sagbi$MnRs_dep <- apply(as.data.frame(mu_sagbi)[ ,c('Lthc_dep', 'Plth_dep', 'Drpn_dep', 'ATC_dep', 'Natr_dep', 'Salc_dep', 'SCTS_dep', 'MRes_dep')], 1, min_modified)
summary(mu_sagbi$MnRs_dep)

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