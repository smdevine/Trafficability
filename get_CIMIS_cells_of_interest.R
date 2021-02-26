climateDir <- 'D:/Allowable_Depletion/model_scaffold/run_model/Mar2018'
resultsDir <- 'D:/PostDoc/Trafficability/climate_runs/CIMIS_cell_selection'
CIMISrawDir <- 'D:/Dissertation/Allowable_Depletion/SpatialCIMIS'
library(raster)
ETo_raster <- raster(file.path(CIMISrawDir, 'ETo', '2009', 'ETo20090115.tif'))
ETo_cellnumbers <- ETo_raster
values(ETo_cellnumbers) <- 1:ncell(ETo_cellnumbers)
writeRaster(ETo_cellnumbers, file.path(resultsDir, 'CIMIS_cell_numbers.tif'), format='GTiff')
#create PET vector for testing general soil differences in trafficability
ETo <- read.csv(file.path(climateDir, 'SpatialCIMIS.ETo.QCpass.csv'), stringsAsFactors = FALSE)
# head(ETo$dates)
# tail(ETo$dates)
# unique(ETo$year)
ETo <- ETo[ETo$year %in% 2004:2017,]
# dim(ETo)
ETo_annual <- aggregate(ETo[,6:ncol(ETo)], by=list(year=ETo$year), FUN = sum)
# dim(ETo_annual)
# ETo_annual$cell_3091
ETo_annual_means <- apply(ETo_annual[,2:ncol(ETo_annual)], 1, mean)
names(ETo_annual_means) <- 2004:2017
ETo_annual_means[order(ETo_annual_means)]
median(ETo_annual_means) #either 2012 or 2016
mean(ETo_annual_means) #1393.3, closest to 2009
# conclusion: 2005=low ET, 2009=avg. ET, 2015=high ET
# 
quantiles_find <- quantile(ETo_annual[ETo_annual$year==2009,], probs = seq(from=5, to=95, by=5)/100)
ETo_annual_ranked <- as.data.frame(t(ETo_annual[,2:ncol(ETo_annual)]))
colnames(ETo_annual_ranked) <- paste0('yr_', 2004:2017)
ETo_annual_ranked <- ETo_annual_ranked[order(ETo_annual_ranked$yr_2009),]
# dim(ETo_annual_ranked)
rows_to_select <- round(c(0.05, seq(from=0.1, to=0.9, by=0.1), 0.95)*13060, digits = 0)
# rows_to_select
ETo_annual_ranked_selection <- ETo_annual_ranked[rows_to_select,]
CIMIS_cells <- rownames(ETo_annual_ranked_selection)
CIMIS_cells
write.csv(data.frame(CIMIS=CIMIS_cells, percentiles=c(0.05, seq(from=0.1, to=0.9, by=0.1), 0.95)), file.path(resultsDir, 'CIMIS_cells_percentiles.csv'), row.names=FALSE)
cellnumbers_of_interest <- as.integer(gsub('cell_', '', CIMIS_cells))
cellnumbers_of_interest_sp <- xyFromCell(ETo_raster, cell = cellnumbers_of_interest, spatial=TRUE)
cellnumbers_of_interest_sp$cellnumber <- cellnumbers_of_interest
plot(ETo_raster)
plot(cellnumbers_of_interest_sp, add=TRUE)
text(cellnumbers_of_interest_sp,  labels=cellnumbers_of_interest, cex=0.9, pos=1)
shapefile(cellnumbers_of_interest_sp, file.path(resultsDir, 'cellnumbers_of_interest.shp'), overwrite=TRUE)
#make summary by years (2005, 2009, and 2015) and months (Jan, Feb, Mar, Apr) used in modeling 
ETo_cells_of_interest <- read.csv(file.path(resultsDir, 'ETo_cells_of_interest.csv'), stringsAsFactors = FALSE)
colnames(ETo_cells_of_interest)
getETo_daily_mn <- function(df, cellname, start_date, days, days_to_flood=4) {
  flood_yr <- unlist(strsplit(start_date, '-'))[1]
  df <- df[, c('dates', cellname)]
  df <- df[which(df$dates==paste0('01_01_', flood_yr)):which(df$dates==paste0('12_31_', flood_yr)),2] / 10
  df[(as.integer(format(as.Date(start_date), format='%j'))+days_to_flood):(as.integer(format(as.Date(start_date), format='%j'))+days+days_to_flood)]
}
mean(getETo_daily_mn(ETo_cells_of_interest, 'cell_5248', '2005-01-15', 30))

#new cells of interest file
climateDir <- 'C:/Users/smdevine/Desktop/PostDoc/trafficability/CIMIS'
ETo <- read.csv(file.path(climateDir, 'SpatialCIMIS.ETo.QCpass.csv'), stringsAsFactors=FALSE)
ETo_annual <- aggregate(ETo[,6:ncol(ETo)], by=list(year=ETo$year), FUN = sum)
ETo_annual_ranked <- as.data.frame(t(ETo_annual[,2:ncol(ETo_annual)]))
colnames(ETo_annual_ranked) <- paste0('yr_', 2004:2018)
ETo_annual_ranked$meanETo <- apply(ETo_annual_ranked, 1, mean)
ETo_annual_ranked <- ETo_annual_ranked[order(ETo_annual_ranked$meanETo),]
write.csv(ETo_annual_ranked, file.path(climateDir, 'annual_ETo_ranked_AOI.csv'), row.names=TRUE)
colnames(ETo_cells_of_interest)
which(rownames(ETo_annual_ranked)=='cell_5248')/nrow(ETo_annual_ranked) #5248 was not in AOI
which(rownames(ETo_annual_ranked)=='cell_18033')/nrow(ETo_annual_ranked)
which(rownames(ETo_annual_ranked)=='cell_102682')/nrow(ETo_annual_ranked)
which(rownames(ETo_annual_ranked)=='cell_119481')/nrow(ETo_annual_ranked)
which(rownames(ETo_annual_ranked)=='cell_130212')/nrow(ETo_annual_ranked)
which(rownames(ETo_annual_ranked)=='cell_155771')/nrow(ETo_annual_ranked)
which(rownames(ETo_annual_ranked)=='cell_168006')/nrow(ETo_annual_ranked)
which(rownames(ETo_annual_ranked)=='cell_174576')/nrow(ETo_annual_ranked)
which(rownames(ETo_annual_ranked)=='cell_181807')/nrow(ETo_annual_ranked)
which(rownames(ETo_annual_ranked)=='cell_187903')/nrow(ETo_annual_ranked)
which(rownames(ETo_annual_ranked)=='cell_239134')/nrow(ETo_annual_ranked)

