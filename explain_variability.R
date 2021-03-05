dataDir <- 'C:/Users/smdevine/Desktop/PostDoc/trafficability'
soil_df <- read.csv(file.path(dataDir, 'soils_of_interest', 'soils_modeled_revised_QCpass_Oct2020.csv'), stringsAsFactors = FALSE)
dim(soil_df)
colnames(soil_df)
list.files(file.path(dataDir, 'climate_run_summaries', 'by_cell'))
cell_130212_results <- read.csv(file.path(dataDir, 'climate_run_summaries', 'by_cell', 'cell_130212_results.csv'), stringsAsFactors = FALSE)
dim(cell_130212_results)
cell_130212_results$textural_class <- NULL
soils_130212_results <- merge(soil_df, cell_130212_results, by='cokey')
tapply(soils_130212_results$fd_2009.01.15, soils_130212_results$textural_class, summary)
tapply(soils_130212_results$fd_2009.02.15, soils_130212_results$textural_class, summary)
tapply(soils_130212_results$fd_2009.03.15, soils_130212_results$textural_class, summary)
tapply(soils_130212_results$fd_2009.04.15, soils_130212_results$textural_class, summary)
