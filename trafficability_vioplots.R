mainDir <- 'C:/Users/smdevine/Desktop/PostDoc/trafficability'
FiguresDir <- file.path(mainDir, 'Figures')
resultsDir <- 'D:/PostDoc/Trafficability/climate_runs'
MonthlyMeansDir <- file.path(resultsDir, 'overall_summaries/by_cell/monthly_means')
library(vioplot)
library(extrafont)
library(extrafontdb)
#font_import() #only needs to be done one time after updating and re-installing R and moving and updating packages
loadfonts(device = 'win')
texture_colors <- data.frame(textures=c('clay', 'silty clay' , 'silty clay loam', 'clay loam', 'silt loam', 'sandy clay loam', 'loam', 'sandy loam', 'loamy sand', 'sand'), texture_labs=c('clay', 'silty\nclay' , 'silty\nclay\nloam', 'clay\nloam', 'silt\nloam', 'sandy\nclay\nloam', 'loam', 'sandy\nloam', 'loamy\nsand', 'sand'), red=c(169, 0, 0, 223, 0, 170, 255, 230, 115, 255), green=c(0, 112, 197, 115, 168, 255, 0, 152, 76, 255), blue=c(230, 255, 255, 255, 132, 0, 0, 0, 0, 0), stringsAsFactors = FALSE)


# fname <- 'cell_181807_Jan.csv'
# ylim_vioplot <- c(0,75)
# month <- 'January'
# col_index <- c(6,8:10)

vioplot_traffic_all <- function(fname,  ylim_vioplot, mar, fig_label, fig_height, fig_width, month, savename, axis_line) {
  df <- read.csv(file.path(MonthlyMeansDir, fname), stringsAsFactors = FALSE)
  df <- df[df$textural_class %in% texture_colors$textures,]
  tiff(file = file.path(FiguresDir, 'violin plots', savename), family = 'Times New Roman', width = fig_width, height = fig_height, pointsize = 12, units = 'in', res=800, compression='lzw')
  par(mar=mar)
  vioplot(df[,4][df$textural_class==texture_colors$textures[1]], df[,4][df$textural_class==texture_colors$textures[2]], df[,4][df$textural_class==texture_colors$textures[3]], df[,4][df$textural_class==texture_colors$textures[4]], df[,4][df$textural_class==texture_colors$textures[5]], df[,4][df$textural_class==texture_colors$textures[6]], df[,4][df$textural_class==texture_colors$textures[7]], df[,4][df$textural_class==texture_colors$textures[8]], df[,4][df$textural_class==texture_colors$textures[9]], df[,4][df$textural_class==texture_colors$textures[10]], col=rgb(texture_colors$red/255, texture_colors$green/255, texture_colors$blue/255), wex=1.2, cex=0.8,  rectCol = 'gray', ylim = ylim_vioplot, ylab = NULL, xaxt = 'n')
  axis(side = 1, at=1:10, texture_colors$texture_labs, line=axis_line, tick = FALSE, cex.axis=0.8)
  mtext('Textural classes (0-10 cm)', side = 1, line = 1)
  mtext(paste0('Days to trafficability (', month, ')'), side=2, line=2.5)
  # text(x=7.5, y=180, labels=fig_label)
  dev.off()
}

vioplot_traffic_all(fname = 'cell_181807_Jan.csv', ylim_vioplot = c(0,100), mar = c(2, 4, 0.5, 0.5), fig_height = 4, fig_width=6.5, month='January', savename = 'cell_181807_Jan.tif', axis_line = -2) #fname,  ylim_vioplot, mar, fig_label, fig_height, month, savename

fname <- 'cell_181807_Feb.csv'
ylim_vioplot <- c(0,50)
month <- 'February'
col_index <- 3:10

vioplot_traffic_loamy_coarse <- function(fname, to_plot, ylim_vioplot, mar, fig_label, fig_height, fig_width, month, savename, axis_line, col_index, lab_ht, ylab) {
  df <- read.csv(file.path(MonthlyMeansDir, fname), stringsAsFactors = FALSE)
  df <- df[df$textural_class %in% texture_colors$textures,]
  tiff(file = file.path(FiguresDir, 'violin plots', savename), family = 'Times New Roman', width = fig_width, height = fig_height, pointsize = 12, units = 'in', res=800, compression='lzw')
  par(mar=mar)
  vioplot(df[,4][df$textural_class==texture_colors$textures[3]], df[,4][df$textural_class==texture_colors$textures[4]], df[,4][df$textural_class==texture_colors$textures[5]], df[,4][df$textural_class==texture_colors$textures[6]], df[,4][df$textural_class==texture_colors$textures[7]], df[,4][df$textural_class==texture_colors$textures[8]], df[,4][df$textural_class==texture_colors$textures[9]], df[,4][df$textural_class==texture_colors$textures[10]], col=rgb(texture_colors$red/255, texture_colors$green/255, texture_colors$blue/255)[col_index], wex=1, cex=0.3,  rectCol = 'gray', pchMed = 1, colMed = 'black', ylim = ylim_vioplot, ylab = NULL, xaxt = 'n')
  axis(side = 1, at=1:length(col_index), texture_colors$texture_labs[col_index], line=axis_line, tick = FALSE, cex.axis=0.7)
  # mtext('Textural classes, <40% clay (0-10 cm)', side = 1, line = 0.5)
  if(ylab) {mtext(paste0('Days to trafficability after flooding'), side=2, line=2.25)}
  text(x=8, y=lab_ht, labels=month)
  dev.off()
}
#Tulare region in Kern County
vioplot_traffic_loamy_coarse(fname = 'cell_181807_Jan.csv', ylim_vioplot = c(0,65), mar = c(0.1, 3.5, 0.1, 0.1), fig_height = 3.25, fig_width=4.75, month='January', savename = 'cell_181807_Jan_loamy_coarse.tif', axis_line = -2, col_index=3:10, lab_ht = 65, ylab=TRUE)
vioplot_traffic_loamy_coarse(fname = 'cell_181807_Feb.csv', ylim_vioplot = c(0,65), mar = c(0.1, 3.5, 0.1, 0.1), fig_height = 3.25, fig_width=4.75, month='February', savename = 'cell_181807_Feb_loamy_coarse.tif', axis_line = -2, col_index=3:10, lab_ht = 65, ylab=FALSE)
vioplot_traffic_loamy_coarse(fname = 'cell_181807_Mar.csv', ylim_vioplot = c(0,45), mar = c(0.1, 3.5, 0.1, 0.1), fig_height = 3.25, fig_width=4.75, month='March', savename = 'cell_181807_Mar_loamy_coarse.tif', axis_line = -2, col_index=3:10, lab_ht = 45, ylab = TRUE)
vioplot_traffic_loamy_coarse(fname = 'cell_181807_Apr.csv', ylim_vioplot = c(0,45), mar = c(0.1, 3.5, 0.1, 0.1), fig_height = 3.25, fig_width=4.75, month='April', savename = 'cell_181807_Apr_loamy_coarse.tif', axis_line = -2, col_index=3:10, lab_ht = 45, ylab = FALSE)

#near Salinas
vioplot_traffic_loamy_coarse(fname = 'cell_174576_Jan.csv', ylim_vioplot = c(0,65), mar = c(0.1, 3.5, 0.1, 0.1), fig_height = 3.25, fig_width=4.75, month='January', savename = 'cell_174576_Jan_loamy_coarse.tif', axis_line = -2, col_index=3:10, lab_ht = 65, ylab=TRUE)
vioplot_traffic_loamy_coarse(fname = 'cell_174576_Feb.csv', ylim_vioplot = c(0,65), mar = c(0.1, 3.5, 0.1, 0.1), fig_height = 3.25, fig_width=4.75, month='February', savename = 'cell_174576_Feb_loamy_coarse.tif', axis_line = -2, col_index=3:10, lab_ht = 65, ylab=FALSE)
vioplot_traffic_loamy_coarse(fname = 'cell_174576_Mar.csv', ylim_vioplot = c(0,45), mar = c(0.1, 3.5, 0.1, 0.1), fig_height = 3.25, fig_width=4.75, month='March', savename = 'cell_174576_Mar_loamy_coarse.tif', axis_line = -2, col_index=3:10, lab_ht = 45, ylab = TRUE)
vioplot_traffic_loamy_coarse(fname = 'cell_174576_Apr.csv', ylim_vioplot = c(0,45), mar = c(0.1, 3.5, 0.1, 0.1), fig_height = 3.25, fig_width=4.75, month='April', savename = 'cell_174576_Apr_loamy_coarse.tif', axis_line = -2, col_index=3:10, lab_ht = 45, ylab = FALSE)

#near Stockton
vioplot_traffic_loamy_coarse(fname = 'cell_119481_Jan.csv', ylim_vioplot = c(0,65), mar = c(0.1, 3.5, 0.1, 0.1), fig_height = 3.25, fig_width=4.75, month='January', savename = 'cell_119481_Jan_loamy_coarse.tif', axis_line = -2, col_index=3:10, lab_ht = 65, ylab=TRUE)
vioplot_traffic_loamy_coarse(fname = 'cell_119481_Feb.csv', ylim_vioplot = c(0,65), mar = c(0.1, 3.5, 0.1, 0.1), fig_height = 3.25, fig_width=4.75, month='February', savename = 'cell_119481_Feb_loamy_coarse.tif', axis_line = -2, col_index=3:10, lab_ht = 65, ylab=FALSE)
vioplot_traffic_loamy_coarse(fname = 'cell_119481_Mar.csv', ylim_vioplot = c(0,45), mar = c(0.1, 3.5, 0.1, 0.1), fig_height = 3.25, fig_width=4.75, month='March', savename = 'cell_119481_Mar_loamy_coarse.tif', axis_line = -2, col_index=3:10, lab_ht = 45, ylab = TRUE)
vioplot_traffic_loamy_coarse(fname = 'cell_119481_Apr.csv', ylim_vioplot = c(0,45), mar = c(0.1, 3.5, 0.1, 0.1), fig_height = 3.25, fig_width=4.75, month='April', savename = 'cell_119481_Apr_loamy_coarse.tif', axis_line = -2, col_index=3:10, lab_ht = 45, ylab = FALSE)

#fine textures
vioplot_traffic_fine <- function(fname, to_plot, ylim_vioplot, mar, fig_label, fig_height, fig_width, month, savename, axis_line,  lab_ht, ylab) {
  df <- read.csv(file.path(MonthlyMeansDir, fname), stringsAsFactors = FALSE)
  df <- df[df$textural_class %in% texture_colors$textures,]
  tiff(file = file.path(FiguresDir, 'violin plots', savename), family = 'Times New Roman', width = fig_width, height = fig_height, pointsize = 12, units = 'in', res=800, compression='lzw')
  par(mar=mar)
  vioplot(df[,4][df$textural_class==texture_colors$textures[1]], df[,4][df$textural_class==texture_colors$textures[2]], col=rgb(texture_colors$red/255, texture_colors$green/255, texture_colors$blue/255)[1:2], wex=1.2, cex=0.8,  rectCol = 'gray', ylim = ylim_vioplot, ylab = NULL, xaxt = 'n')
  axis(side = 1, at=1:2, texture_colors$texture_labs[1:2], line=axis_line, tick = FALSE, cex.axis=0.7)
  # mtext('Textural classes, <40% clay (0-10 cm)', side = 1, line = 0.5)
  if(ylab) {mtext(paste0('Days to trafficability after flooding'), side=2, line=2.25)}
  text(x=1.5, y=lab_ht, labels=month)
  dev.off()
}
vioplot_traffic_fine(fname = 'cell_181807_Jan.csv', ylim_vioplot = c(0,95), mar = c(0.1, 3.5, 0.1, 0.1), fig_height = 3.25, fig_width=2, month='January', savename = 'cell_181807_Jan_fine.tif', axis_line = -2, lab_ht = 95, ylab=TRUE)
vioplot_traffic_fine(fname = 'cell_181807_Feb.csv', ylim_vioplot = c(0,95), mar = c(0.1, 3.5, 0.1, 0.1), fig_height = 3.25, fig_width=2, month='February', savename = 'cell_181807_Feb_fine.tif', axis_line = -2, lab_ht = 95, ylab=FALSE)
vioplot_traffic_fine(fname = 'cell_181807_Mar.csv', ylim_vioplot = c(0,95), mar = c(0.1, 3.5, 0.1, 0.1), fig_height = 3.25, fig_width=2, month='March', savename = 'cell_181807_Mar_fine.tif', axis_line = -2, lab_ht = 95, ylab = TRUE)
vioplot_traffic_fine(fname = 'cell_181807_Apr.csv', ylim_vioplot = c(0,95), mar = c(0.1, 3.5, 0.1, 0.1), fig_height = 3.25, fig_width=2, month='April', savename = 'cell_181807_Apr_fine.tif', axis_line = -2, lab_ht = 95, ylab = FALSE)

#texture specific violin plot with 4 months
vioplot_traffic_texture <- function(cellname, texture, to_plot, ylim_vioplot, mar, fig_label, fig_height, fig_width, savename, axis_line,  lab_ht, ylab) {
  df_Jan <- read.csv(file.path(MonthlyMeansDir, paste0(cellname, '_Jan.csv')), stringsAsFactors = FALSE)
  df_Jan <- df_Jan[df_Jan$textural_class == texture,]
  df_Feb <- read.csv(file.path(MonthlyMeansDir, paste0(cellname, '_Feb.csv')), stringsAsFactors = FALSE)
  df_Feb <- df_Feb[df_Feb$textural_class == texture,]
  df_Mar <- read.csv(file.path(MonthlyMeansDir, paste0(cellname, '_Mar.csv')), stringsAsFactors = FALSE)
  df_Mar <- df_Mar[df_Mar$textural_class == texture,]
  df_Apr <- read.csv(file.path(MonthlyMeansDir, paste0(cellname, '_Apr.csv')), stringsAsFactors = FALSE)
  df_Apr <- df_Apr[df_Apr$textural_class == texture,]
  tiff(file = file.path(FiguresDir, 'violin plots', savename), family = 'Times New Roman', width = fig_width, height = fig_height, pointsize = 12, units = 'in', res=800, compression='lzw')
  par(mar=mar)
  vioplot(df_Jan[,4][df_Jan$textural_class==texture], df_Feb[,4][df_Feb$textural_class==texture], df_Mar[,4][df_Mar$textural_class==texture], df_Apr[,4][df_Apr$textural_class==texture], col=rgb(texture_colors$red/255, texture_colors$green/255, texture_colors$blue/255)[texture_colors$textures==texture], wex=1, cex=0.3,  rectCol = 'gray', pchMed = 1, colMed = 'black', ylim = ylim_vioplot, ylab = NULL, xaxt = 'n')
  axis(side = 1, at=1:4, c('Jan', 'Feb', 'Mar', 'Apr'), line=axis_line, tick = FALSE, cex.axis=1)
  # mtext('Textural classes, <40% clay (0-10 cm)', side = 1, line = 0.5)
  if(ylab) {mtext(paste0('Days to trafficability after flooding'), side=2, line=2.25)}
  text(x=2.5, y=lab_ht, labels=texture)
  dev.off()
}
vioplot_traffic_texture(cellname = 'cell_181807', texture = 'sandy loam', ylim_vioplot = c(0,42), mar = c(0.1, 3.5, 0.1, 0.1), fig_height = 3, fig_width=3, savename = 'cell_181807_sandy_loam.tif', axis_line = -2, lab_ht = 42, ylab = TRUE)
vioplot_traffic_texture(cellname = 'cell_181807', texture = 'loam', ylim_vioplot = c(0,60), mar = c(0.1, 3.5, 0.1, 0.1), fig_height = 3, fig_width=3, savename = 'cell_181807_loam.tif', axis_line = -2, lab_ht = 60, ylab = TRUE)
vioplot_traffic_texture(cellname = 'cell_181807', texture = 'loamy sand', ylim_vioplot = c(0,20), mar = c(0.1, 3.5, 0.1, 0.1), fig_height = 3, fig_width=3, savename = 'cell_181807_loamy_sand.tif', axis_line = -2, lab_ht = 20, ylab = TRUE)
vioplot_traffic_texture(cellname = 'cell_181807', texture = 'sand', ylim_vioplot = c(0,22), mar = c(0.1, 3.5, 0.1, 0.1), fig_height = 3, fig_width=3, savename = 'cell_181807_sand.tif', axis_line = -2, lab_ht = 22, ylab = TRUE)
vioplot_traffic_texture(cellname = 'cell_181807', texture = 'clay', ylim_vioplot = c(0,90), mar = c(0.1, 3.5, 0.1, 0.1), fig_height = 3, fig_width=3, savename = 'cell_181807_clay.tif', axis_line = -2, lab_ht = 90, ylab = TRUE)
vioplot_traffic_texture(cellname = 'cell_181807', texture = 'silty clay', ylim_vioplot = c(0,95), mar = c(0.1, 3.5, 0.1, 0.1), fig_height = 3, fig_width=3, savename = 'cell_181807_silty_clay.tif', axis_line = -2, lab_ht = 95, ylab = TRUE)
vioplot_traffic_texture(cellname = 'cell_181807', texture = 'clay loam', ylim_vioplot = c(0,50), mar = c(0.1, 3.5, 0.1, 0.1), fig_height = 3, fig_width=3, savename = 'cell_181807_clay_loam.tif', axis_line = -2, lab_ht = 50, ylab = TRUE)
vioplot_traffic_texture(cellname = 'cell_181807', texture = 'silty clay loam', ylim_vioplot = c(0,70), mar = c(0.1, 3.5, 0.1, 0.1), fig_height = 3, fig_width=3, savename = 'cell_181807_silty_clay_loam.tif', axis_line = -2, lab_ht = 70, ylab = TRUE)
vioplot_traffic_texture(cellname = 'cell_181807', texture = 'silt loam', ylim_vioplot = c(0,70), mar = c(0.1, 3.5, 0.1, 0.1), fig_height = 3, fig_width=3, savename = 'cell_181807_silt_loam.tif', axis_line = -2, lab_ht = 70, ylab = TRUE)
vioplot_traffic_texture(cellname = 'cell_181807', texture = 'sandy clay loam', ylim_vioplot = c(0,27), mar = c(0.1, 3.5, 0.1, 0.1), fig_height = 3, fig_width=3, savename = 'cell_181807_sandy_clay_loam.tif', axis_line = -2, lab_ht = 27, ylab = TRUE)
