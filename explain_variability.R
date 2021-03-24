laptop <- TRUE
if (laptop) {
  dataDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/trafficability study/'
  workDir <- 'C:/Users/smdevine/Desktop/post doc/Dahlke/data from Stathis'
} else {dataDir <- 'C:/Users/smdevine/Desktop/PostDoc/trafficability'}

soil_df <- read.csv(file.path(dataDir, 'soils_of_interest', 'soils_modeled_revised_QCpass_Oct2020.csv'), stringsAsFactors = FALSE)
dim(soil_df)
colnames(soil_df)
list.files(file.path(dataDir, 'climate_run_summaries', 'by_cell'))
cell_130212_results <- read.csv(file.path(dataDir, 'climate_run_summaries', 'by_cell', 'cell_130212_results.csv'), stringsAsFactors = FALSE)
dim(cell_130212_results)
cell_130212_results <- cell_130212_results[cell_130212_results$textural_class!='sandy clay',]
cell_130212_results$textural_class <- NULL
soils_130212_results <- merge(soil_df, cell_130212_results, by='cokey')

#read in soil profile data
#read in soil data
mod_database <- read.csv(file.path(workDir, 'modelling_database.csv'), stringsAsFactors = FALSE, na.strings = '-9.9')
mod_database <- mod_database[which(mod_database$hzn_top<200),] #delete horizons that start at or below 200 cm depth
mod_database <- mod_database[mod_database$cokey %in% soil_df$cokey,]
length(unique(mod_database$cokey)) #2911
mod_database_hzn_des <- mod_database[!mod_database$hzn_desgn %in% c('H', 'H1', 'H2', 'H3', 'H4', 'H5', 'H6'),] 
length(unique(mod_database_hzn_des$cokey)) #only leaves 981
unique(mod_database_hzn_des$hzn_desgn)[grepl('t', unique(mod_database$hzn_desgn))] #"Ct"          "Crt"

soils_130212_results_hzn_desgn <- soils_130212_results[soils_130212_results$cokey %in% mod_database_hzn_des$cokey,]
dim(soils_130212_results_hzn_desgn)
table(soils_130212_results_hzn_desgn$textural_class)

#check results by date and textural class
tapply(soils_130212_results$fd_2009.01.15, soils_130212_results$textural_class, summary)
tapply(soils_130212_results$fd_2009.02.15, soils_130212_results$textural_class, summary)
tapply(soils_130212_results$fd_2009.03.15, soils_130212_results$textural_class, summary)
tapply(soils_130212_results$fd_2009.04.15, soils_130212_results$textural_class, summary)

#textural classes (ordered) for parsing results
textural_classes <- unique(soil_df$textural_class)
textural_classes <- textural_classes[textural_classes != 'sandy clay']
textural_classes <- textural_classes[c(6,3,2,7,5,8,1,9,10,4)]
textural_classes

#summary table of 0-10 cm soil properties by textural class
soil_properties_by_texture <- data.frame(texture=textural_classes, sand=NA, silt=NA, clay=NA, BD=NA, ksat=NA, theta_r=NA, theta_S=NA, alpha=NA, n=NA, theta_wp=NA, theta_fc=NA, h_fc=NA, stringsAsFactors = FALSE)
for (i in seq_along(textural_classes)) {
  soil_properties_by_texture[i,2:ncol(soil_properties_by_texture)] <- c(median(soil_df$sand_10cm[soil_df$textural_class==textural_classes[i]]), median(soil_df$silt_10cm[soil_df$textural_class==textural_classes[i]]), median(soil_df$clay_10cm[soil_df$textural_class==textural_classes[i]]), median(soil_df$bd_13b_10cm[soil_df$textural_class==textural_classes[i]], na.rm = TRUE), median(soil_df$ksat_10cm[soil_df$textural_class==textural_classes[i]]), median(soil_df$theta_r_10cm[soil_df$textural_class==textural_classes[i]]), median(soil_df$theta_s_10cm[soil_df$textural_class==textural_classes[i]]), median(soil_df$alpha_10cm[soil_df$textural_class==textural_classes[i]]), median(soil_df$n_10cm[soil_df$textural_class==textural_classes[i]]),  median(soil_df$theta_wp_10cm[soil_df$textural_class==textural_classes[i]]), median(soil_df$theta_fc_10cm[soil_df$textural_class==textural_classes[i]]), median(soil_df$h_fc_10cm[soil_df$textural_class==textural_classes[i]]))
}
write.csv(soil_properties_by_texture, file.path(dataDir, 'soils_of_interest', 'median_0_10cm_soil_props_by_texture.csv'), row.names = FALSE)

lapply(soils_130212_results[,33:44], function(x) summary(lm(x ~ soils_130212_results$textural_class)))
lapply(soils_130212_results[,33:44], function(x) summary(lm(x ~ soils_130212_results$clay_10cm + soils_130212_results$silt_10cm)))

#exploratory analysis function
exploratory_analysis <- function(df, varname, result) {
  for (i in seq_along(textural_classes)) {
    plot(df[[varname]][df$textural_class==textural_classes[i]], df[[result]][df$textural_class==textural_classes[i]], main=textural_classes[i])
    print(textural_classes[i])
    print(summary(lm(df[[result]][df$textural_class==textural_classes[i]] ~ df[[varname]][df$textural_class==textural_classes[i]])))
  }
}
exploratory_analysis(df = soils_130212_results, varname = 'bd_13b_10cm', result = 'fd_2009.02.15')
exploratory_analysis(df = soils_130212_results, varname = 'ksat_10cm', result = 'fd_2009.02.15')
exploratory_analysis(df = soils_130212_results, varname = 'theta_s_10cm', result = 'fd_2009.02.15')
exploratory_analysis(df = soils_130212_results, varname = 'theta_r_10cm', result = 'fd_2009.02.15')
exploratory_analysis(df = soils_130212_results, varname = 'alpha_10cm', result = 'fd_2009.02.15')
plot(soils_130212_results$alpha_10cm, soils_130212_results$fd_2009.02.15)
plot(soils_130212_results$alpha_10cm, soils_130212_results$sand_10cm)
plot(soils_130212_results$alpha_10cm, soils_130212_results$theta_0.33b_10cm)
exploratory_analysis(df = soils_130212_results, varname = 'alpha_10cm', result = 'theta_0.33b_10cm')
exploratory_analysis(df = soils_130212_results, varname = 'alpha_10cm', result = 'h_fc_10cm')
exploratory_analysis(df = soils_130212_results, varname = 'h_fc_10cm', result = 'theta_0.33b_10cm')
exploratory_analysis(df = soils_130212_results, varname = 'n_10cm', result = 'fd_2009.02.15')
exploratory_analysis(df = soils_130212_results, varname = 'h_fc_10cm', result = 'fd_2009.03.15')
exploratory_analysis(df = soils_130212_results, varname = 'h_fc_10cm', result = 'fd_2009.02.15')
exploratory_analysis(df = soils_130212_results, varname = 'h_fc_10cm', result = 'fd_2009.01.15')
exploratory_analysis(df = soils_130212_results, varname = 'theta_fc_10cm', result = 'fd_2009.02.15')
exploratory_analysis(df = soils_130212_results, varname = 'theta_wp_10cm', result = 'fd_2009.02.15')
exploratory_analysis(df = soils_130212_results, varname = 'theta_0.33b_10cm', result = 'fd_2009.02.15')
exploratory_analysis(df = soils_130212_results, varname = 'theta_15b_10cm', result = 'fd_2009.02.15')
exploratory_analysis(df = soils_130212_results, varname = 'clay_10cm', result = 'fd_2009.02.15')
exploratory_analysis(df = soils_130212_results, varname = 'silt_10cm', result = 'fd_2009.02.15')
exploratory_analysis(df = soils_130212_results, varname = 'sand_10cm', result = 'fd_2009.02.15')
exploratory_analysis(df = soils_130212_results, varname = 'awc_10cm', result = 'fd_2009.02.15')
exploratory_analysis(df = soils_130212_results, varname = 'awc_v2_10cm', result = 'fd_2009.02.15')

#get 25th, 50th, and 75th results' percentiles soil cokeys in each textural class across dates for further investigation of water retention function and fluxes across different depths
colnames(soils_130212_results)
median_soils <- data.frame(texture=textural_classes, fd_2005.01.15=NA, fd_2005.02.15=NA, fd_2005.03.15=NA, fd_2005.04.15=NA, fd_2009.01.15=NA, fd_2009.02.15=NA, fd_2009.03.15=NA, fd_2009.04.15=NA, fd_2015.01.15=NA, fd_2015.02.15=NA, fd_2015.03.15=NA, fd_2015.04.15=NA, stringsAsFactors = FALSE)
Q15_soils <- data.frame(texture=textural_classes, fd_2005.01.15=NA, fd_2005.02.15=NA, fd_2005.03.15=NA, fd_2005.04.15=NA, fd_2009.01.15=NA, fd_2009.02.15=NA, fd_2009.03.15=NA, fd_2009.04.15=NA, fd_2015.01.15=NA, fd_2015.02.15=NA, fd_2015.03.15=NA, fd_2015.04.15=NA, stringsAsFactors = FALSE)
Q85_soils <- data.frame(texture=textural_classes, fd_2005.01.15=NA, fd_2005.02.15=NA, fd_2005.03.15=NA, fd_2005.04.15=NA, fd_2009.01.15=NA, fd_2009.02.15=NA, fd_2009.03.15=NA, fd_2009.04.15=NA, fd_2015.01.15=NA, fd_2015.02.15=NA, fd_2015.03.15=NA, fd_2015.04.15=NA, stringsAsFactors = FALSE)

identify_soil <- function(input_df) {
  for (i in seq_along(textural_classes)) {
    input_df_texture <- input_df[input_df$textural_class==textural_classes[i],]
    for (j in 2:13) {
      final_df_texture <- input_df_texture[!is.na(input_df_texture[,colnames(median_soils)[j]]),]
      crit_vals <- quantile(final_df_texture[,colnames(median_soils)[j]], probs = c(0.15,0.5,0.85))
      min_Q15 <- min(abs(final_df_texture[,colnames(Q15_soils)[j]] - crit_vals[1]))
      min_median <- min(abs(final_df_texture[,colnames(median_soils)[j]] - crit_vals[2]))
      min_Q85 <- min(abs(final_df_texture[,colnames(Q85_soils)[j]] - crit_vals[3]))
      if(length(final_df_texture$cokey[which(final_df_texture[,colnames(Q15_soils)[j]]==crit_vals[1])])==0) {
        Q15_soils[i,j] <- paste(final_df_texture$cokey[which(abs(final_df_texture[,colnames(Q15_soils)[j]] - crit_vals[1])==min_Q15)], collapse = ',')
      } else {
        Q15_soils[i,j] <- paste(final_df_texture$cokey[which(final_df_texture[,colnames(Q15_soils)[j]]==crit_vals[1])], collapse = ',')
        }
      if (length(final_df_texture$cokey[which(final_df_texture[,colnames(median_soils)[j]]==crit_vals[2])])==0) {
        median_soils[i,j] <- paste(final_df_texture$cokey[which(abs(final_df_texture[,colnames(median_soils)[j]] - crit_vals[2])==min_median)], collapse = ',')
      } else {
        median_soils[i,j] <- paste(final_df_texture$cokey[which(final_df_texture[,colnames(median_soils)[j]]==crit_vals[2])], collapse = ',')
      }
      if (length(final_df_texture$cokey[which(final_df_texture[,colnames(Q85_soils)[j]]==crit_vals[3])])==0) {
        Q85_soils[i,j] <- paste(final_df_texture$cokey[which(abs(final_df_texture[,colnames(Q85_soils)[j]] - crit_vals[3])==min_Q85)], collapse = ',')
      } else {
        Q85_soils[i,j] <- paste(final_df_texture$cokey[which(final_df_texture[,colnames(Q85_soils)[j]]==crit_vals[3])], collapse = ',')
      }
    }
  }
  list(Q15=Q15_soils, median=median_soils, Q85=Q85_soils)
}
soils_130212_quantiles <- identify_soil(soils_130212_results)

soils_130212_quantiles$Q15[which(textural_classes=='silty clay'),]
identify_most_common <- function(texture, df, q) {
  test <- paste(df[[q]][which(textural_classes==texture),2:13], collapse = ',')
  test <- unlist(strsplit(test, ','))
  test <- test[test != '']
  table(test)[order(table(test))]
}

textural_classes
check <- 'clay'
test <- identify_most_common(check, soils_130212_quantiles, 'Q85')
test <- as.data.frame(test)
colnames(test) <- c('cokey', 'freq')
test$cokey <- as.integer(as.character(test$cokey))
test <- test[order(test$freq, decreasing = TRUE),]
head(test, 10)
test[which(test$cokey %in% mod_database_hzn_des$cokey),]
test2 <- head(which(test$cokey %in% mod_database_hzn_des$cokey), 10)
test3 <- test2[1] #2,3,4 ok
test3 <- 1
mod_database_hzn_des[test$cokey[test3] == mod_database_hzn_des$cokey,]
mod_database[test$cokey[test3] == mod_database$cokey,]
cell_130212_results$fd_2009.01.15[cell_130212_results$cokey==test$cokey[test3]]
summary(soils_130212_results$fd_2009.01.15[soils_130212_results$textural_class==check])
cell_130212_results$fd_2009.02.15[cell_130212_results$cokey==test$cokey[test3]]
summary(soils_130212_results$fd_2009.02.15[soils_130212_results$textural_class==check])
cell_130212_results$fd_2009.03.15[cell_130212_results$cokey==test$cokey[test3]]
summary(soils_130212_results$fd_2009.03.15[soils_130212_results$textural_class==check])
cell_130212_results$fd_2009.04.15[cell_130212_results$cokey==test$cokey[test3]]
summary(soils_130212_results$fd_2009.04.15[soils_130212_results$textural_class==check])

#check 85th vs. 15th manually
cokey_15 <- 19909 #15th percentile
cokey_85 <- test$cokey[test3] #85th percentiles
soil_df[soil_df$cokey==cokey_15,]
soil_df[soil_df$cokey==cokey_85,] 
lapply(cell_130212_results[,3:14], function(x) x[cell_130212_results$cokey==cokey_85] - x[cell_130212_results$cokey==cokey_15])
mod_database[cokey_15==mod_database$cokey,]
mod_database[cokey_85==mod_database$cokey,]

textural_classes
lapply(textural_classes, function(x) identify_most_common(x, soils_130212_quantiles, 'Q15'))
lapply(textural_classes, function(x) identify_most_common(x, soils_130212_quantiles, 'median'))
lapply(textural_classes, function(x) identify_most_common(x, soils_130212_quantiles, 'Q85'))


property_by_cokey <- function(cokey, property) {
  soil <- mod_database[mod_database$cokey==cokey,]
  soil <- soil[!is.na(soil[[property]]),]
  profile_depth <- max(soil$hzn_bot) - min(soil$hzn_top)
  soil$hzn_thickness <- soil$hzn_bot - soil$hzn_top
  soil$wtd.factor <- soil$hzn_thickness/profile_depth
  sum(soil$wtd.factor * soil[[property]])
}
mod_database[mod_database$cokey==soils_130212_results$cokey[9],]
property_by_cokey(2691, 'Ks..cm.d.')
property_by_cokey(2691, 'clay')
#check clay
8.6*(25/153) + 31.3*(36/153) + 23.1*(36/153) + 13.9*(55/153) + 9.6*(1/153) #19.26471

t_horizons_by_cokey <- function(cokey) {
  soil <- mod_database[mod_database$cokey==cokey,]
  hzn_desgns <- soil$hzn_desgn
  hzn_desgns <- hzn_desgns[hzn_desgns != 'C1 to C4']
  hzn_desgns <- hzn_desgns[hzn_desgns != 'Ct']
  hzn_desgns <- hzn_desgns[hzn_desgns != 'Crt']
  hzn_desgns <- hzn_desgns[hzn_desgns != 'R2t']
  # print(hzn_desgns)
  sum(grepl('t', hzn_desgns)) > 0
}

n_horizons_by_cokey <- function(cokey) {
  soil <- mod_database[mod_database$cokey==cokey,]
  length(soil$hzn_desgn)
}

surface_thickness_by_cokey <- function(cokey) {
  soil <- mod_database[mod_database$cokey==cokey,]
  soil <- soil[order(soil$hzn_top),]
  soil$hzn_bot[1] - soil$hzn_top[1] 
}
surface_thickness_by_cokey(soils_130212_results$cokey[1])
  
t_horizon_depth_by_cokey <- function(cokey) {
  soil <- mod_database[mod_database$cokey==cokey,]
  hzn_desgns <- soil[,c('hzn_desgn', 'hzn_top', 'hzn_bot')]
  hzn_desgns <- hzn_desgns[hzn_desgns$hzn_desgn != 'C1 to C4',]
  hzn_desgns <- hzn_desgns[hzn_desgns$hzn_desgn != 'Ct',]
  hzn_desgns <- hzn_desgns[hzn_desgns$hzn_desgn != 'Crt',]
  hzn_desgns <- hzn_desgns[hzn_desgns$hzn_desgn != 'R2t',]
  if(sum(grepl('t', hzn_desgns$hzn_desgn)) > 0) {
    hzn_desgns <- hzn_desgns[grepl('t', hzn_desgns$hzn_desgn),]
    min(hzn_desgns$hzn_top)
  } else {NA}
}
t_horizon_depth_by_cokey(11681)

mod_database[mod_database$cokey==soils_130212_results_hzn_desgn$cokey[5],]
t_horizons_by_cokey(soils_130212_results_hzn_desgn$cokey[5])
t_horizon_depth_by_cokey(soils_130212_results_hzn_desgn$cokey[5])
soils_130212_results_hzn_desgn[5,]
soils_130212_results_hzn_desgn$t_horizon <- sapply(soils_130212_results_hzn_desgn$cokey, t_horizons_by_cokey)
soils_130212_results_hzn_desgn$t_horizon_depth <- sapply(soils_130212_results_hzn_desgn$cokey, t_horizon_depth_by_cokey)
tapply(soils_130212_results_hzn_desgn$t_horizon, soils_130212_results_hzn_desgn$textural_class, table)

#compare stats by textural class and presence of Bt
tapply(soils_130212_results_hzn_desgn$fd_2005.02.15[soils_130212_results_hzn_desgn$textural_class=='loamy sand'], soils_130212_results_hzn_desgn$t_horizon[soils_130212_results_hzn_desgn$textural_class=='loamy sand'], summary) #median is 1 day slower; mean is 1.5 days lower

tapply(soils_130212_results_hzn_desgn$fd_2005.02.15[soils_130212_results_hzn_desgn$textural_class=='sandy loam'], soils_130212_results_hzn_desgn$t_horizon[soils_130212_results_hzn_desgn$textural_class=='sandy loam'], summary) #mean is 1 day slower; median 0.5 day slower

tapply(soils_130212_results_hzn_desgn$fd_2005.02.15[soils_130212_results_hzn_desgn$textural_class=='loam'], soils_130212_results_hzn_desgn$t_horizon[soils_130212_results_hzn_desgn$textural_class=='loam'], summary) #practically the same

tapply(soils_130212_results_hzn_desgn$fd_2005.02.15[soils_130212_results_hzn_desgn$textural_class=='silt loam'], soils_130212_results_hzn_desgn$t_horizon[soils_130212_results_hzn_desgn$textural_class=='silt loam'], summary) #mean is 5 days slower; median 11 days slower

tapply(soils_130212_results_hzn_desgn$fd_2005.02.15[soils_130212_results_hzn_desgn$textural_class=='clay loam'], soils_130212_results_hzn_desgn$t_horizon[soils_130212_results_hzn_desgn$textural_class=='clay loam'], summary) #median 0.4 day faster; mean is 0.5 day slower

tapply(soils_130212_results_hzn_desgn$fd_2005.02.15[soils_130212_results_hzn_desgn$textural_class=='silty clay loam'], soils_130212_results_hzn_desgn$t_horizon[soils_130212_results_hzn_desgn$textural_class=='silty clay loam'], summary) #median is 5 days faster; mean is 7 days faster

tapply(soils_130212_results_hzn_desgn$fd_2005.02.15[soils_130212_results_hzn_desgn$textural_class=='sandy clay loam'], soils_130212_results_hzn_desgn$t_horizon[soils_130212_results_hzn_desgn$textural_class=='sandy clay loam'], summary)

tapply(soils_130212_results_hzn_desgn$fd_2005.02.15[soils_130212_results_hzn_desgn$textural_class=='clay'], soils_130212_results_hzn_desgn$t_horizon[soils_130212_results_hzn_desgn$textural_class=='clay'], summary)

tapply(soils_130212_results_hzn_desgn$fd_2005.02.15[soils_130212_results_hzn_desgn$textural_class=='sand'], soils_130212_results_hzn_desgn$t_horizon[soils_130212_results_hzn_desgn$textural_class=='sand'], summary)

#now do same but with loams but with reference to where Bt horizon starts
hist(soils_130212_results_hzn_desgn$t_horizon_depth)
hist(soils_130212_results_hzn_desgn$t_horizon_depth[soils_130212_results_hzn_desgn$textural_class=='loam'])
summary(soils_130212_results_hzn_desgn$fd_2009.02.15[which(soils_130212_results_hzn_desgn$textural_class=='loam' & soils_130212_results_hzn_desgn$t_horizon_depth<=10)])
summary(soils_130212_results_hzn_desgn$fd_2009.02.15[which(soils_130212_results_hzn_desgn$textural_class=='loam' & soils_130212_results_hzn_desgn$t_horizon_depth>10 & soils_130212_results_hzn_desgn$t_horizon_depth<=30)])
summary(soils_130212_results_hzn_desgn$fd_2009.02.15[which(soils_130212_results_hzn_desgn$textural_class=='loam' & soils_130212_results_hzn_desgn$t_horizon_depth>30)])
summary(soils_130212_results_hzn_desgn$fd_2009.02.15[which(soils_130212_results_hzn_desgn$textural_class=='loam' & is.na(soils_130212_results_hzn_desgn$t_horizon_depth))])

length(which(soils_130212_results_hzn_desgn$textural_class=='loam' & soils_130212_results_hzn_desgn$t_horizon_depth<=10)) #27
length(which(soils_130212_results_hzn_desgn$textural_class=='loam' & soils_130212_results_hzn_desgn$t_horizon_depth>10&soils_130212_results_hzn_desgn$t_horizon_depth<=30)) #71
length(which(soils_130212_results_hzn_desgn$textural_class=='loam' & soils_130212_results_hzn_desgn$t_horizon_depth>30))#51
length(which(soils_130212_results_hzn_desgn$textural_class=='loam' & is.na(soils_130212_results_hzn_desgn$t_horizon_depth))) #135
    
#and with sandy loam
hist(soils_130212_results_hzn_desgn$t_horizon_depth[soils_130212_results_hzn_desgn$textural_class=='sandy loam'])
summary(soils_130212_results_hzn_desgn$fd_2009.02.15[which(soils_130212_results_hzn_desgn$textural_class=='sandy loam' & soils_130212_results_hzn_desgn$t_horizon_depth<=10)])
summary(soils_130212_results_hzn_desgn$fd_2009.02.15[which(soils_130212_results_hzn_desgn$textural_class=='sandy loam' & soils_130212_results_hzn_desgn$t_horizon_depth>10 & soils_130212_results_hzn_desgn$t_horizon_depth<=30)])
summary(soils_130212_results_hzn_desgn$fd_2009.02.15[which(soils_130212_results_hzn_desgn$textural_class=='sandy loam' & soils_130212_results_hzn_desgn$t_horizon_depth>30)])
summary(soils_130212_results_hzn_desgn$fd_2009.02.15[which(soils_130212_results_hzn_desgn$textural_class=='sandy loam' & is.na(soils_130212_results_hzn_desgn$t_horizon_depth))])

length(which(soils_130212_results_hzn_desgn$textural_class=='sandy loam' & soils_130212_results_hzn_desgn$t_horizon_depth<=10)) #10
length(which(soils_130212_results_hzn_desgn$textural_class=='sandy loam' & soils_130212_results_hzn_desgn$t_horizon_depth>10&soils_130212_results_hzn_desgn$t_horizon_depth<=30)) #56
length(which(soils_130212_results_hzn_desgn$textural_class=='sandy loam' & soils_130212_results_hzn_desgn$t_horizon_depth>30))#68
length(which(soils_130212_results_hzn_desgn$textural_class=='sandy loam' & is.na(soils_130212_results_hzn_desgn$t_horizon_depth))) #198

#and clay loam
hist(soils_130212_results_hzn_desgn$t_horizon_depth[soils_130212_results_hzn_desgn$textural_class=='clay loam'])
summary(soils_130212_results_hzn_desgn$fd_2009.02.15[which(soils_130212_results_hzn_desgn$textural_class=='clay loam' & soils_130212_results_hzn_desgn$t_horizon_depth<=10)])
summary(soils_130212_results_hzn_desgn$fd_2009.02.15[which(soils_130212_results_hzn_desgn$textural_class=='clay loam' & soils_130212_results_hzn_desgn$t_horizon_depth>10 & soils_130212_results_hzn_desgn$t_horizon_depth<=30)])
summary(soils_130212_results_hzn_desgn$fd_2009.02.15[which(soils_130212_results_hzn_desgn$textural_class=='clay loam' & soils_130212_results_hzn_desgn$t_horizon_depth>30)])
summary(soils_130212_results_hzn_desgn$fd_2009.02.15[which(soils_130212_results_hzn_desgn$textural_class=='clay loam' & is.na(soils_130212_results_hzn_desgn$t_horizon_depth))])

length(which(soils_130212_results_hzn_desgn$textural_class=='clay loam' & soils_130212_results_hzn_desgn$t_horizon_depth<=10)) #5
length(which(soils_130212_results_hzn_desgn$textural_class=='clay loam' & soils_130212_results_hzn_desgn$t_horizon_depth>10&soils_130212_results_hzn_desgn$t_horizon_depth<=30)) #12
length(which(soils_130212_results_hzn_desgn$textural_class=='clay loam' & soils_130212_results_hzn_desgn$t_horizon_depth>30))#5
length(which(soils_130212_results_hzn_desgn$textural_class=='clay loam' & is.na(soils_130212_results_hzn_desgn$t_horizon_depth))) #45

#and silt loam
hist(soils_130212_results_hzn_desgn$t_horizon_depth[soils_130212_results_hzn_desgn$textural_class=='silt loam'])
summary(soils_130212_results_hzn_desgn$fd_2009.02.15[which(soils_130212_results_hzn_desgn$textural_class=='silt loam' & soils_130212_results_hzn_desgn$t_horizon_depth<=10)])
summary(soils_130212_results_hzn_desgn$fd_2009.02.15[which(soils_130212_results_hzn_desgn$textural_class=='silt loam' & soils_130212_results_hzn_desgn$t_horizon_depth>10 & soils_130212_results_hzn_desgn$t_horizon_depth<=30)])
summary(soils_130212_results_hzn_desgn$fd_2009.02.15[which(soils_130212_results_hzn_desgn$textural_class=='silt loam' & soils_130212_results_hzn_desgn$t_horizon_depth>30)])
summary(soils_130212_results_hzn_desgn$fd_2009.02.15[which(soils_130212_results_hzn_desgn$textural_class=='silt loam' & is.na(soils_130212_results_hzn_desgn$t_horizon_depth))])

length(which(soils_130212_results_hzn_desgn$textural_class=='silt loam' & soils_130212_results_hzn_desgn$t_horizon_depth<=10)) #3
length(which(soils_130212_results_hzn_desgn$textural_class=='silt loam' & soils_130212_results_hzn_desgn$t_horizon_depth>10&soils_130212_results_hzn_desgn$t_horizon_depth<=30)) #3
length(which(soils_130212_results_hzn_desgn$textural_class=='silt loam' & soils_130212_results_hzn_desgn$t_horizon_depth>30))#9
length(which(soils_130212_results_hzn_desgn$textural_class=='silt loam' & is.na(soils_130212_results_hzn_desgn$t_horizon_depth))) #33

#compare stats with profile weighted ksat and clay
soils_130212_results$clay_profile <- sapply(soils_130212_results$cokey, property_by_cokey, property='clay')
soils_130212_results$silt_profile <- sapply(soils_130212_results$cokey, property_by_cokey, property='silt')
soils_130212_results$sand_profile <- sapply(soils_130212_results$cokey, property_by_cokey, property='sand')
soils_130212_results$ksat_profile <- sapply(soils_130212_results$cokey, property_by_cokey, property='Ks..cm.d.')
soils_130212_results$theta_s_profile <- sapply(soils_130212_results$cokey, property_by_cokey, property='teta_s')
soils_130212_results$theta_r_profile <- sapply(soils_130212_results$cokey, property_by_cokey, property='teta_r')
soils_130212_results$alpha_profile <- sapply(soils_130212_results$cokey, property_by_cokey, property='alpha..1.cm.')
soils_130212_results$n_horizons <- sapply(soils_130212_results$cokey, n_horizons_by_cokey)
soils_130212_results$surface_thickness <- sapply(soils_130212_results$cokey, surface_thickness_by_cokey)
summary(soils_130212_results$clay_profile)
summary(soils_130212_results$ksat_profile)
tapply(soils_130212_results$n_horizons, soils_130212_results$textural_class, table)

#clay
plot(soils_130212_results$clay_10cm, soils_130212_results$clay_profile)
plot(soils_130212_results$clay_10cm, soils_130212_results$fd_2009.02.15)
plot(soils_130212_results$clay_profile, soils_130212_results$fd_2009.02.15)

#sand
plot(soils_130212_results$sand_10cm, soils_130212_results$sand_profile)
plot(soils_130212_results$sand_10cm, soils_130212_results$fd_2009.02.15)
plot(soils_130212_results$sand_profile, soils_130212_results$fd_2009.02.15)

#ksat
plot(log(soils_130212_results$ksat_10cm), log(soils_130212_results$ksat_profile))
plot(log(soils_130212_results$ksat_10cm), soils_130212_results$fd_2009.02.15)
plot(log(soils_130212_results$ksat_profile), soils_130212_results$fd_2009.02.15)

#linear models by textural class
exploratory_analysis(df = soils_130212_results, varname = 'n_horizons', result = 'fd_2009.02.15')
exploratory_analysis(df = soils_130212_results, varname = 'surface_thickness', result = 'fd_2009.02.15')
exploratory_analysis(df = soils_130212_results, varname = 'clay_profile', result = 'fd_2009.02.15')
exploratory_analysis(df = soils_130212_results, varname = 'sand_profile', result = 'fd_2009.02.15')
exploratory_analysis(df = soils_130212_results, varname = 'ksat_profile', result = 'fd_2009.02.15')
exploratory_analysis(df = soils_130212_results, varname = 'theta_s_profile', result = 'fd_2009.02.15')
exploratory_analysis(df = soils_130212_results, varname = 'theta_r_profile', result = 'fd_2009.02.15')
exploratory_analysis(df = soils_130212_results, varname = 'alpha_profile', result = 'fd_2009.02.15')

#additional overall relationships
plot(soils_130212_results$awc_v2_10cm, soils_130212_results$fd_2009.02.15)
plot(soils_130212_results$awc_10cm, soils_130212_results$fd_2009.02.15)
