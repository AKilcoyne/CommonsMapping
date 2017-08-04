#######################################################################################################
### Automated object-based classification using RF for commons mapping ##############################
#
# MASKED NORMAL SEGMENTS


library(rgdal)
library(raster)
library(rgeos)
library(randomForest)
library(impute) # install zip file from: www.bioconductor.org/packages/release/bioc/html/impute.html

setwd("C:/Users/Bertie/Documents/RPA_Commons_Mapping/CommonsMapping")

source("Scripts/zonal_stats.r")
source("Scripts/user_producer_accuracy.r")

#Training data
#training.data.habitat.shp <- readOGR("Training_Data/FieldData_BNG_edits.shp", "FieldData_BNG_edits")

#Edited training data
training.data.habitat.shp <- readOGR("Training_Data/FieldData_BNG_edits_bsh_20170801_largesegments.shp", "FieldData_BNG_edits_bsh_20170801_largesegments")

#Small test area
#segmentation.shp <- readOGR("Segmentation/results/Test_RPA_COMMONS_MAPPING_SEGMENTS.shp", "Test_RPA_COMMONS_MAPPING_SEGMENTS", useC=T)

#Commons only
#segmentation.shp <- readOGR("Segmentation/results/RPA_COMMONS_MAPPING_SEGMENTS_COMMONSONLY.shp", "RPA_COMMONS_MAPPING_SEGMENTS_COMMONSONLY", useC=T)

#normal
segmentation.shp <- readOGR("Segmentation/results/RPA_COMMONS_MAPPING_SEGMENTS.shp", "RPA_COMMONS_MAPPING_SEGMENTS", useC=T)

#LARGE
#segmentation.shp <- readOGR("Segmentation/results/RPA_COMMONS_MAPPING_SEGMENTS_LARGE.shp", "RPA_COMMONS_MAPPING_SEGMENTS_LARGE", useC=T)


### Zonal statistics layers ###########################################################################
#
#

S2_winter <- "S2/s2_20151229_137_mask_uwm.tif"
WV2_summer <- "WV2/SHER_8BAND_MS_20150409_mask_uwm.tif"
SPOT_summer <- "SPOT/SPOT_SHER_4BAND_MS_20150611_mask_uwm.tif"
height <- "Height/IHM_05m_AOI_BNG_uwm.tif"
slope <- "Height/IHM_05m_AOI_BNG_slope_uwm.tif"
aspect <- "Height/IHM_05m_AOI_BNG_aspect_uwm.tif"
sar_summer <- "S1/S1A_IW_GRDH_1SDV_20150503T061321_20150503T061526_005753_007632_FD40_AOI_uwm.tif"
sar_winter <- "S1/S1A_IW_GRDH_1SDV_20151229T061338_20151229T061518_009253_00D59B_F24A_AOI_uwm.tif"
LSU_summer <- "LSU/SHER_8BAND_MS_20150409_mask_unmixed_uwm.img"
LSU_winter <- "LSU/s2_20151229_137_mask_unmixed_uwm.tif"
OS_dist_building <- "OS/Proximity_Building_uwm.tif"
OS_dist_road <- "OS/Proximity_Road_uwm.tif"
OS_dist_surfacewaterline <- "OS/Proximity_SurfaceWaterLine_uwm.tif"
OS_dist_surfacewaterarea <- "OS/Proximity_SurfaceWaterArea_uwm.tif"
OS_dist_woodland <- "OS/Proximity_Woodland_uwm.tif"
bioclim_max_temp <- "Bioclim/bio2_30s_05_maxtemp_BNG27700_uwm.tif"
bioclim_min_temp <- "Bioclim/bio2_30s_06_mintemp_BNG27700_uwm.tif"
bioclim_annual_rainfall <- "Bioclim/bio2_30s_12_annaulrainfall_BNG27700_uwm.tif"
moorlandline_dist <- "MoorlandLine/Proximity_MoorlandLineRaster_uwm.tif"

list.rasters <- list(S2_winter_blue=c(S2_winter, 1),
                     S2_winter_green=c(S2_winter, 2),
                     S2_winter_red=c(S2_winter, 3),
                     S2_winter_rededge5=c(S2_winter, 4),
                     S2_winter_rededge6=c(S2_winter, 5),
                     S2_winter_rededge7=c(S2_winter, 6),
                     S2_winter_rededge8a=c(S2_winter, 7),
                     S2_winter_nir=c(S2_winter, 8),
                     S2_winter_swir1=c(S2_winter, 9),
                     S2_winter_swir2=c(S2_winter, 10),
                     S2_winter_blue_median=c(S2_winter, 1, "median"),
                     S2_winter_green_median=c(S2_winter, 2, "median"),
                     S2_winter_red_median=c(S2_winter, 3, "median"),
                     S2_winter_rededge5_median=c(S2_winter, 4, "median"),
                     S2_winter_rededge6_median=c(S2_winter, 5, "median"),
                     S2_winter_rededge7_median=c(S2_winter, 6, "median"),
                     S2_winter_rededge8a_median=c(S2_winter, 7, "median"),
                     S2_winter_nir_median=c(S2_winter, 8, "median"),
                     S2_winter_swir1_median=c(S2_winter, 9, "median"),
                     S2_winter_swir2_median=c(S2_winter, 10, "median"),
                     S2_winter_blue_sd=c(S2_winter, 1, "sd"),
                     S2_winter_green_sd=c(S2_winter, 2, "sd"),
                     S2_winter_red_sd=c(S2_winter, 3, "sd"),
                     S2_winter_rededge5_sd=c(S2_winter, 4, "sd"),
                     S2_winter_rededge6_sd=c(S2_winter, 5, "sd"),
                     S2_winter_rededge7_sd=c(S2_winter, 6, "sd"),
                     S2_winter_rededge8a_sd=c(S2_winter, 7, "sd"),
                     S2_winter_nir_sd=c(S2_winter, 8, "sd"),
                     S2_winter_swir1_sd=c(S2_winter, 9, "sd"),
                     S2_winter_swir2_sd=c(S2_winter, 10, "sd"),
                     
                     WV2_summer_coastal=c(WV2_summer, 1),
                     WV2_summer_blue=c(WV2_summer, 2),
                     WV2_summer_green=c(WV2_summer, 3),
                     WV2_summer_yellow=c(WV2_summer, 4),
                     WV2_summer_red=c(WV2_summer, 5),
                     WV2_summer_rededge=c(WV2_summer, 6),
                     WV2_summer_nir1=c(WV2_summer, 7),
                     WV2_summer_nir2=c(WV2_summer, 8),
                     WV2_summer_coastal_median=c(WV2_summer, 1, "median"),
                     WV2_summer_blue_median=c(WV2_summer, 2, "median"),
                     WV2_summer_green_median=c(WV2_summer, 3, "median"),
                     WV2_summer_yellow_median=c(WV2_summer, 4, "median"),
                     WV2_summer_red_median=c(WV2_summer, 5, "median"),
                     WV2_summer_rededge_median=c(WV2_summer, 6, "median"),
                     WV2_summer_nir1_median=c(WV2_summer, 7, "median"),
                     WV2_summer_nir2_median=c(WV2_summer, 8, "median"),
                     WV2_summer_coastal_sd=c(WV2_summer, 1, "sd"),
                     WV2_summer_blue_sd=c(WV2_summer, 2, "sd"),
                     WV2_summer_green_sd=c(WV2_summer, 3, "sd"),
                     WV2_summer_yellow_sd=c(WV2_summer, 4, "sd"),
                     WV2_summer_red_sd=c(WV2_summer, 5, "sd"),
                     WV2_summer_rededge_sd=c(WV2_summer, 6, "sd"),
                     WV2_summer_nir1_sd=c(WV2_summer, 7, "sd"),
                     WV2_summer_nir2_sd=c(WV2_summer, 8, "sd"),
                     
                     SPOT_summer_blue=c(SPOT_summer, 1),
                     SPOT_summer_green=c(SPOT_summer, 2),
                     SPOT_summer_red=c(SPOT_summer, 3),
                     SPOT_summer_nir=c(SPOT_summer, 4),
                     SPOT_summer_blue_median=c(SPOT_summer, 1, "median"),
                     SPOT_summer_green_median=c(SPOT_summer, 2, "median"),
                     SPOT_summer_red_median=c(SPOT_summer, 3, "median"),
                     SPOT_summer_nir_median=c(SPOT_summer, 4, "median"),
                     SPOT_summer_blue_sd=c(SPOT_summer, 1, "sd"),
                     SPOT_summer_green_sd=c(SPOT_summer, 2, "sd"),
                     SPOT_summer_red_sd=c(SPOT_summer, 3, "sd"),
                     SPOT_summer_nir_sd=c(SPOT_summer, 4, "sd"),
                     
                     sar_summer_vh=c(sar_summer,1),
                     sar_summer_vv=c(sar_summer,2),
                     sar_summer_vh_median=c(sar_summer,1, "median"),
                     sar_summer_vv_median=c(sar_summer,2, "median"),
                     sar_summer_vh_sd=c(sar_summer,1, "sd"),
                     sar_summer_vv_sd=c(sar_summer,2, "sd"),
                     
                     sar_winter_vh=c(sar_winter,1),
                     sar_winter_vv=c(sar_winter,2),
                     sar_winter_vh_median=c(sar_winter,1, "median"),
                     sar_winter_vv_median=c(sar_winter,2, "median"),
                     sar_winter_vh_sd=c(sar_winter,1, "sd"),
                     sar_winter_vv_sd=c(sar_winter,2, "sd"),
                     
                     LSU_summer_PV=c(LSU_summer,1),
                     LSU_summer_NPV=c(LSU_summer,2),
                     LSU_summer_S=c(LSU_summer,3),
                     LSU_winter_PV=c(LSU_winter,1),
                     LSU_winter_NPV=c(LSU_winter,2),
                     LSU_winter_S=c(LSU_winter,3),
                     
                     dist_building=c(OS_dist_building,1),
                     dist_road=c(OS_dist_road,1),
                     dist_surfacewater_l=c(OS_dist_surfacewaterline,1),
                     dist_surfacewater_a=c(OS_dist_surfacewaterarea,1),
                     dist_woodland=c(OS_dist_woodland,1),
                     
                     max_temp=c(bioclim_max_temp,1),
                     min_temp=c(bioclim_min_temp,1),
                     annual_rainfall=c(bioclim_annual_rainfall,1),
                     
                     height=c(height,1),
                     slope=c(slope,1),
                     aspect=c(aspect, 1),
                     
                     moorlandline_dist=c(moorlandline_dist,1))



### Zonal Stats for Segmented Polygons ################################################################
#
#
#segmentation.raster <-raster("Segmentation/Segmentation_LARGE_Rasterised.tif")
segmentation.raster <-raster("Segmentation/Segmentation_Rasterised.tif")
#segmentation.raster <-raster("Segmentation/Test_Segmentation_Rasterised.tif")
#segmentation.raster <-raster("Segmentation/Segmentation_Rasterised_Commons.tif")


# Calculate the zonal stats for each segmented polygon.
start <- proc.time()
zonal_stats_seg <- zonal_stats_raster(segmentation.raster, list.rasters, clusters=10, tiles=5)
proc.time()-start

# Save the results as an intermediate file
zonal_stats_seg <- write.table(zonal_stats_seg, "Zonal_stats/zonal_stats_seg_20170801_MASKED_NORMAL.txt", sep="\t")
zonal_stats_seg <- read.table("Zonal_stats/zonal_stats_seg_20170801_MASKED_NORMAL.txt", sep="\t", header=T)

# Append area and perimeter from shapefile if not already calculated
if (!"area_ratio1" %in% names(zonal_stats_seg))
{
  zonal_stats_seg <- merge(zonal_stats_seg, segmentation.shp, by="ID")
  zonal_stats_seg$area_ratio1 <- with(zonal_stats_seg, Area_m2/Perim_m)
  zonal_stats_seg$area_ratio2 <- with(zonal_stats_seg, Perim_m/sqrt(Area_m2))
}

# Ensure that catagorical data doesn't having any missing or inf values
zonal_stats_seg[is.na(zonal_stats_seg)] <- 0
zonal_stats_seg[sapply(zonal_stats_seg, is.infinite)] <- 0

# Impute missing values for all S1 and S2 columns, excluding max and min statistics
#impute.cols <- grepl("S2|sar|LSU|dist_building",colnames(zonal_stats_seg)) & !grepl("max|min",colnames(zonal_stats_seg))
#zonal_stats.imputed <- impute.knn(as.matrix(zonal_stats_seg[,impute.cols]))
#zonal_stats_seg <- cbind(zonal_stats_seg[,!impute.cols], zonal_stats.imputed$data[,colnames(zonal_stats_seg)[impute.cols]])


## Indices ###############

# Calculate NDVI and NDWI
zonal_stats_seg$WV2_summer_ndvi <- with(zonal_stats_seg, (WV2_summer_nir1 - WV2_summer_red)/(WV2_summer_nir1 + WV2_summer_red))

zonal_stats_seg$SPOT_summer_ndvi <- with(zonal_stats_seg, (SPOT_summer_nir - SPOT_summer_red)/(SPOT_summer_nir + SPOT_summer_red))

zonal_stats_seg$S2_winter_ndvi <- with(zonal_stats_seg, (S2_winter_nir - S2_winter_red)/(S2_winter_nir + S2_winter_red))

zonal_stats_seg$S2_winter_ndwi <- with(zonal_stats_seg, (S2_winter_nir - S2_winter_swir1)/(S2_winter_nir + S2_winter_swir1))

# Top indices identified as important - edit 
# zonal_stats_seg$sar_winter_ndvhvvi <- with(zonal_stats_seg, (sar_winter_vh_median - sar_winter_vv_median)/(sar_winter_vh_median + sar_winter_vv_median))
# 
# zonal_stats_seg$s2_summer_greenmrededge5m_di <- with(zonal_stats_seg, (S2_summer_green_median - S2_summer_rededge5_median))
# zonal_stats_seg$s2_summer_bluerededge5_di <- with(zonal_stats_seg, (S2_summer_blue - S2_summer_rededge5))
# zonal_stats_seg$s2_summer_redmrededge5m_di <- with(zonal_stats_seg, (S2_summer_red_median - S2_summer_rededge5_median))
# zonal_stats_seg$s2_summer_rededge6rededge7_ndi <- with(zonal_stats_seg, (S2_summer_rededge6 - S2_summer_rededge7)/(S2_summer_rededge6 + S2_summer_rededge7))
# zonal_stats_seg$s2_summer_redmSWIR2m_di <- with(zonal_stats_seg, (S2_summer_red_median - S2_summer_swir2_median))
# zonal_stats_seg$S2_summer_rededge6_S2_summer_rededge7_ri <- with(zonal_stats_seg, (S2_summer_rededge6 / S2_summer_rededge7))
# zonal_stats_seg$s2_summer_rededge6mrededge7m_ndi <- with(zonal_stats_seg, (S2_summer_rededge6_median - S2_summer_rededge7_median)/(S2_summer_rededge6_median + S2_summer_rededge7_median))
# zonal_stats_seg$s2_summer_redrededge5_di <- with(zonal_stats_seg, (S2_summer_red - S2_summer_rededge5))
# zonal_stats_seg$s2_summer_redswir2_di <- with(zonal_stats_seg, (S2_summer_red - S2_summer_swir2))
# zonal_stats_seg$s2_summer_bluegreen_di <- with(zonal_stats_seg, (S2_summer_blue - S2_summer_green))
# zonal_stats_seg$S2_summer_blue_S2_summer_red_ri <- with(zonal_stats_seg, (S2_summer_blue / S2_summer_red))
# zonal_stats_seg$s2_summer_blueswir1_di <- with(zonal_stats_seg, (S2_summer_blue - S2_summer_swir1))
# zonal_stats_seg$s2_summer_redswir1_di <- with(zonal_stats_seg, (S2_summer_red - S2_summer_swir1))
# zonal_stats_seg$s2_summer_greenrededge5_di <- with(zonal_stats_seg, (S2_summer_green - S2_summer_rededge5))
# zonal_stats_seg$s2_summer_redmSWIR1m_di <- with(zonal_stats_seg, (S2_summer_red_median - S2_summer_swir1_median))
# zonal_stats_seg$S2_summer_green_S2_summer_swir1_ri <- with(zonal_stats_seg, (S2_summer_green / S2_summer_swir1))
# 
# # Seasonal difference indices
# zonal_stats_seg$S2_summer_swir2m_S2_winter_rededge5m_di <- with(zonal_stats_seg, (S2_summer_swir2_median - S2_winter_rededge5_median))
# zonal_stats_seg$S2_summer_swir1_S2_winter_swir1_di <- with(zonal_stats_seg, (S2_summer_swir1 - S2_winter_swir1))
# zonal_stats_seg$S2_summer_rededge6_S2_winter_swir1_di <- with(zonal_stats_seg, (S2_summer_rededge6 - S2_winter_swir1))
# zonal_stats_seg$S2_summer_swir1m_S2_winter_swir1m_di <- with(zonal_stats_seg, (S2_summer_swir1_median - S2_winter_swir1_median))
# zonal_stats_seg$S2_summer_redm_S2_winter_rededge5m_di <- with(zonal_stats_seg, (S2_summer_red_median - S2_winter_rededge5_median))
# zonal_stats_seg$S2_summer_swir2sd_S2_winter_redsd_di <- with(zonal_stats_seg, (S2_summer_swir2_sd - S2_winter_red_sd))
# zonal_stats_seg$LSU_summer_PV_LSU_winter_S_ndi <- with(zonal_stats_seg, (LSU_summer_PV - LSU_winter_S)/(LSU_summer_PV + LSU_winter_S))
# zonal_stats_seg$S2_summer_nir_S2_winter_rededge6_di <- with(zonal_stats_seg, (S2_summer_nir - S2_winter_rededge6))
# zonal_stats_seg$LSU_summer_NPV_LSU_winter_NPV_di <- with(zonal_stats_seg, (LSU_summer_NPV - LSU_winter_NPV))
# zonal_stats_seg$S2_summer_nirsd_S2_winter_rededge7sd_di <- with(zonal_stats_seg, (S2_summer_nir_sd - S2_winter_rededge7_sd))
# zonal_stats_seg$S2_summer_redsd_S2_winter_bluesd_di <- with(zonal_stats_seg, (S2_summer_red_sd - S2_winter_blue_sd))


# Ensure that catagorical data doesn't having any missing or inf values
zonal_stats_seg[is.na(zonal_stats_seg)] <- 0
zonal_stats_seg[sapply(zonal_stats_seg, is.infinite)] <- 0

write.table(zonal_stats_seg, "zonal_stats/zonal_stats_seg_20170801_MASKED_NORMAL.txt", sep="\t")

### Training data (Zonal Stats) #######################################################################
#
#

nmax <- 30 #### Number of training points per class

if (!exists("zonal_stats_seg"))
{
  zonal_stats_seg <- read.table("zonal_stats/zonal_stats_seg_20170801_MASKED_NORMAL.txt", sep="\t", header=T, as.is=T)
}
#segmentation.raster <-raster("Segmentation/Segmentation_LARGE_Rasterised.tif")
segmentation.raster <-raster("Segmentation/Segmentation_Rasterised.tif")
#segmentation.raster <-raster("Segmentation/Test_Segmentation_Rasterised.tif")

training.data.shp <- training.data.habitat.shp[c("PRIHAB", "DOMINANCE")]

# Identify the segmented polygons the training points fall within and extract the zonal statistics from these
training.data.ids <- extract(segmentation.raster, training.data.shp)

training.data <- merge(data.frame(ID=1:nrow(training.data.shp), training.data.shp, seg.id=training.data.ids), zonal_stats_seg, by.x="seg.id", by.y="ID")

training.data.all <- training.data[c(2:4,8:ncol(training.data))]

# Remove rows with mising values
training.data.all <- training.data.all[complete.cases(training.data.all),]
training.data.all$DOMINANCE <- as.numeric(training.data.all$DOMINANCE)

# Select the training data for points with accurate spatial mapping and for mappable habitat classes
training.data.all <- subset(training.data.all, DOMINANCE<=2)

# Split into stratified training and test datasets
training.data <- NULL
training.data.test <- NULL

# Loop through all the classes
set.seed(1) # Set a random seed to ensure consistent training and test datasets
for(c in unique(training.data.all$PRIHAB))
{
  # Select the subset of rows for the current class
  training.data.sub <- subset(training.data.all, PRIHAB==c)
  
  # Select a sample prioritising the training points from the highest tier
  n <- nrow(training.data.sub)
  prb <- ifelse(training.data.sub$DOMINANCE == 1,0.85, ifelse(training.data.sub$DOMINANCE == 2, 0.75, 0.25))
  training.data.sub <- training.data.sub[sample(n, min(n,nmax), prob=prb, replace=F),]
  
  # Split the data using a random sample
  subset <- random.subset(training.data.sub, 0.8)
  training.data <- rbind(training.data, training.data.sub[subset,])
  training.data.test <- rbind(training.data.test, training.data.sub[-subset,])
  rownames(training.data.test) <- NULL
}

# Remove duplicates from test dataset
training.data.test <- training.data.test[!duplicated(training.data.test),]

#duplicated(rbind(training.data.test, training.data))

## Write training_data to text file (training_data.txt)
write.table(training.data, "training_data/training_data_20170801_MASKED_NORMAL.txt", sep="\t")
write.table(training.data.test, "training_data/training_data_test_20170801_MASKED_NORMAL.txt", sep="\t")

### Classify training points using random forest ######################################################
#
#

if (!exists("zonal_stats_seg"))
{
  zonal_stats_seg <- read.table("zonal_stats/zonal_stats_seg_20170801_MASKED_NORMAL.txt", sep="\t", header=T, as.is=T)
}

# Read in training and test datasets
training.data <- read.table("training_data/training_data_20170801_MASKED_NORMAL.txt", sep="\t", header=T)
training.data.test <- read.table("training_data/training_data_test_20170801_MASKED_NORMAL.txt", sep="\t", header=T)

training.data$PRIHAB <- as.factor(as.character(training.data$PRIHAB))

# Predict detailed habitats using random forest
#Run for only the top 42 most important variables
M.rf.detailed.all <- randomForest(PRIHAB ~ ., data=training.data[c(2,4:ncol(training.data))], na.action=na.omit)
i <- colnames(training.data) %in% c(rownames(M.rf.detailed.all$importance)[order(M.rf.detailed.all$importance, decreasing=T)][1:90],"PRIHAB")
M.rf.detailed <- randomForest(PRIHAB ~ ., data=training.data[i], na.action=na.omit)

# Calculate confusion matrix
p <- predict(M.rf.detailed, training.data.test, type="response")
confusion.matrix(training.data.test$PRIHAB, p)


# Predict classes for all polygons
results.detailed.probs <- predict(M.rf.detailed, zonal_stats_seg,
                                  type="vote", norm.votes=TRUE,
                                  progress="text")

responseNFromProbs <- function(df, n=1) {
  columns <- colnames(df)
  response <- apply(df, MARGIN=1, FUN=function(x) {columns[order(x, decreasing=TRUE)[n]]})
  return (response)
}

probNFromProbs <- function(df, n=1) {
  response <- apply(df, MARGIN=1, FUN=function(x) {sort(x, decreasing=TRUE)[n]})
  return (response)
}

results.detailed.response1 <- responseNFromProbs(results.detailed.probs, n=1)
results.detailed.prob1 <- probNFromProbs(results.detailed.probs, n=1)*100
results.detailed.response2 <- responseNFromProbs(results.detailed.probs, n=2)
results.detailed.prob2 <- probNFromProbs(results.detailed.probs, n=2)*100


# Combine results with segmentation polygons and save to new shapefile
results.rf <- data.frame(ID=zonal_stats_seg$ID,
                         A_pred=results.detailed.response1,
                         A_prob=results.detailed.prob1,
                         B_pred=results.detailed.response2,
                         B_prob=results.detailed.prob2)
segmentation.p <- merge(segmentation.shp, results.rf, by="ID")
writeOGR(segmentation.p,
         "Outputs/CommonsMapping_20170801_MASKED_NORMAL.shp",
         "CommonsMapping_20170801_MASKED_NORMAL",
         driver="ESRI Shapefile",
         overwrite=T)
rm(segmentation.p)


### Graphs ###########################################################################################
#
#

#Create a variable called confusion matrix - broad classes
cm1 <- broadclass.confusion.matrix(training.data.test$PRIHAB, training.data.test$PRIHAB, p)

C_graph1 <- barplot.confusion.matrix(cm1)

#Adjust margins
C_graph1 <- C_graph1 + theme(plot.margin = unit(c(1,1,1,1), "cm"))

C_graph1

#Create a variable called confusion matrix 2 - detailed classes
cm2 <- broadclass.confusion.matrix(training.data.test$PRIHAB, training.data.test$PRIHAB, p)

C_graph2 <- barplot.confusion.matrix(cm2)

#Adjust margins
C_graph2 <- C_graph2 + theme(plot.margin = unit(c(1,1,1,3), "cm"))

C_graph2

## Percentage plots ###############


# Create a simple percentage plot - broad classes
C_graph3 <- barplot.percent(cm1)

#Adjust margins
C_graph3 <- C_graph3 + theme(plot.margin = unit(c(1,1,1,1), "cm"))

C_graph3


# Create a simple percentage plot - detailed classes
C_graph4 <- barplot.percent(cm2)

#Adjust margins
C_graph4 <- C_graph4 + theme(plot.margin = unit(c(1,1,1,1), "cm"))

C_graph4