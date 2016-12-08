#Image classification using Random Forest
#Written by Christine Swanson
##November 22, 2016

##Load libraries
library(rgdal)
library(raster)
library(RStoolbox)
library(caret)
library(randomForest)
library(e1071)
library(snow)

##Create a files list for each image
rasters_2013217 <- list.files(path = "./Data", 
                              pattern = "LC82220672013217LGN00_B.*.TIF",
                              full.names = TRUE)

##Stack raster images
l8_2013217 <- stack(rasters_2013217)

##Subset raster list for tasseled cap analysis
l8_2013217_tca_stack <- stack(rasters_2013217[c(2, 3, 4, 5, 6, 7)])

##Load specific NIR and red bands for NDVI calculation
l8_2013217_nir <- stack("./Data/LC82220672013217LGN00_B5.TIF")
l8_2013217_red <- stack("./Data/LC82220672013217LGN00_B4.TIF")

##Rename bands in raster stacks
names(l8_2013217) <- c(paste0("B", 1:7, coll = ""), "B9")
names(l8_2013217_tca_stack) <- c(paste0("B", 2:7, coll = ""))


##Calculate NDVI from red band (band 1) and NIR band (band 2)
l8_2013217_ndvi <- overlay(l8_2013217_red, l8_2013217_nir, fun = function(Red, NIR) {
  (NIR-Red) / (NIR+Red)
})

##Rename NDVI layer
names(l8_2013217_ndvi) <- c("NDVI")

##Run a tasseled cap analysis on the raster stack
l8_2013217_tca <- tasseledCap(l8_2013217_tca_stack, "Landsat8OLI")

##Rename TCA layers
names(l8_2013217_tca) <- c("Brightness", "Greenness", "Wetness")

##Load the elevation dataset
elev <- raster("./Data/elev_mos_clp1.tif")

##Re-project the elev raster to match the Landsat data (for some reason crs = l8_2013217 wasn't working)
new_projection <- "+proj=utm +zone=22 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
elev_utm22 <- projectRaster(elev, crs=new_projection)

##Resample the elevation raster to match the Landsat rasters
elev_utm22_resamp <- resample(elev_utm22, l8_2013217, method = "bilinear")

##Rename the elevation raster
names(elev_utm22_resamp) <- c("elev")

##Calculate the slope
l8_2013217_slope <- terrain(elev_utm22_resamp, opt=c('slope'), unit = 'degrees')

##Stack all of the rasters for the analysis
l8_2013217_rf_data <- stack(l8_2013217, l8_2013217_ndvi, elev_utm22_resamp, l8_2013217_slope, l8_2013217_tca)

##Load the training data
trainData <- readOGR(dsn = "./Data/Training_Points", layer = "TO_training_points", pointDropZ = TRUE)
trainData_utm22 <- spTransform(trainData, crs(l8_2013217))
responseCol <- "Name"

extract_pixels <- function(img, training_data){
  ##Extract training pixel values
  dfAll = data.frame(matrix(vector(), nrow = 0, ncol = length(names(img)) + 1))
  for (i in 1:length(unique(training_data[[responseCol]]))){
    category <- unique(training_data[[responseCol]])[i]
    categorymap <- training_data[training_data[[responseCol]] == category,]
    dataSet <- extract(img, categorymap)
  
  if(is(training_data, "SpatialPointsDataFrame")){
    dataSet <- cbind(dataSet, class = as.numeric(category))
    dfAll <- rbind(dfAll, dataSet)
    }
  if(is(training_data, "SpatialPolygonsDataFrame")){
    dataSet <- lapply(dataSet, function(x){cbind(x, class = as.numeric(rep(category, nrow(x))))})
    df <- do.call("rbind", dataSet)
    dfAll <- rbind(dfAll, df)
    }
  }
  return(dfAll)
}

##Run the extract_pixels function for the data
l8_2013217_pixels <- extract_pixels(l8_2013217_rf_data, trainData_utm22)

l8_2013217_rf <- randomForest(as.factor(class) ~ B1 + B2 + B3 + B4 + B5 + B6 + B7 + elev + slope + NDVI + Brightness + Wetness + Greenness,
                              data = l8_2013217_pixels, mtry = 2, importance = TRUE, confusion = TRUE)

rf_predictions <- function(img, model_name){
  ##Create a raster with predictions from the fitted model project
  beginCluster()
  preds_rf <- clusterR(img, raster::predict, args = list(model = model_name))
  endCluster()
  return(preds_rf)
}

##Run random forest predictions for the RF model
l8_2013217_preds <- rf_predictions(l8_2013217_rf_data, l8_2013217_rf)

##Plot the variable importance
varImpPlot(l8_2013217_rf, main = "Variable Importance")

##Write out the rasters
writeRaster(l8_2013217_preds, filename = "l8_2013217_model_predictions", format = "GTiff", overwrite = TRUE)
writeRaster(l8_2013217_tca, filename = "l8_2013217_tca", format = "GTiff", overwrite = TRUE)
writeRaster(l8_2013217_slope, filename = "l8_2013217_slope", format = "GTiff", overwrite = TRUE)
writeRaster(elev_utm22_resamp, filename = "l8_2013217_elev", format = "GTiff", overwrite = TRUE)
writeRaster(l8_2013217_ndvi, filename = "l8_2013217_ndvi", format = "GTiff", overwrite = TRUE)
