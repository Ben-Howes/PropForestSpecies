############################################################
## Calculate forest loss at 1 degree bands
## Based on calculations in Betts et al 2017
############################################################

library(terra)
library(tidyverse)

gpath = "/home/ben/Documents/PhD/PropForestSpecies/"
setwd(gpath)

## Load in forest loss from Betts et al 2017
loss = rast(paste0(gpath, "Data/historical_forest_loss.tif"))

## Convert to behrmann equal area projection
behr = "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
template_raster = rast() %>% project(behr)
loss_rast = project(loss, template_raster, method = "bilinear")

## Load land raster
land = rast(paste0(gpath, "Data/Ranges/forest_rasts.tif"))[[2]]

## Divide deforested area by land area to get a % of land area that has been deforested
## Note that some of these pixels will be above 1 due to issues with areas with particularly small areas of land
## and due to the different initial resolutions of these raster compared to the hansen raster
## E.g some areas have Inf forest loss as the land area is 0 but the forest loss is 1x10-4
## Most of these shouldn't be an issue as we remove pixels with under 0.1 prop land area
loss_rast = loss_rast/land
names(loss_rast) = "prop_land_area_deforested"

writeRaster(loss_rast, paste0(gpath, "Data/Ranges/loss_rast.tif"), overwrite = TRUE)
