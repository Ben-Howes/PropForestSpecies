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

## Multiply land area by forest loss to get a % of land area that has been deforested
loss_rast = loss_rast*land
names(loss_rast) = "prop_land_area_deforested"

writeRaster(loss_rast, paste0(gpath, "Data/Ranges/loss_rast.tif"), overwrite = TRUE)
