#######################################################################
## Transform data on alpha diversity of plants from Sabatini paper
## to behr for use as explanatory variable
#######################################################################

library(tidyverse)
library(sf)
library(raster)
library(maps)
library(terra)

gpath = "/home/ben/Documents/PhD/PropForestSpecies/"
setwd(gpath)

alpha = rast(paste0(gpath, "Data/alpha_plant_diversity.tif"))

## Transform all into the behrmann projection and then rasterize to template
behr = "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
template_raster = rast() %>% project(behr)

alphaBehr = project(alpha, template_raster)

writeRaster(alphaBehr, paste0(gpath, "Data/Ranges/plantAlpha.tif"), overwrite = TRUE)
