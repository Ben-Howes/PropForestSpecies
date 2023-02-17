###### Calculate proportion of forest cover
###### and land area for each cell

library(tidyverse)
library(terra)

gpath = "/home/ben/Documents/PhD/matrix_Response/"
setwd(gpath)

## Load forest map
hansen_map = rast(paste0(gpath, "Data/Ranges/forest_map.tif"))
forest = hansen_map["treecover"] 
land = hansen_map["datamask"]

## Choose only forest cells > 70% tree cover
forest = classify(forest, cbind(0, 69, 0), include.lowest = T) %>% classify(., cbind(70, 100, 1), include.lowest = T)

## Choose only land cells which == 1 which means land
land = classify(land, cbind(2, 0))

## Convert both to behrmann projection and to template raster resolution
behr = "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
template_raster = rast() %>% project(behr)

forest = project(forest, behr, method = "near")
forest = terra::resample(forest, template_raster, method = "average")

land = project(land, behr, method = "near")
land = terra::resample(land, template_raster, method = "average")

prop_forest_area = forest/land

## Around 15 pixels had values over 1, nearly all of those are small rounding erros leading to values of 1.02 etc
## But one has a value of 2 because the land area is so small, something around 0.0001
## We will filter by minimum % of pixel made up by land in analysis, so make them all == 1 if they're over 1

prop_forest_area = classify(prop_forest_area, cbind(1, Inf, 1))
forest_rasts = c(forest, land, prop_forest_area)
names(forest_rasts) = c("forest", "land", "prop_forest_area")

writeRaster(forest_rasts, paste0(gpath, "Data/Ranges/forest_rasts.tif"))
