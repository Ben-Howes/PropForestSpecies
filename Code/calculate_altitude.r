#######################################################################
## Transform data on altitude of an area into behrmann
#######################################################################

library(tidyverse)
# library(raster)
library(maps)
library(terra)

gpath = "/home/ben/Documents/PhD/PropForestSpecies/"
setwd(gpath)

altitude = rast(paste0(gpath, "Data/ETOPO_2022_v1_60s_N90W180_bed.tif"))
names(altitude) = "altitude"

## Transform all into the behrmann projection and then rasterize to template
behr = "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
template_raster = rast() %>% project(behr)

altitude = project(altitude, template_raster)

writeRaster(altitude, paste0(gpath, "Data/Ranges/altitude.tif"), overwrite = TRUE)
