#######################################################################
## Calculate the % of disturbance occuring at each latitudinal band
## from Betts et al 2019
#######################################################################

library(tidyverse)
library(sf)
library(raster)
library(maps)
library(terra)

gpath = "/home/ben/Documents/PhD/matrix_Response/"
setwd(gpath)

## Load in the disturbance shapefiles
fires = st_read(paste0(gpath, "Data/Disturbances/Fires.shp")) %>% st_make_valid()
glaciers = st_read(paste0(gpath, "Data/Disturbances/Glaciers.shp"))
storms = st_read(paste0(gpath, "Data/Disturbances/Storms.shp"))

## Transform all into the behrmann projection and then rasterize to template
behr = "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
template_raster = rast() %>% project(behr)

fires = st_transform(fires, behr) %>% terra::rasterize(template_raster, touches = TRUE, background = 0)
glaciers = st_transform(glaciers, behr) %>% terra::rasterize(template_raster, touches = TRUE, background = 0)
storms = st_transform(storms, behr) %>% terra::rasterize(template_raster, touches = TRUE, background = 0)

dist_list = list(fires, glaciers, storms)
names = c("fires", "glaciers", "storms")

dist_rasts = rast(dist_list)
names(dist_rasts) = names

writeRaster(dist_rasts, paste0(gpath, "Data/Ranges/dist_rasts.tif"), overwrite = TRUE)
