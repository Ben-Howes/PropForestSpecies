##########################################################
## Plot proportion of forest vs matrix species
## from latitudes and habitat usage taken from 
## IUCN for mammals, amphibs, and reptiles
## and from AVONET for birds
##########################################################

library(tidyverse)
library(janitor)
library(sf)
library(lemon)
library(mgcv)
library(lme4)
library(ggeffects)
library(MuMIn)
library(terra)
library(maptools)
library(colorspace)
library(tidyterra)

gpath = "/home/ben/Documents/PhD/PropForestSpecies/"
setwd(gpath)

## Make behrmann template raster for use in transforming in future
behr = "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
template_raster = rast() %>% project(behr)

## Load basic world outline from GSHHS to use in cropping/masking
coast = st_read("../Raw_Data/gshhg-shp-2.3.7/GSHHS_shp/c/GSHHS_c_L1.shp")
coast = st_transform(coast, behr)
coast_outline = st_cast(coast, "MULTILINESTRING")

## Load all latitude and habitat data
forest_prop_maps = rast(paste0(gpath, "Data/Ranges/richness_rasts_resident_breeding.tif"))

ggplot() + geom_spatraster(data = forest_prop_maps[[5:8]]) + geom_sf(data = coast_outline) + theme_classic() +
  scale_fill_viridis_c(na.value = "white") + labs(fill = "Proportion Forest Species") +
  facet_wrap(~lyr, ncol = 2) 

################################################
## Are these trends explained by forest area?
################################################

## Load forest area data calculated from Hansen 2013
## Forest categorised at 30, 50, and 70% tree cover cut-offs
## For now just use 70% cut off values and later use others for sensitivity analysis

forest_rasts = rast(paste0(gpath, "Data/Ranges/forest_rasts.tif"))

forest_rasts = mask(forest_rasts, coast)

ggplot() + geom_spatraster(data = forest_rasts[[3]]) + geom_sf(data = coast_outline) + theme_classic() +
  scale_fill_viridis_c(na.value = "white") + labs(fill = "Proportion Forest Area")

## Add this raster to the stack of proportion forest species rasters
forest_prop_maps = c(forest_prop_maps, forest_rasts)

######################################################
## Could differences be due to historic forest loss?
######################################################
## Load in data from Betts et al tranformed into behr projection and % calculated at 1 cell
## area 96x96km
loss_rast = rast(paste0(gpath, "Data/Ranges/loss_rast.tif"))
loss_rast = mask(loss_rast, coast)

ggplot() + geom_spatraster(data = loss_rast) + geom_sf(data = coast_outline) + scale_fill_viridis_c(na.value = "white") +
  theme_classic() + labs(fill = "Proportion of Land Area Deforested")

## Add this raster to the stack of proportion forest species rasters
forest_prop_maps = c(forest_prop_maps, loss_rast)

######################################################
## Could differences be due to historic disturbances?
######################################################
## Load in data from Betts et al transformed into behrmann and categorised as 1 (disturbed) or 0 (undisturbed)
dist_rast = rast(paste0(gpath, "Data/Ranges/dist_rast.tif"))
dist_rast = mask(dist_rast, coast)

ggplot() + geom_spatraster(data = dist_rast) + geom_sf(data = coast_outline) + scale_fill_viridis_c(na.value = "white") +
  theme_classic() + labs(fill = "Disturbances")

## Add this raster to the stack of proportion forest species rasters
forest_prop_maps = c(forest_prop_maps, dist_rast)

######################################################
## How about how long an area has had that same type 
## of forest?
######################################################

eocene = st_read(paste0(gpath, "Data/Ranges/geological_forest_type/eocene.shp")) %>% st_transform("EPSG:4326")
oligocene = st_read(paste0(gpath, "Data/Ranges/geological_forest_type/oligocene.shp")) %>% st_transform("EPSG:4326")
miocene = st_read(paste0(gpath, "Data/Ranges/geological_forest_type/miocene.shp")) %>% st_transform("EPSG:4326")
pliocene = st_read(paste0(gpath, "Data/Ranges/geological_forest_type/pliocene.shp")) %>% st_transform("EPSG:4326")
holocene = st_read(paste0(gpath, "Data/Ranges/geological_forest_type/holocene.shp"))

eocene = project(as(eocene, 'SpatVector'), behr)
oligocene = project(as(oligocene, 'SpatVector'), behr)
miocene = project(as(miocene, 'SpatVector'), behr)
pliocene = project(as(pliocene, 'SpatVector'), behr)
holocene = project(as(holocene, 'SpatVector'), behr)

eocene = terra::rasterize(eocene, template_raster, field = "type")
oligocene = terra::rasterize(oligocene, template_raster, field = "type")
miocene = terra::rasterize(miocene, template_raster, field = "type")
pliocene = terra::rasterize(pliocene, template_raster, field = "type")
holocene = terra::rasterize(holocene, template_raster, field = "type")

forest_type = c(eocene, oligocene, miocene, pliocene, holocene)
forest_type = as.data.frame(forest_type, xy = T) %>% rename("eocene" = 3, "oligocene" = 4, "miocene" = 5, "pliocene" = 6, "holocene" = 7)
forest_type = forest_type %>% mutate(geological_forest_time = case_when(eocene == oligocene & oligocene == miocene & miocene == pliocene & pliocene == holocene ~ 5, 
oligocene == miocene & miocene == pliocene & pliocene == holocene ~ 4 , miocene == pliocene & pliocene == holocene ~ 3, 
pliocene == holocene ~ 2)) %>% mutate(geological_forest_time = ifelse(is.na(geological_forest_time), 1, geological_forest_time))

forest_type$geological_forest_stability = 0

for (i in 1:nrow(forest_type)) {
  type = forest_type[i,7]
  count = 1
  if(!is.na(type)) {
    for (j in 3:6) {
      if(!is.na(forest_type[i,j])) {
        if(type == forest_type[i,j]) {count = count + 1}
      }
    }
  }
  forest_type[i,"geological_forest_stability"] = count
}

## Select just x y and geological forest, where the larger the value of geological forest the longer
## the forest in that area has continued to be the same forest type
forest_type_rast = rast(dplyr::select(forest_type, x,y,geological_forest_time, geological_forest_stability), crs = behr)
forest_type_rast = mask(forest_type_rast, coast) 

ggplot() + geom_spatraster(data = forest_type_rast) +
  geom_sf(data = coast_outline) +
  theme_classic() + scale_fill_viridis_c(na.value = "white") + 
  labs(fill = "Geological Forest Time")

## Add this raster to the stack of proportion forest species rasters
forest_prop_maps = c(forest_prop_maps, forest_type_rast)

######################################################
## Add data on the alpha diversity of plants
## which is a good proxy for structural diversity
######################################################

alpha = rast(paste0(gpath, "Data/Ranges/plantAlpha.tif"))
forest_prop_maps = c(forest_prop_maps, alpha)

######################################################
## Add data on the mean altitude of the pixels
######################################################

altitude = rast(paste0(gpath, "Data/Ranges/altitude.tif"))
forest_prop_maps = c(forest_prop_maps, altitude)

######################################################
## Mask by forest ecoregions as we only expect forest
## species in forest ecoregions
## as well as geological forest type and
## structural complexity only make sense in
## forest ecoregions
#######################################################

## Load in ecoregions as we just want to show data for forested biomes
ecoregions = st_read(paste0(gpath, "../Forest_thresholds/Data/Ecoregions2017/Ecoregions2017.shp"))

## Convert to behrmann projection
ecoregions = st_transform(ecoregions, behr)

## Just keep shapefiles with the word "forest" in them
ecoregions = ecoregions %>%
    mutate(forest = ifelse(grepl("forest", BIOME_NAME, ignore.case = TRUE), 1, 0)) %>%
    filter(forest == 1)

## Create a raster for the ecoregions so we can give each pixel
## a factor value of which ecoregion they were in - could be useful for spatial autocorrelation
behr = "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
template_raster = rast() %>% project(behr)

ecoNameRast = rasterize(ecoregions, template_raster, field = "ECO_NAME", touches = TRUE)
realmRast = rasterize(ecoregions, template_raster, field = "REALM", touches = TRUE)
biomeNameRast = rasterize(ecoregions, template_raster, field = "BIOME_NAME", touches = TRUE)

forest_prop_maps = c(forest_prop_maps, ecoNameRast, realmRast, biomeNameRast)

## Then we need to crop our predicted threshold type raster
## by these ecoregions, so we only keep data from the forest ecoregions
forest_prop_maps = mask(forest_prop_maps, ecoregions)

######################################################
## Analysis
######################################################

## Turn these rasters into a data frame for use in models
total_df = as.data.frame(forest_prop_maps, xy = T) %>% dplyr::select(-prop_amphibians, -prop_mammals, -prop_reptiles, -prop_birds) %>% 
  pivot_longer(cols = c(total_amphibians, total_mammals, total_reptiles, total_birds), names_to = "taxa", values_to = "n_spec") %>%
  mutate(taxa = sub("^[^_]*_", "", taxa))

prop_forest_df = as.data.frame(forest_prop_maps, xy = T) %>% dplyr::select(-total_amphibians, -total_mammals, -total_reptiles, -total_birds) %>%
  pivot_longer(cols = c(prop_amphibians, prop_mammals, prop_reptiles, prop_birds), names_to = "taxa", values_to = "prop_forest") %>%
  mutate(taxa = sub("^[^_]*_", "", taxa))

prop_forest_df = prop_forest_df %>% left_join(total_df) %>% relocate(x,y,prop_forest,n_spec,prop_forest_area,taxa)

## Remove rows that have no species and no land area
## Also use a cut off for what cells are included based on the total land area
## In this case a minimum of 10% land area
prop_forest_df = prop_forest_df %>% filter(!is.na(n_spec) & !is.na(land) & land > 0.1) %>%
  mutate(prop_forest = ifelse(is.na(prop_forest), 0, prop_forest), 
  prop_forest_area = ifelse(is.na(prop_forest_area), 0, prop_forest_area), 
  prop_land_area_deforested = ifelse(is.na(prop_land_area_deforested), 0, ifelse(prop_land_area_deforested > 1, 1, prop_land_area_deforested)),
  dist_equator_1000km = abs(y/1000000),
  disturbances = disturbances, alpha_plant_diversity = ifelse(is.na(alpha_plant_diversity), 0, alpha_plant_diversity), 
  geological_forest_time = case_when(geological_forest_time == 5 ~ 55, geological_forest_time == 4 ~ 33, geological_forest_time == 3 ~ 23,
  geological_forest_time == 2 ~ 5, geological_forest_time == 1 ~ 0), geological_forest_time = ifelse(is.na(geological_forest_time), 0, geological_forest_time))

## Remove areas with small proportion of forest species for plotting
prop_forest_df_plot = prop_forest_df %>% filter(n_spec >= 10 & land > 0.1)

write_csv(prop_forest_df, paste0(gpath, "Data/proportion_forest_species_analysis_data.csv"))
write_csv(prop_forest_df_plot, paste0(gpath, "Data/proportion_forest_species_plot_data.csv"))
