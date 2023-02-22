######################################################################################
## Plot proportion forest species of each taxa
## Small land area cells (<100km) have been removed to aid visualisation
## But no minimum number of species used
######################################################################################

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

## Read in data
prop_forest_df = read_csv(paste0(gpath, "Data/proportion_forest_species_analysis_data.csv"))

## Make behrmann template raster for use in transforming in future
behr = "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
template_raster = rast() %>% project(behr)

## Load basic world outline from GSHHS to use in cropping/masking
coast = st_read("../Raw_Data/gshhg-shp-2.3.7/GSHHS_shp/c/GSHHS_c_L1.shp")
coast = st_transform(coast, behr)
coast_outline = st_cast(coast, "MULTILINESTRING")

figure1_data = prop_forest_df %>% filter(!is.na(n_spec) & n_spec >= 10 & !is.na(land) & land > 0.1) %>%
  mutate(prop_forest = ifelse(is.na(prop_forest), 0, prop_forest))

figure1_data = dplyr::select(figure1_data, x, y, prop_forest, taxa) %>% group_split(taxa) %>% 
  map(., ~rast(dplyr::select(., x,y,prop_forest), crs = behr) %>% 
  project(template_raster)) %>% rast()

names(figure1_data) = c("amphibians", "birds", "mammals", "reptiles")

ggplot() + geom_spatraster(data = figure1_data) + facet_wrap(.~lyr, ncol = 2) + 
  geom_sf(data = coast_outline) +
  theme_classic() + scale_fill_viridis_c(na.value = "white") + 
  labs(fill = "Proportion Forest Species")
