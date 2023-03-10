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
library(ggimage)

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

coast_rast = coast %>% terra::rasterize(template_raster, touchess = TRUE, background = NA)

figure1_data = prop_forest_df

figure1_map_data = dplyr::select(figure1_data, x, y, prop_forest, taxa) %>% group_split(taxa) %>% 
  map(., ~rast(dplyr::select(., x,y,prop_forest), crs = behr) %>% 
  project(template_raster)) %>% rast()

## Set names of each raster layer
names(figure1_map_data) = c("Amphibians", "Birds", "Mammals", "Reptiles")

## Change facet labels
facet_labels =  c(
  "Amphibians"="A)",
  "Birds"="B)",
  "Mammals"="C)",
  "Reptiles"="D)"
)

## Create plot
maps = ggplot() + 
        geom_sf(data = coast, fill = "grey70") + 
        geom_spatraster(data = figure1_map_data) + 
        geom_sf(data = coast_outline, linewidth = 0.75) +
        facet_wrap(.~lyr, ncol = 2, labeller = as_labeller(facet_labels)) + 
        theme_classic() + 
        scale_fill_viridis_c(na.value = "transparent", breaks = c(0, 0.25, 0.5, 0.75),
        limits = c(0,0.75)) + 
        labs(fill = "Proportion of Forest Species") +
        theme(text = element_text(size = 30),
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
          # strip.text.x = element_text(hjust = 0, margin=margin(l=0))) +
        theme(legend.position="bottom") +
        guides(fill = guide_colourbar(barwidth=35, 
                      ticks = TRUE, frame.colour = "black", 
                      frame.linewidth = 0.5, title.position = "top",
                      title.hjust = 0.5, title.vjust = -0.5)) 

maps

ggsave(paste0(gpath, "Paper/Figures/prop_map.png"), plot = maps, width = 20, height = 10)
