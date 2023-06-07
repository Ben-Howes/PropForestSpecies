#####################################
## Make figures showing
## prop forest cover, historic loss
## disturbance and forest status
#####################################

library(tidyverse)
library(janitor)
library(sf)
library(lemon)
library(terra)
library(maptools)
library(colorspace)
library(tidyterra)

gpath = "/home/ben/Documents/PhD/PropForestSpecies/"
setwd(gpath)

## Load raster files
forest = rast(paste0(gpath, "Data/Ranges/forest_rasts.tif"))
loss = rast(paste0(gpath, "Data/Ranges/loss_rast.tif"))
dist = rast(paste0(gpath, "Data/Ranges/dist_rasts.tif"))

## Combine different dist types
dist = terra::mosaic(dist[[1]], dist[[2]], dist[[3]], fun = "max")
names(dist) = "disturbances"

## Load basic world outline from GSHHS to use in cropping/masking
coast = st_read("../Raw_Data/gshhg-shp-2.3.7/GSHHS_shp/c/GSHHS_c_L1.shp")
behr = "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
coast = st_transform(coast, behr)
coast_outline = st_cast(coast, "MULTILINESTRING")

forest = mask(forest, coast, touches = F)
loss = mask(loss, coast, touches = F)
dist = mask(dist, coast, touches = F)

forest_map = ggplot() + 
        geom_spatraster(data = forest$prop_forest_area) + 
        geom_sf(data = coast_outline, linewidth = 0.75) +
        theme_classic() + 
        scale_fill_viridis_c(na.value = "transparent", breaks = c(0, 0.25, 0.5, 0.75, 1),
        limits = c(0,1)) + 
        labs(fill = "Proportion of Forested Land") +
        theme(text = element_text(size = 30)) +
        theme(legend.position="bottom") +
        guides(fill = guide_colourbar(barwidth=35, 
                      ticks = TRUE, frame.colour = "black", 
                      frame.linewidth = 0.5, title.position = "top",
                      title.hjust = 0.5, title.vjust = -0.5)) 

forest_map

ggsave(paste0(gpath, "Paper/Figures/forest_map.png"), plot = forest_map, width = 20, height = 10)

forest_map = ggplot() + 
        geom_spatraster(data = forest$prop_forest_area) + 
        geom_sf(data = coast_outline, linewidth = 0.75) +
        theme_classic() + 
        scale_fill_viridis_c(na.value = "transparent", breaks = c(0, 0.25, 0.5, 0.75, 1),
        limits = c(0,1)) + 
        labs(fill = "Proportion of Forested Land") +
        theme(text = element_text(size = 30)) +
        theme(legend.position="bottom") +
        guides(fill = guide_colourbar(barwidth=35, 
                      ticks = TRUE, frame.colour = "black", 
                      frame.linewidth = 0.5, title.position = "top",
                      title.hjust = 0.5, title.vjust = -0.5)) 

forest_map

ggsave(paste0(gpath, "Paper/Figures/forest_map.png"), plot = forest_map, width = 20, height = 10)

loss_map = ggplot() + 
        geom_spatraster(data = loss$prop_land_area_deforested) + 
        geom_sf(data = coast_outline, linewidth = 0.75) +
        theme_classic() + 
        scale_fill_viridis_c(na.value = "transparent", breaks = c(0, 0.25, 0.5, 0.75, 1),
        limits = c(0,1.05)) + 
        labs(fill = "Proportion of Land Area Deforested") +
        theme(text = element_text(size = 30)) +
        theme(legend.position="bottom") +
        guides(fill = guide_colourbar(barwidth=35, 
                      ticks = TRUE, frame.colour = "black", 
                      frame.linewidth = 0.5, title.position = "top",
                      title.hjust = 0.5, title.vjust = -0.5)) 

loss_map

ggsave(paste0(gpath, "Paper/Figures/loss_map.png"), plot = loss_map, width = 20, height = 10)

dist_map = ggplot() + 
        geom_spatraster(data = dist$disturbances) + 
        geom_sf(data = coast_outline, linewidth = 0.75) +
        theme_classic() + 
        scale_fill_viridis_c(na.value = "transparent", breaks = c(0, 0.25, 0.5, 0.75, 1),
        limits = c(0,1.05)) +
        theme(text = element_text(size = 30)) +
        theme(legend.position="none") +
        guides(fill = guide_colourbar(barwidth=35, 
                      ticks = TRUE, frame.colour = "black", 
                      frame.linewidth = 0.5, title.position = "top",
                      title.hjust = 0.5, title.vjust = -0.5)) 

dist_map

ggsave(paste0(gpath, "Paper/Figures/dist_map.png"), plot = dist_map, width = 20, height = 10)
