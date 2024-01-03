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
library(cowplot)

gpath = "/home/ben/Documents/PhD/PropForestSpecies/"
setwd(gpath)

## Read in data
prop_forest_df = read_csv(paste0(gpath, "Data/proportion_forest_species_plot_data.csv"))

## Make behrmann template raster for use in transforming in future
behr = "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
template_raster = rast() %>% project(behr)

## Load basic world outline from GSHHS to use in cropping/masking
coast = st_read("../Raw_Data/gshhg-shp-2.3.7/GSHHS_shp/c/GSHHS_c_L1.shp")
coast = st_transform(coast, behr)
coast_outline = st_cast(coast, "MULTILINESTRING")

coast_rast = coast %>% terra::rasterize(template_raster, touchess = TRUE, background = NA)

figure1_map_data = dplyr::select(prop_forest_df, x, y, prop_forest, taxa) %>% group_split(taxa) %>% 
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

################################################
## Plot with all four taxa maps together
################################################

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

# ggsave(paste0(gpath, "Paper/Figures/prop_map.png"), plot = maps, width = 20, height = 10)

################################################
## Plot individual maps and then grid them together
################################################

## Make amphibian map
amphib_map = ggplot() + 
            geom_sf(data = coast, fill = "grey70") + 
            geom_spatraster(data = figure1_map_data$Amphibians) + 
            geom_sf(data = coast_outline, linewidth = 0.75) +
            theme_void() + 
            scale_fill_viridis_c(na.value = "transparent", breaks = c(0, 0.25, 0.5, 0.75),
            limits = c(0,0.75)) + 
            labs(fill = "Proportion of Forest Species") +
            theme(legend.position = "none") +
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0)) +
            geom_hline(yintercept = 0, linetype = "dashed") +
            # geom_hline(yintercept = max(prop_forest_df[prop_forest_df$taxa == "amphibians",]$y)) +
            # geom_hline(yintercept = min(prop_forest_df[prop_forest_df$taxa == "amphibians",]$y)) +
            theme(plot.margin = unit(c(0.25,0,0.68,0), "cm"))

amphib_map

amphib_density = ggplot(filter(prop_forest_df, taxa == "amphibians"), aes(y, prop_forest)) + 
  stat_smooth(aes(col = after_stat(y)), se = F, linewidth = 5, n = 2500) +
  coord_flip(ylim = c(0, 0.4)) +
  theme_classic() +
  scale_colour_viridis_c(limits = c(0 ,0.75)) +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL) +
  theme(text = element_text(size = 20),
  axis.text= element_blank(),
  axis.ticks = element_blank(),
  axis.line = element_line(linewidth = 1),
  axis.line.x = element_blank()) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  # geom_vline(xintercept = max(prop_forest_df[prop_forest_df$taxa == "amphibians",]$y)) +
  # geom_vline(xintercept = min(prop_forest_df[prop_forest_df$taxa == "amphibians",]$y)) +
  theme(plot.margin = unit(c(0.5,0.5,0.86,0), "cm")) +
  xlim(-7367884, 7342230)

amphib_density

a = plot_grid(amphib_map, amphib_density, rel_widths = c(1, 0.2), labels = c("A1", "A2"), label_size = 35)

ggsave(paste0(gpath, "Paper/Figures/TaxaMaps/amphibian_map.png"), plot = a, width = 20, height = 8.3)

## Make bird map
bird_map = ggplot() + 
            geom_sf(data = coast, fill = "grey70") + 
            geom_spatraster(data = figure1_map_data$Birds) + 
            geom_sf(data = coast_outline, linewidth = 0.75) +
            theme_void() + 
            scale_fill_viridis_c(na.value = "transparent", breaks = c(0, 0.25, 0.5, 0.75),
            limits = c(0,0.75)) + 
            labs(fill = "Proportion of Forest Species") +
            theme(legend.position = "none") +
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0)) +
            geom_hline(yintercept = 0, linetype = "dashed") +
            # geom_hline(yintercept = max(prop_forest_df[prop_forest_df$taxa == "birds",]$y)) +
            # geom_hline(yintercept = min(prop_forest_df[prop_forest_df$taxa == "birds",]$y)) +
            theme(plot.margin = unit(c(0.25,0,0.68,0), "cm"))

bird_map

bird_density = ggplot(filter(prop_forest_df, taxa == "birds"), aes(y, prop_forest)) + 
  stat_smooth(aes(col = after_stat(y)), se = F, linewidth = 5, n = 2500) +
  coord_flip(ylim = c(0, 0.4)) +
  theme_classic() +
  scale_colour_viridis_c(limits = c(0, 0.75)) +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL) +
  theme(text = element_text(size = 20),
  axis.text= element_blank(),
  axis.ticks = element_blank(),
  axis.line = element_line(linewidth = 1),
  axis.line.x = element_blank()) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  # geom_vline(xintercept = max(prop_forest_df[prop_forest_df$taxa == "birds",]$y)) +
  # geom_vline(xintercept = min(prop_forest_df[prop_forest_df$taxa == "birds",]$y)) +
  theme(plot.margin = unit(c(0.5,0.5,0.86,0), "cm")) +
  xlim(-7367884, 7342230)

bird_density

b = plot_grid(bird_map, bird_density, rel_widths = c(1, 0.2), labels = c("B1", "B2"), label_size = 35)

ggsave(paste0(gpath, "Paper/Figures/TaxaMaps/bird_map.png"), plot = b, width = 20, height = 8.3)

## Make mammal map
mammal_map = ggplot() + 
            geom_sf(data = coast, fill = "grey70") + 
            geom_spatraster(data = figure1_map_data$Mammals) + 
            geom_sf(data = coast_outline, linewidth = 0.75) +
            theme_void() + 
            scale_fill_viridis_c(na.value = "transparent", breaks = c(0, 0.25, 0.5, 0.75),
            limits = c(0,0.75)) + 
            labs(fill = "Proportion of Forest Species") +
            theme(legend.position = "none") +
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0)) +
            geom_hline(yintercept = 0, linetype = "dashed") +
            # geom_hline(yintercept = max(prop_forest_df[prop_forest_df$taxa == "mammals",]$y)) +
            # geom_hline(yintercept = min(prop_forest_df[prop_forest_df$taxa == "mammals",]$y)) +
            theme(plot.margin = unit(c(0.25,0,0.68,0), "cm"))

mammal_map

mammal_density = ggplot(filter(prop_forest_df, taxa == "mammals"), aes(y, prop_forest)) + 
  stat_smooth(aes(col = after_stat(y)), se = F, linewidth = 5, n = 2500) +
  coord_flip(ylim = c(0, 0.4)) +
  theme_classic() +
  scale_colour_viridis_c(limits = c(0, 0.75)) +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL) +
  theme(text = element_text(size = 20),
  axis.text= element_blank(),
  axis.ticks = element_blank(),
  axis.line = element_line(linewidth = 1),
  axis.line.x = element_blank()) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  # geom_vline(xintercept = max(prop_forest_df[prop_forest_df$taxa == "mammals",]$y)) +
  # geom_vline(xintercept = min(prop_forest_df[prop_forest_df$taxa == "mammals",]$y)) +
  theme(plot.margin = unit(c(0.5,0.5,0.86,0), "cm")) +
  xlim(-7367884, 7342230)

mammal_density

c = plot_grid(mammal_map, mammal_density, rel_widths = c(1, 0.2), labels = c("C1", "C2"), label_size = 35)

ggsave(paste0(gpath, "Paper/Figures/TaxaMaps/mammal_map.png"), plot = c, width = 20, height = 8.3)

## Make reptile map
reptile_map = ggplot() + 
            geom_sf(data = coast, fill = "grey70") + 
            geom_spatraster(data = figure1_map_data$Reptiles) + 
            geom_sf(data = coast_outline, linewidth = 0.75) +
            theme_void() + 
            scale_fill_viridis_c(na.value = "transparent", breaks = c(0, 0.25, 0.5, 0.75),
            limits = c(0,0.75)) + 
            labs(fill = "Proportion of Forest Species") +
            theme(legend.position = "none") +
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0)) +
            geom_hline(yintercept = 0, linetype = "dashed") +
            # geom_hline(yintercept = max(prop_forest_df[prop_forest_df$taxa == "reptiles",]$y)) +
            # geom_hline(yintercept = min(prop_forest_df[prop_forest_df$taxa == "reptiles",]$y)) +
            theme(plot.margin = unit(c(0.25,0,0.68,0), "cm"))

reptile_map

reptile_density = ggplot(filter(prop_forest_df, taxa == "reptiles"), aes(y, prop_forest)) + 
  stat_smooth(aes(col = after_stat(y)), se = F, linewidth = 5, n = 2500) +
  coord_flip(ylim = c(0, 0.4)) +
  theme_classic() +
  scale_colour_viridis_c(limits = c(0, 0.75)) +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL) +
  theme(text = element_text(size = 20),
  axis.text= element_blank(),
  axis.ticks = element_blank(),
  axis.line = element_line(linewidth = 1),
  axis.line.x = element_blank()) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  # geom_vline(xintercept = max(prop_forest_df[prop_forest_df$taxa == "reptiles",]$y)) +
  # geom_vline(xintercept = min(prop_forest_df[prop_forest_df$taxa == "reptiles",]$y)) +
  theme(plot.margin = unit(c(0.5,0.5,0.86,0), "cm")) +
  xlim(-7367884, 7342230)

reptile_density

d = plot_grid(reptile_map, reptile_density, rel_widths = c(1, 0.2), labels = c("D1", "D2"), label_size = 35)

ggsave(paste0(gpath, "Paper/Figures/TaxaMaps/reptile_map.png"), plot = d, width = 20, height = 8.3)
