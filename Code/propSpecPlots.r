##########################################
## Plot of propoprtion forest species
## against latitude and species richness
## for manuscript
##########################################

library(tidyverse)
library(janitor)
library(lemon) ## facet_rep
library(sf)
library(ggpubr) ## ggarrange
library(mgcv) ## bam
library(terra) ## rast

gpath = "/home/ben/Documents/PhD/PropForestSpecies/"
setwd(gpath)

prop_forest_df = read_csv(paste0(gpath, "Data/proportion_forest_species_analysis_data.csv")) %>%
    mutate(taxa = str_to_title(taxa))

## Plot of proportion of forest species against latitude
## Need to convert from behrmann to latitude
behr = "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"

prop_forest_df = st_as_sf(prop_forest_df, coords = c("x", "y")) %>% st_set_crs(behr) %>%
    st_transform(., 4326) %>% mutate(x = st_coordinates(.)[,1], y = st_coordinates(.)[,2]) %>%
    as.data.frame() %>%
    filter(n_spec > 10 & land > 0.1)

## Scale species richness by dividing n_spec by
## maximum n_spec per data
lat_forest_df = prop_forest_df %>%
    group_by(taxa) %>%
    mutate(prop_n_spec = n_spec/max(n_spec)) %>%
    pivot_longer(cols = c("prop_forest", "prop_n_spec"), names_to = "prop_type", values_to = "prop")

latDotPlot = ggplot(lat_forest_df, aes(abs(y), prop, col = prop_type)) +
    geom_point(size = 2, alpha = 0.2) +
    theme_classic() +
    labs(x = "Absolute Latitude", y = NULL, col = NULL) +
    theme(legend.position = c(0.3885, 0.92),
    legend.spacing.y = unit(0.25, 'cm'),
    legend.box.background = element_rect(colour = "black", 
    linewidth = 1.1), legend.margin = margin(t=0,r=20,b=10,l=15)) +
    theme(text = element_text(size = 30), 
    legend.text = element_text(size = 25)) + 
    guides(col = guide_legend(byrow = TRUE)) +
    facet_rep_wrap(. ~ taxa, scales = "free_x") +
    geom_smooth(aes(group = prop_type), linewidth = 4.5, col = "black", method = "gam", se = FALSE) +
    geom_smooth(linewidth = 2.5, method = "gam", se = FALSE) +
    scale_y_continuous(expand = c(0, 0), limits = c(-0.01, 1), breaks = c(0, 0.25, 0.5, 0.75)) +
    scale_x_continuous(expand = c(0, 0), limits = c(-0.01, 71)) +
    scale_color_manual(values = c("#3b528b", "#f8e621"), labels = c("Proportion\nForest Species", "Proportion Max\nSpecies Richness"))

latDotPlot

ggsave(paste0(gpath, "Paper/Figures/latitudeDotPlot.png"), plot = latDotPlot, width = 20, height = 12)

## Plot of proportion of forest species against species richness
propVsRich = ggplot(prop_forest_df, aes(n_spec, prop_forest)) +
    geom_point(size = 2, alpha = 0.2, col = "#3b528b") +
    theme_classic() +
    labs(x = "Total Species Richness", y = "Proportion Forest Species") +
    theme(text = element_text(size = 30)) +
    facet_rep_wrap(. ~ taxa, scales = "free_x") +
    geom_smooth(linewidth = 4.5, col = "black", method = "gam", se = F) +
    geom_smooth(linewidth = 2.5, col = "#3b528b", method = "gam") +
    scale_y_continuous(expand = c(0, 0), limits = c(-0.01, 1), breaks = c(0, 0.25, 0.5, 0.75)) +
    scale_x_continuous(expand = c(0.01, 0)) 

propVsRich

ggsave(paste0(gpath, "Paper/Figures/propVsRich.png"), plot = propVsRich, width = 20, height = 12)

ggarrange(latDotPlot + theme(legend.position = c(0.315, 0.935), legend.text = element_text(size = 15)), 
    propVsRich + labs(y = NULL), ncol = 2, nrow = 1, labels = c("A", "B"), font.label = list(size = 25))
    
ggsave(paste0(gpath, "Paper/Figures/combinedFigure1.png"), width = 20, height = 12)

#######################################
## Run model
#######################################

prop_forest_df = prop_forest_df %>% mutate(ECO_NAME = as.factor(ECO_NAME))

totalSpecRes = list()

for(i in unique(prop_forest_df$taxa)) {
    print(i)
    dat = prop_forest_df %>% filter(taxa == i)
    mod = bam(prop_forest ~ n_spec, data = dat, 
        family = binomial(link = "logit"), weights = n_spec, discrete = TRUE, nthreads = 6)
    totalSpecRes[[i]] = dat %>% mutate(res = residuals(mod))
}

dat = totalSpecRes %>% bind_rows()

## Make behrmann template raster for use in transforming in future
behr = "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
template_raster = rast() %>% project(behr)

## Load basic world outline from GSHHS to use in cropping/masking
coast = st_read("../Raw_Data/gshhg-shp-2.3.7/GSHHS_shp/c/GSHHS_c_L1.shp")
coast = st_transform(coast, behr)
coast_outline = st_cast(coast, "MULTILINESTRING")

coast_rast = coast %>% terra::rasterize(template_raster, touchess = TRUE, background = NA)

figure1_map_data = dplyr::select(dat, x, y, res, taxa) %>% group_split(taxa) %>% 
  map(., ~rast(dplyr::select(., x, y, res), crs = behr) %>% 
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


####

## Make amphibian map
amphib_map = ggplot() + 
            geom_sf(data = coast, fill = "grey70") + 
            geom_spatraster(data = figure1_map_data$Amphibians) + 
            geom_sf(data = coast_outline, linewidth = 0.75) +
            theme_void() + 
            scale_fill_continuous_diverging(na.value = "transparent", p1 = 0.1,
            rev = T) + 
            labs(fill = "Residuals") +
            # theme(legend.position = "none") +
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0)) +
            theme(plot.margin = unit(c(0.25,0,0.68,0), "cm"))

amphib_map

## Make bird map
bird_map = ggplot() + 
            geom_sf(data = coast, fill = "grey70") + 
            geom_spatraster(data = figure1_map_data$Birds) + 
            geom_sf(data = coast_outline, linewidth = 0.75) +
            theme_void() + 
            scale_fill_continuous_diverging(na.value = "transparent", p1 = 0.1,
            rev = T) + 
            labs(fill = "Residuals") +
            # theme(legend.position = "none") +
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0)) +
            theme(plot.margin = unit(c(0.25,0,0.68,0), "cm"))

bird_map

## Make mammal map
mammal_map = ggplot() + 
            geom_sf(data = coast, fill = "grey70") + 
            geom_spatraster(data = figure1_map_data$Mammals) + 
            geom_sf(data = coast_outline, linewidth = 0.75) +
            theme_void() + 
            scale_fill_continuous_diverging(na.value = "transparent", p1 = 0.1,
            rev = T) + 
            labs(fill = "Residuals") +
            # theme(legend.position = "none") +
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0)) +
            theme(plot.margin = unit(c(0.25,0,0.68,0), "cm"))

mammal_map

## Make reptile map
reptile_map = ggplot() + 
            geom_sf(data = coast, fill = "grey70") + 
            geom_spatraster(data = figure1_map_data$Reptiles) + 
            geom_sf(data = coast_outline, linewidth = 0.75) +
            theme_void() + 
            scale_fill_continuous_diverging(na.value = "transparent", p1 = 0.1,
            rev = T) + 
            labs(fill = "Residuals") +
            # theme(legend.position = "none") +
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0)) +
            theme(plot.margin = unit(c(0.25,0,0.68,0), "cm"))

reptile_map
