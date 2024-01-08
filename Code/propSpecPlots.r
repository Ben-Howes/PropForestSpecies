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
    geom_point(size = 2, alpha = 0.3) +
    theme_classic() +
    labs(x = "Absolute Latitude", y = NULL, col = NULL) +
    theme(legend.position = c(0.3885, 0.92),
    legend.spacing.y = unit(0.25, 'cm'),
    legend.box.background = element_rect(colour = "black", 
    linewidth = 1.1), legend.margin = margin(t=0,r=20,b=10,l=15)) +
    theme(text = element_text(size = 30), legend.text=element_text(size = 20)) + 
    guides(col = guide_legend(byrow = TRUE)) +
    facet_rep_wrap(. ~ taxa, scales = "free_x") +
    geom_smooth(aes(group = prop_type), linewidth = 3.5, col = "black", method = "gam", se = FALSE) +
    geom_smooth(linewidth = 2.5, method = "gam", se = FALSE) +
    scale_y_continuous(expand = c(0, 0), limits = c(-0.01, 1)) +
    scale_x_continuous(expand = c(0, 0), limits = c(-0.01, 71)) +
    scale_color_manual(values = c("#3b528b", "#f8e621"), labels = c("Proportion\nForest Species", "Proportion Max\nSpecies Richness"))

ggsave(paste0(gpath, "Paper/Figures/latitudeDotPlot.png"), plot = latDotPlot, width = 20, height = 12)

## Plot of proportion of forest species against species richness
propVsRich = ggplot(prop_forest_df, aes(n_spec, prop_forest)) +
    geom_point(size = 2, alpha = 0.3, col = "#3b528b") +
    theme_classic() +
    labs(x = "Total Species Richness", y = "Proportion Forest Species") +
    theme(text = element_text(size = 30)) +
    facet_rep_wrap(. ~ taxa, scales = "free_x") +
    geom_smooth(linewidth = 3.5, col = "black", method = "gam", se = F) +
    geom_smooth(linewidth = 2.5, col = "#3b528b", method = "gam") +
    scale_y_continuous(expand = c(0, 0), limits = c(-0.01, 1)) +
    scale_x_continuous(expand = c(0.01, 0)) 

ggsave(paste0(gpath, "Paper/Figures/propVsRich.png"), plot = propVsRich, width = 20, height = 12)

ggarrange(latDotPlot + theme(legend.position = c(0.315, 0.935), legend.text = element_text(size = 15)), 
    propVsRich + labs(y = NULL), ncol = 2, nrow = 1)
    
ggsave(paste0(gpath, "Paper/Figures/combinedFigure1.png"), width = 20, height = 12)
