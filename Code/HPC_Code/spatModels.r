##########################################################
## Plot proportion of forest vs matrix species
## from latitude only
## IUCN for mammals, amphibs, and reptiles
## and from AVONET for birds
##########################################################

library(dplyr)
library(mgcv)
library(spdep)

gpath = "/rds/general/user/bh719/home/prop_forest_species/"
setwd(gpath)

run = commandArgs(trailingOnly = TRUE)

filtered_taxa = case_when(run == 1 ~ "amphibians", 
                          run == 2 ~ "birds", 
                          run == 3 ~ "mammals",
                          run == 4 ~ "reptiles")


# Run model accounting for spatial autocorrelation
mod_data = read.csv(paste0(gpath, "Data/proportion_forest_species_analysis_data.csv"))

mod_data = mod_data %>% mutate(scaled_prop_forest_area = scale(prop_forest_area),
                        scaled_prop_land_area_deforested = scale(prop_land_area_deforested),
                        scaled_dist_equator_1000km = scale(dist_equator_1000km),
                        scaled_geological_forest_time = scale(geological_forest_time),
                        scaled_geological_forest_stability = scale(geological_forest_stability),
                        scaled_alpha_plant_diversity= scale(alpha_plant_diversity),
                        scaled_disturbances = scale(disturbances))

mod_data = mod_data %>% filter(taxa == filtered_taxa) %>% mutate(id = as.factor(1:nrow(.)))

mod = bam(prop_forest ~ scaled_prop_forest_area*scaled_prop_land_area_deforested + scaled_disturbances + scaled_dist_equator_1000km + scaled_geological_forest_time + 
  scaled_geological_forest_stability + scaled_alpha_plant_diversity + s(x, y, bs = "gp") + s(id, bs = "re"), 
  data = mod_data, family = binomial, weights = n_spec, discrete = TRUE, nthreads = 64)

save(mod, file = paste0(gpath, "Results/", filtered_taxa, "_spat_mod.RData"))
