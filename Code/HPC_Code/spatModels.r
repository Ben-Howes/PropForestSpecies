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

size = case_when(run == 1 ~ 4031043, 
                        run == 2 ~ 2464440, 
                          run == 3 ~ 2745890,
                          run == 4 ~ 2735834)

# Run model accounting for spatial autocorrelation
load(paste0(gpath, "Data/sfDat.RData"))
mod_data = dat

mod_data = mod_data %>% mutate(scaled_prop_forest_area = scale(prop_forest_area),
                        scaled_prop_land_area_deforested = scale(prop_land_area_deforested),
                        scaled_dist_equator_1000km = scale(dist_equator_1000km),
                        scaled_geological_forest_time = scale(geological_forest_time),
                        scaled_geological_forest_stability = scale(geological_forest_stability),
                        scaled_alpha_plant_diversity= scale(alpha_plant_diversity),
                        scaled_disturbances = scale(disturbances))

mod_data = mod_data %>% filter(taxa == filtered_taxa)

mod_data = mutate(mod_data, id = as.factor(seq(1: nrow(mod_data))))
rook = dnearneigh(mod_data, d1 = 0, d2 = size, row.names = mod_data$id)
names(rook) = attr(rook, "region.id")

mod = bam(prop_forest ~ scaled_prop_forest_area*scaled_prop_land_area_deforested + scaled_disturbances + scaled_dist_equator_1000km + scaled_geological_forest_time + 
  scaled_geological_forest_stability + scaled_alpha_plant_diversity + s(id, bs = "mrf", xt = list(nb = rook)), 
  data = mod_data, family = binomial, weights = n_spec, discrete = TRUE, nthreads = 64)

save(mod, file = paste0(gpath, "Results/", filtered_taxa, "_spat_mod.RData"))
