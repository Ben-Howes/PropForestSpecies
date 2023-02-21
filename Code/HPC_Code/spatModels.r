##########################################################
## Plot proportion of forest vs matrix species
## from latitudes and habitat usage taken from 
## IUCN for mammals, amphibs, and reptiles
## and from AVONET for birds
##########################################################

library(tidyverse)
library(glmmTMB)

gpath = "/rds/general/user/bh719/home/prop_forest_species/"
setwd(gpath)

filtered_taxa = commandArgs(trailingOnly = TRUE)

filtered_taxa = case_when(filtered_taxa == 1 ~ "amphibians", 
                          filtered_taxa == 2 ~ "birds", 
                          filtered_taxa == 3 ~ "mammals",
                          filtered_taxa == 4 ~ "reptiles")

# Run model accounting for spatial autocorrelation
mod_data = read_csv(paste0(gpath, "Data/proportion_forest_species_analysis_data.csv"))

mod_data = mod_data %>% filter(taxa == filtered_taxa)

## Need to change the format a bit for spatial autocorrelation glmmTMB models
mod_data$pos = numFactor(scale(mod_data$x), scale(mod_data$y))
# then create a dummy group factor to be used as a random term
mod_data$ID = factor(rep(1, nrow(mod_data)))

mod = glmmTMB(prop_forest ~ prop_forest_area*historical_forest_loss + disturbances + dist_equator_1000km + geological_forest + mat(pos + 0 | ID), 
  data = mod_data, weights = n_spec, family = "binomial")

out = list(taxa = filtered_taxa, data = mod_data, model = mod)

save(out, file = paste0(gpath, "Results/", filtered_taxa,"_spat_mod.Rdata"))
