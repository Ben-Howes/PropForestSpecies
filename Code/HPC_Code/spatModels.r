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

run = commandArgs(trailingOnly = TRUE)

filtered_taxa = case_when(run %in% seq(1, 1000, by = 1) ~ "amphibians", 
                          run %in% seq(1001, 2000, by = 1) ~ "birds", 
                          run %in% seq(2001, 3000, by = 1) ~ "mammals",
                          run %in% seq(3001, 4000, by = 1) ~ "reptiles")

# Run model accounting for spatial autocorrelation
mod_data = read_csv(paste0(gpath, "Data/proportion_forest_species_analysis_data.csv"))

mod_data = mod_data %>% filter(taxa == filtered_taxa) %>% slice_sample(prop = 0.1)

## Need to change the format a bit for spatial autocorrelation glmmTMB models
mod_data$pos = numFactor(scale(mod_data$x), scale(mod_data$y))
# then create a dummy group factor to be used as a random term
mod_data$ID = factor(rep(1, nrow(mod_data)))

mod = glmmTMB(prop_forest ~ prop_forest_area*historical_forest_loss + disturbances + dist_equator_1000km + geological_forest + mat(pos + 0 | ID), 
  data = mod_data, weights = n_spec, family = "binomial")

out = list(PBS = run, taxa = filtered_taxa, data = mod_data, model = mod)

save(out, file = paste0(gpath, "Results/", filtered_taxa, run, "_spat_mod.Rdata"))
