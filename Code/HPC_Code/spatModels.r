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

mod_data = mod_data %>% mutate(prop_forest_area = scale(prop_forest_area),
                                prop_land_area_deforested = scale(prop_land_area_deforested),
                                dist_equator_1000km = scale(dist_equator_1000km),
                                geological_forest = scale(geological_forest),
                                disturbances = scale(disturbances))

mod_data = mod_data %>% filter(taxa == filtered_taxa) %>% slice_sample(prop = 0.1, replace = TRUE)

## Need to change the format a bit for spatial autocorrelation glmmTMB models
mod_data$pos = numFactor(scale(mod_data$x), scale(mod_data$y))
# then create a dummy group factor to be used as a random term
mod_data$ID = factor(rep(1, nrow(mod_data)))

mod = glmmTMB(prop_forest ~ prop_forest_area*prop_land_area_deforested + disturbances + dist_equator_1000km + geological_forest + mat(pos + 0 | ID), 
  data = mod_data, weights = n_spec, family = "binomial")

## Get coefficients
coef = summary(mod)$coef$cond %>%
        as.data.frame() %>%
        rename("est" = 1, "std_err" = 2, "z_val" = 3, "p_val" = 4)

CI = confint(mod, parm = rownames(coef)) %>% 
      as.data.frame() %>% 
      rename("low_2.5" = 1,"up_97.5" = 2, "est" = 3) %>%
      dplyr::select(-est)

final = data.frame(predictor = rownames(coef), coef, CI, taxa = filtered_taxa, run = run)

write_csv(final, paste0(gpath, "Results/", filtered_taxa, run, "_spat_mod.csv"))
