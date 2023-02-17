##########################################################
## Plot proportion of forest vs matrix species
## from latitudes and habitat usage taken from 
## IUCN for mammals, amphibs, and reptiles
## and from AVONET for birds
##########################################################

install.packages("spaMM", repos='http://cran.us.r-project.org')

library(tidyverse)
library(spaMM)

gpath = "/rds/general/user/bh719/home/prop_forest_species/"
setwd(gpath)

# Run model accounting for spatial autocorrelation
mod_data = read_csv(paste0(gpath, "Data/spaMM_data.csv"))

spat_forest_amount_mod = fitme(cbind(n_forest, n_nonforest) ~ prop_forest_area*historical_forest_loss + disturbances + dist_equator_1000km + geological_forest + Matern(1|x+y) +
  (prop_forest_area*historical_forest_loss|taxa) + (disturbances|taxa) + (dist_equator_1000km|taxa) + (geological_forest|taxa), 
  data = mod_data, family = "binomial")

save(spat_forest_amount_mod, file = paste0(gpath, "Results/spaMM_model.Rdata"))