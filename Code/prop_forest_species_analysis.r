##########################################################
## Plot proportion of forest vs matrix species
## from latitudes and habitat usage taken from 
## IUCN for mammals, amphibs, and reptiles
## and from AVONET for birds
##########################################################

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
library(glmmTMB)
library(DHARMa)
library(dotwhisker)

gpath = "/home/ben/Documents/PhD/PropForestSpecies/"
setwd(gpath)

prop_forest_df = read_csv(paste0(gpath, "Data/proportion_forest_species_analysis_data.csv"))

## Scale all continuous variables so their effect sizes are comparable
prop_forest_df = prop_forest_df %>% mutate(scaled_prop_forest_area = scale(prop_forest_area),
                                            scaled_historical_forest_loss = scale(historical_forest_loss),
                                            scaled_dist_equator_1000km = scale(dist_equator_1000km),
                                            scaled_geological_forest = scale(geological_forest),
                                            scaled_disturbances = scale(disturbances))

##########################################################################################
## How much of the variation is explained by just distance to equator (latitude)?
##########################################################################################

lat_mod = glmmTMB(prop_forest ~ scaled_dist_equator_1000km + (scaled_dist_equator_1000km|taxa), zi = ~scaled_dist_equator_1000km, weights = n_spec, data = prop_forest_df, family = "binomial")
amphibian_lat_mod = glmmTMB(prop_forest ~ scaled_dist_equator_1000km, weights = n_spec, data = filter(prop_forest_df, taxa == "amphibians"), family = "binomial")
bird_lat_mod = glmmTMB(prop_forest ~ scaled_dist_equator_1000km, weights = n_spec, data = filter(prop_forest_df, taxa == "birds"), family = "binomial")
mammal_lat_mod = glmmTMB(prop_forest ~ scaled_dist_equator_1000km, weights = n_spec, data = filter(prop_forest_df, taxa == "mammals"), family = "binomial")
reptile_lat_mod = glmmTMB(prop_forest ~ scaled_dist_equator_1000km , weights = n_spec, data = filter(prop_forest_df, taxa == "reptiles"), family = "binomial")

## Check diagnostics of random effect model, others likely similar
## We saw zero inflation issues so adding zi call
## Probably still issues with spatial autocorrelation, but diagnostics not looking too bad considering
testResiduals(lat_mod, plot = T)
simulateResiduals(lat_mod, plot = T)
testZeroInflation(lat_mod, plot = T)

## Look at results
summary(lat_mod)
r.squaredGLMM(lat_mod)

## All taxa follow the same negative pattern
ranef(lat_mod)

##########################################################################################
## Run model with all variables
##########################################################################################

## Run full model
## The most complex structure I can use without getting convergence or singularity warnings includes prop forest area and historical forest area as random slopes
full_mod = glmmTMB(prop_forest ~ scaled_prop_forest_area*scaled_historical_forest_loss + scaled_disturbances + scaled_dist_equator_1000km + scaled_geological_forest + 
  (scaled_prop_forest_area|taxa), zi = ~scaled_prop_forest_area,
  weights = n_spec, data = prop_forest_df, family = "binomial")

## Diagnostics
testResiduals(full_mod, plot = T)
simulateResiduals(full_mod, plot = T)
testZeroInflation(full_mod, plot = T)

summary(full_mod)

## I have successfully run a random effect binomial model, with all variables are predictors and then random intercepts for the taxa
## This however does not look great when assessing diagnostics with DHARMa, with a key issue being spatial autocorrelation - tested on individually models
## I can run spatial autocorrelation models to account for this - from running sample data they seem to fix nearly all diagnostic issues
## But this will need to be done on the HPC as the computing time is too long
## I will run models for each taxa individually

#################################################################################
## For now without spatial autocorrelation 
## We run basic models without any spatial autocorrelation here
################################################################################

################
## Amphibians
################
amphibian_data = prop_forest_df %>% filter(taxa == "amphibians")

amphibian_mod = glmmTMB(prop_forest ~ scaled_prop_forest_area*scaled_historical_forest_loss + scaled_disturbances + scaled_dist_equator_1000km + scaled_geological_forest, 
    weights = n_spec, data = amphibian_data, family = "binomial")

## Test model diagnostics
## Note these are pretty bad as we aren't accounting for spatial autocorrelation
## Accounting for this solves all of the issues in my experience
testResiduals(amphibian_mod, plot = T)
simulateResiduals(amphibian_mod, plot = T)
testZeroInflation(amphibian_mod, plot = T)

## Model summary
summary(amphibian_mod)

################
## Reptiles
################
reptile_data = prop_forest_df %>% filter(taxa == "reptiles")

reptile_mod = glmmTMB(prop_forest ~ scaled_prop_forest_area*scaled_historical_forest_loss + scaled_disturbances + scaled_dist_equator_1000km + scaled_geological_forest, 
    weights = n_spec, data = reptile_data, family = "binomial")

## Test model diagnostics
## Note these are pretty bad as we aren't accounting for spatial autocorrelation
## Accounting for this solves all of the issues in my experience
testResiduals(reptile_mod, plot = T)
simulateResiduals(reptile_mod, plot = T)
testZeroInflation(reptile_mod, plot = T)

## Model summary
summary(reptile_data)

################
## Mammals
################
mammal_data = prop_forest_df %>% filter(taxa == "mammals")

mammal_mod = glmmTMB(prop_forest ~ scaled_prop_forest_area*scaled_historical_forest_loss + scaled_disturbances + scaled_dist_equator_1000km + scaled_geological_forest, 
    weights = n_spec, data = mammal_data, family = "binomial")

## Test model diagnostics
## Note these are pretty bad as we aren't accounting for spatial autocorrelation
## Accounting for this solves all of the issues in my experience
testResiduals(mammal_mod, plot = T)
simulateResiduals(mammal_mod, plot = T)
testZeroInflation(mammal_mod, plot = T)

## Model summary
summary(mammal_mod)

################
## Birds
################
bird_data = prop_forest_df %>% filter(taxa == "birds")

bird_mod = glmmTMB(prop_forest ~ scaled_prop_forest_area*scaled_historical_forest_loss + scaled_disturbances + scaled_dist_equator_1000km + scaled_geological_forest, 
    weights = n_spec, data = bird_data, family = "binomial")

## Test model diagnostics
## Note these are pretty bad as we aren't accounting for spatial autocorrelation
## Accounting for this solves all of the issues in my experience
testResiduals(bird_mod, plot = T)
simulateResiduals(bird_mod, plot = T)
testZeroInflation(bird_mod, plot = T)

## Model summary
summary(bird_mod)

## Plot model outputs for all taxa models
dwplot(list(amphibian_mod, bird_mod, mammal_mod, reptile_mod), 
whisker_args = list(size = 1.5),
dot_args = list(size = 5),
model_order = c("Model 1", "Model 2", "Model 3", "Model 4")) %>% 
relabel_predictors(
    c(
            scaled_prop_forest_area = "Proporion Forested Area",
            scaled_historical_forest_loss = "Proportion of Historial Forest Loss",
            "scaled_prop_forest_area:scaled_historical_forest_loss" = 
            "Interaction Term between Proportion Forested Area\nand Hostiral Fores Loss",
            scaled_disturbances = "Historical Disturbances",
            scaled_geological_forest = "Geological Time Area has been Forest",
            scaled_dist_equator_1000km = "Distance to Equator (1000km)"
        )
    ) + 
  theme_classic() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(text = element_text(size = 20)) +
  labs(x = "Coefficient of Effect on\nProportion of Forest Species (95% CI)", y = "Scaled Variables", col = "Taxa") +
  scale_colour_viridis_d(labels = c("Reptiles", "Mammals", "Birds", "Amphibians"))


# y_partial = remef(full_mod, fix = c("prop_forest_area", "disturbances", "dist_equator_1000km"), ran = "all")
# test = prop_forest_df %>% add_column(partial_historic_res = y_partial)

## Interaction between forest area and historia forest loss: 
## If you used to have lots of forest then the effect of losing forest is greater
## E.g if you have 0.25 forest now, and you used to have 0.5 (+0.25) then the difference will be smaller than
## if you have 0.5 forest now and you used to have 0.75 (+0.25)
## Interaction between forest area and disturbance:
## Disturbance has a greater effect at low levels of forest area, and little effect at high levels
## Interaction between forest area and distance from equator:
## Forest area has a greater effect far from the equator as opposed to close to the equator when it has a smaller effect

# ggplot(test) + geom_raster(aes(x, y, fill = partial_historic_res)) + theme_classic() +
#   scale_fill_continuous_diverging(mid = 0, palette = "Purple-Green", guide = "colourbar", na.value = "black") + 
#   labs(x = NULL, y = NULL) +
#   theme(axis.title = element_blank(),axis.text = element_blank(),axis.ticks.x = element_blank(),
#         axis.line=element_blank()) + facet_rep_wrap(. ~ taxa, ncol = 2)
