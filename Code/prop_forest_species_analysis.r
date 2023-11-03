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
library(piecewiseSEM)
library(gstat) # variogram
library(sp) ## spatialDataframe
library(nlme) ## corSpatial
library(MASS) ## glmmPQL

gpath = "/home/ben/Documents/PhD/PropForestSpecies/"
setwd(gpath)

prop_forest_df = read_csv(paste0(gpath, "Data/proportion_forest_species_analysis_data.csv"))

## Scale all continuous variables so their effect sizes are comparable
prop_forest_df = prop_forest_df %>% mutate(scaled_prop_forest_area = scale(prop_forest_area),
                                            scaled_prop_land_area_deforested = scale(prop_land_area_deforested),
                                            scaled_dist_equator_1000km = scale(dist_equator_1000km),
                                            scaled_geological_forest_time = scale(geological_forest_time),
                                            scaled_geological_forest_stability = scale(geological_forest_stability),
                                            scaled_alpha_plant_diversity= scale(alpha_plant_diversity),
                                            scaled_disturbances = scale(disturbances))

##########################################################################################
## How much of the variation is explained by just distance to equator (latitude)?
##########################################################################################

## Create function to run models with and without random effect for autocorrelation
## and save DHARMA results, as well as model outputs, R2, and the models of course

runLatMods = function(x) {

    taxa = x$taxa[[1]]

    pdf(paste0(gpath, "Paper/Figures/Supplementary/", taxa, "LatDiagnosticPlots.pdf"),width = 14)

    mod = glmmTMB(prop_forest ~ scaled_dist_equator_1000km, weights = n_spec, data = x, family = "binomial")

    par()
    testResiduals(mod, plot = T)
    testQuantiles(mod, plot = T)
    testZeroInflation(mod, plot = T)

    ## Test for spatial autocorrelation using semivariogram
    dat = mutate(x, res = residuals(mod))
    coordinates(dat)=~x+y ## Convert to spatialDataFrame
    ## Create semi-variance data, looking at pixels up to 5000km apart, and in 20 bins
    vario = variogram(res ~ x+y, data = dat, cutoff = 5000000, width = 5000000/20, cloud = F)

    randMod = glmmTMB(prop_forest ~ scaled_dist_equator_1000km + (1|ECO_NAME),
        weights = n_spec, data = x, family = "binomial")

    testResiduals(randMod, plot = T)
    testQuantiles(randMod, plot = T)
    testZeroInflation(randMod, plot = T)

    ## Test for spatial autocorrelation using semivariogram
    dat = mutate(x, res = residuals(randMod))
    coordinates(dat)=~x+y ## Convert to spatialDataFrame
    ## Create semi-variance data, looking at pixels up to 5000km apart, and in 20 bins
    spatVario = variogram(res ~ x+y, data = dat, cutoff = 5000000, width = 5000000/20, cloud = F)

    datVario = vario %>%
        dplyr::select(np, dist, gamma) %>%
        mutate(type = "No Random Effect") %>%
        bind_rows(mutate(dplyr::select(spatVario, np, dist, gamma), type = "Random Effect"))
    
    varioPlot = ggplot(varioDat, aes(dist/1000, gamma, col = type)) +
    geom_point(size = 5) +
    geom_smooth(linewidth = 3) +
    theme_classic() +
    labs(x = "Distance (km)", y = "Semi-variance", col = NULL) +
    theme(text = element_text(size = 30)) +
    scale_colour_viridis_d(labels = c("No Random Effect", "Random Effect"))

    print(varioPlot)
    dev.off()

    save(mod, randMod, file = paste0(gpath, "Results/", taxa, "LatMods.Rdata"))


}

runLatMods(filter(prop_forest_df, taxa == "amphibians"))

amphibian_lat_mod = glmmTMB(prop_forest ~ scaled_dist_equator_1000km + (1|ECO_NAME), weights = n_spec, data = filter(prop_forest_df, taxa == "amphibians"), family = "binomial")
bird_lat_mod = glmmTMB(prop_forest ~ scaled_dist_equator_1000km, weights = n_spec, data = filter(prop_forest_df, taxa == "birds"), family = "binomial")
mammal_lat_mod = glmmTMB(prop_forest ~ scaled_dist_equator_1000km, weights = n_spec, data = filter(prop_forest_df, taxa == "mammals"), family = "binomial")
reptile_lat_mod = glmmTMB(prop_forest ~ scaled_dist_equator_1000km , weights = n_spec, data = filter(prop_forest_df, taxa == "reptiles"), family = "binomial")

##########################################################################################
## Run model with all variables
## Using econame as a random effect to account for spatial autocorrelation
##########################################################################################

################
## Amphibians
################

## Without accounting for spatial autocorrelation

amphibian_data = prop_forest_df %>% filter(taxa == "amphibians")

amphibian_mod = glmmTMB(prop_forest ~ scaled_prop_forest_area*scaled_prop_land_area_deforested + scaled_disturbances + scaled_dist_equator_1000km + scaled_geological_forest_time + 
    scaled_geological_forest_stability + scaled_alpha_plant_diversity + (1|ECO_NAME),
    zi = ~scaled_prop_forest_area,  
    weights = n_spec, data = amphibian_data, family = "binomial")

amphibian_data = amphibian_data %>% mutate(res = residuals(amphibian_mod))
correlog1 = with(dplyr::select(amphibian_data, x, y, res), correlog(x, y, res, increment = 93694.99, resamp = 0, na.rm = TRUE))

plot(correlog)
lines(correlation ~ mean.of.class, data=data.frame(correlog1[1:3]), type='o', col='red')

## Test model diagnostics
## Note these are pretty bad as we aren't accounting for spatial autocorrelation
## Accounting for this solves all of the issues in my experience
testResiduals(amphibian_mod, plot = T) ## Awful, disperson issues, outliers, lack of normality
testQuantiles(amphibian_mod, plot = T) ## Heteroscedasticity issues
testZeroInflation(amphibian_mod, plot = T) ## Zero inflation issues
testSpatialAutocorrelation(amphibian_mod, x = amphibian_data$x, y = amphibian_data$y)

## Test for spatial autocorrelation using semivariogram
amphibian_data = mutate(amphibian_data, modRes = residuals(amphibian_mod))
coordinates(amphibian_data)=~x+y ## Convert to spatialDataFrame
## Create semi-variance data, looking at pixels up to 5000km apart, and in 20 bins
amphibianVario = variogram(modRes ~ x+y, data = amphibian_data, cutoff = 5000000, width = 5000000/20, cloud = F)
plot(amphibianVario) ## Autocorrelation up to 4e+06, x values from 0.002 - 0.007
## Fit a variogram
## Fit a semi-variogram model
varioFit = fit.variogram(amphibianVario, model = vgm(psill = 2, model="Sph", range = 1000000, nugget = 0.5)) 
## Look at the result of the fit
varioFit
plot(amphibianVario, varioFit)

## Make correlation structure to account for the autocorrelation
correl = corSpatial(value = c(400000), form = ~x + y, nugget = FALSE, fixed = FALSE, type = "spherical")

## Account for spatial autocorrelation using the ecoregion name as a random effect
amphibian_data = prop_forest_df %>% filter(taxa == "amphibians")

amphibian_spat_mod = glmmTMB(prop_forest ~ scaled_prop_forest_area*scaled_prop_land_area_deforested + 
    scaled_disturbances + scaled_dist_equator_1000km + scaled_geological_forest_time + 
    scaled_geological_forest_stability + scaled_alpha_plant_diversity + (scaled_prop_forest_area*scaled_prop_land_area_deforested|ECO_NAME),
    zi = ~scaled_prop_forest_area, 
    weights = n_spec, data = amphibian_data, family = "binomial")

## Test model diagnostics
## Note these are pretty bad as we aren't accounting for spatial autocorrelation
## Accounting for this solves all of the issues in my experience
testResiduals(amphibian_spat_mod, plot = T) ## Much better, no dispersion issues and outliers
## failed KS test but not a big problem
testQuantiles(amphibian_spat_mod, plot = T) ## Quantiles look pretty good
testZeroInflation(amphibian_spat_mod, plot = T) ## Mild zero inflation issues, not a big problem
res2 = recalculateResiduals(amphibian_spat_mod, rotation = TRUE)
testSpatialAutocorrelation(res2, x = amphibian_data$x, y = amphibian_data$y)

## Test spatial autocorrelation in a semivariogram
amphibian_data = mutate(amphibian_data, res = residuals(amphibian_spat_mod))
coordinates(amphibian_data)=~x+y ## Convert to spatialDataFrame
## Create vario data
amphibianSpatVario = variogram(res ~ x+y, data = amphibian_data, cutoff = 5000000, width = 5000000/20, cloud = F)
plot(vario) ## Much less spatial autocorrelation, x axis from 0.001 - 0.0015

summary(amphibian_spat_mod)

## Plot varigrams with and without random effect
varioDat = amphibianVario %>%
    dplyr::select(np, dist, gamma) %>%
    mutate(type = "No Random Effect") %>%
    bind_rows(mutate(dplyr::select(amphibianSpatVario, np, dist, gamma), type = "Random Effect"))

ggplot(varioDat, aes(dist/1000, gamma, col = type)) +
    geom_point(size = 5) +
    geom_smooth(linewidth = 3) +
    theme_classic() +
    labs(x = "Distance (km)", y = "Semi-variance", col = NULL) +
    theme(text = element_text(size = 30)) +
    scale_colour_viridis_d(labels = c("No Random Effect", "Random Effect"))


### glmmPQL
library(MASS)
amphibian_data = prop_forest_df %>% filter(taxa == "amphibians")
amphibian_data = mutate(amphibian_data, dummy = 1)

test = glmmPQL(prop_forest ~ scaled_prop_forest_area*scaled_prop_land_area_deforested + 
    scaled_disturbances + scaled_dist_equator_1000km + scaled_geological_forest_time + 
    scaled_geological_forest_stability + scaled_alpha_plant_diversity, 
    random = list(dummy = ~1),
    correlation = correl,
    weights = n_spec, data = amphibian_data, family = "binomial")

testSpatialAutocorrelation(test, x = amphibian_data$x, y = amphibian_data$y)

amphibian_data = mutate(amphibian_data, res = resid(test, type = "normalized"))
coordinates(amphibian_data)=~x+y ## Convert to spatialDataFrame
## Create vario data
amphibianSpatVario = variogram(res ~ x+y, data = amphibian_data, cutoff = 5000000, width = 5000000/20, cloud = F)
plot(amphibianSpatVario) ## Much less spatial autocorrelation, x axis from 0.001 - 0.0015

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
