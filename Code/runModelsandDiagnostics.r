##########################################################
## Plot proportion of forest vs matrix species
## from latitudes and habitat usage taken from 
## IUCN for mammals, amphibs, and reptiles
## and from AVONET for birds
##########################################################

library(tidyverse)
library(janitor)
library(glmmTMB)
library(DHARMa)
library(piecewiseSEM)
library(ncf) ## correlogram
library(performance) ## check_collinearity
library(mgcv) ## bam

gpath = "/home/ben/Documents/PhD/PropForestSpecies/"
setwd(gpath)

prop_forest_df = read_csv(paste0(gpath, "Data/proportion_forest_species_analysis_data.csv"))

## Scale all continuous variables so their effect sizes are comparable
prop_forest_df = prop_forest_df %>% mutate(scaled_prop_forest_area = scale(prop_forest_area)[,1],
                                            scaled_prop_land_area_deforested = scale(prop_land_area_deforested)[,1],
                                            scaled_dist_equator_1000km = scale(dist_equator_1000km)[,1],
                                            scaled_geological_forest_time = scale(geological_forest_time)[,1],
                                            scaled_geological_forest_stability = scale(geological_forest_stability)[,1],
                                            scaled_alpha_plant_diversity= scale(alpha_plant_diversity)[,1],
                                            scaled_disturbances = scale(disturbances)[,1],
                                            ECO_NAME = as.factor(ECO_NAME), 
                                            taxa = as.factor(taxa)) 

##########################################################################################
## How much of the variation is explained by just distance to equator (latitude)?
##########################################################################################

## Save all model diagnostics to the same pdf
pdf(paste0(gpath, "Paper/Figures/Supplementary/modelDiagnostics.pdf"), width = 10, height = 10)

## Run model predicting prop forest species
## with latitude as the only predictor
## per taxa

latModels = list()

for(i in unique(prop_forest_df$taxa)) {
    dat = filter(prop_forest_df, taxa == i)
    mod = bam(prop_forest ~ scaled_dist_equator_1000km, weights = n_spec, data = dat, family = "binomial", discrete = TRUE, nthreads = 6)
    latModels[[i]] = mod
}

save(latModels, file = paste0(gpath, "Results/simpleLatModels.RData"))

## Get model diagnostics

for(i in names(latModels)) {

    mod = latModels[[i]]
    taxa = names(latModels[i])

    par(mfrow=c(2,2), oma = c(0, 0, 8, 0))
    testUniformity(mod, plot = T)
    testDispersion(mod, plot = T)
    testQuantiles(mod, plot = T)
    mtext(paste0("Diagnostics for latitude only model for ", i, " without accounting for spatial autocorrelation"), side = 3, line = 3, outer = TRUE, font = 2)

    ## Spatial autocorrelation using correlogram
    dat = filter(prop_forest_df, taxa == i)
    dat = mutate(dat, res = residuals(mod))
    correlog = with(dat, correlog(x, y, res, increment = 93694.99, resamp = 0, na.rm = TRUE))
    cor = as.data.frame(correlog[1:3]) %>% filter(n > 100 )
    plot(correlation ~ mean.of.class, data = cor, type = "o", xlab = "Distance (m)", ylab = "Correlation")

}

##########################################################################################
## Run model with all variables
##########################################################################################

fullModels = list()

for(i in unique(prop_forest_df$taxa)) {
    dat = filter(prop_forest_df, taxa == i)
    mod = bam(prop_forest ~ scaled_prop_forest_area*scaled_prop_land_area_deforested + scaled_disturbances + scaled_dist_equator_1000km + scaled_geological_forest_time + 
        scaled_geological_forest_stability + scaled_alpha_plant_diversity, weights = n_spec, data = dat, family = "binomial", discrete = TRUE, nthreads = 6)
    fullModels[[i]] = mod
}

save(fullModels, file = paste0(gpath, "Results/simpleFullModels.RData"))

## Get model diagnostics

for(i in names(fullModels)) {

    mod = fullModels[[i]]
    taxa = names(fullModels[i])

    par(mfrow=c(2,2), oma = c(0, 0, 8, 0))
    testUniformity(mod, plot = T)
    testDispersion(mod, plot = T)
    testQuantiles(mod, plot = T)

    col = as.data.frame(check_collinearity(mod)) %>%
    mutate(term = gsub("_", " ", Term)) %>%
    mutate(term = gsub("scaled", "", term)) %>%
    mutate(term = gsub(":", " x \n", term))

    mtext(paste0("Diagnostics for full model for ", i, " without accounting for spatial autocorrelation"), side = 3, line = 3, outer = TRUE, font = 2)

    ## Spatial autocorrelation using correlogram
    dat = filter(prop_forest_df, taxa == i)
    dat = mutate(dat, res = residuals(mod))
    correlog = with(dat, correlog(x, y, res, increment = 93694.99, resamp = 0, na.rm = TRUE))
    cor = as.data.frame(correlog[1:3]) %>% filter(n > 100 )
    plot(correlation ~ mean.of.class, data = cor, type = "o", xlab = "Distance (m)", ylab = "Correlation")

}

##########################################################################################
## Run latitude models accounting for spatial autocorrelation
##########################################################################################

randLatModels = list()

for(i in unique(prop_forest_df$taxa)) {
    dat = filter(prop_forest_df, taxa == i) %>% mutate(id = as.factor(1:nrow(.)))
    mod = bam(prop_forest ~ scaled_dist_equator_1000km + s(x, y, bs = "gp") + s(ECO_NAME, bs = "re"), 
    weights = n_spec, data = dat, family = binomial(), discrete = TRUE, nthreads = 6)
    randLatModels[[i]] = mod
}

save(randLatModels, file = paste0(gpath, "Results/randLatModels.RData"))

## Get model diagnostics

for(i in names(randLatModels)) {

    mod = randLatModels[[i]]
    taxa = names(randLatModels[i])

    par(mfrow=c(2,2), oma = c(0, 0, 8, 0))
    testUniformity(mod, plot = T)
    testDispersion(mod, plot = T)
    testQuantiles(mod, plot = T)
    mtext(paste0("Diagnostics for latitude only model for ", i, " accounting for spatial autocorrelation\nand including ecoregion as a random effect"), 
    side = 3, line = 3, outer = TRUE, font = 2)

    ## Spatial autocorrelation using correlogram
    dat = filter(prop_forest_df, taxa == i)
    dat = mutate(dat, res = residuals(mod))
    correlog = with(dat, correlog(x, y, res, increment = 93694.99, resamp = 0, na.rm = TRUE))
    cor = as.data.frame(correlog[1:3]) %>% filter(n > 100 )
    plot(correlation ~ mean.of.class, data = cor, type = "o", xlab = "Distance (m)", ylab = "Correlation")

}


##########################################################################################
## Run full models accounting for spatial autocorrelation
##########################################################################################

randFullModels = list()

for(i in unique(prop_forest_df$taxa)) {
    dat = filter(prop_forest_df, taxa == i)
    mod = bam(prop_forest ~ scaled_prop_forest_area*scaled_prop_land_area_deforested + scaled_disturbances + scaled_dist_equator_1000km + scaled_geological_forest_time + 
    scaled_geological_forest_stability + scaled_alpha_plant_diversity + s(x, y, bs = "gp") + s(ECO_NAME, bs = "re"),
    weights = n_spec, data = dat, family = binomial(), discrete = TRUE, nthreads = 6)
    randFullModels[[i]] = mod
}

save(randFullModels, file = paste0(gpath, "Results/randFullModels.RData"))

## Get model diagnostics

for(i in names(randFullModels)) {

    mod = randFullModels[[i]]
    taxa = names(randFullModels[i])

    par(mfrow=c(2,2), oma = c(0, 0, 8, 0))
    testUniformity(mod, plot = T)
    testDispersion(mod, plot = T)
    testQuantiles(mod, plot = T)

    mtext(paste0("Diagnostics for full model for ", i, " accounting for spatial autocorrelation\nand including ecoregion as a random effect"), 
    side = 3, line = 3, outer = TRUE, font = 2)

    ## Spatial autocorrelation using correlogram
    dat = filter(prop_forest_df, taxa == i)
    dat = mutate(dat, res = residuals(mod))
    correlog = with(dat, correlog(x, y, res, increment = 93694.99, resamp = 0, na.rm = TRUE))
    cor = as.data.frame(correlog[1:3]) %>% filter(n > 100 )
    plot(correlation ~ mean.of.class, data = cor, type = "o", xlab = "Distance (m)", ylab = "Correlation")

}

## Let's also run a model with all taxa (taxa as random intercept)

taxaLatModel = bam(prop_forest ~ scaled_dist_equator_1000km + s(x, y, bs = "gp") + 
    s(ECO_NAME, bs = "re") + s(taxa, bs = "re"),
    weights = n_spec, data = prop_forest_df, family = "binomial", discrete = TRUE, nthreads = 6)

save(taxaLatModel, file = paste0(gpath, "Results/taxaLatModel.RData"))

par(mfrow=c(2,2), oma = c(0, 0, 8, 0))
testUniformity(taxaLatModel, plot = T)
testDispersion(taxaLatModel, plot = T)
testQuantiles(taxaLatModel, plot = T)

mtext(paste0("Diagnostics for lat model for all taxa accounting for spatial autocorrelation\nand ecoregion as a random effect"), 
side = 3, line = 3, outer = TRUE, font = 2)

## Spatial autocorrelation using correlogram
prop_forest_df = mutate(prop_forest_df, res = residuals(taxaLatModel))
correlog = with(prop_forest_df, correlog(x, y, res, increment = 93694.99, resamp = 0, na.rm = TRUE))
cor = as.data.frame(correlog[1:3]) %>% filter(n > 100 )
plot(correlation ~ mean.of.class, data = cor, type = "o")

taxaFullModel = bam(prop_forest ~ scaled_prop_forest_area*scaled_prop_land_area_deforested + scaled_disturbances + 
    scaled_dist_equator_1000km + scaled_geological_forest_time + 
    scaled_geological_forest_stability + scaled_alpha_plant_diversity + 
    s(x, y, bs = "gp") + s(ECO_NAME, bs = "re") + s(taxa, bs = "re"),
    weights = n_spec, data = prop_forest_df, family = binomial(), discrete = TRUE, nthreads = 6)

save(taxaFullModel, file = paste0(gpath, "Results/taxaFullModel.RData"))

par(mfrow=c(2,2), oma = c(0, 0, 8, 0))
testUniformity(taxaFullModel, plot = T)
testDispersion(taxaFullModel, plot = T)
testQuantiles(taxaFullModel, plot = T)

mtext(paste0("Diagnostics for full model for all taxa accounting for spatial autocorrelation\nand ecoregions as a random effect"), 
side = 3, line = 3, outer = TRUE, font = 2)

## Spatial autocorrelation using correlogram
prop_forest_df = mutate(prop_forest_df, res = residuals(taxaFullModel))
correlog = with(prop_forest_df, correlog(x, y, res, increment = 93694.99, resamp = 0, na.rm = TRUE))
cor = as.data.frame(correlog[1:3]) %>% filter(n > 100 )
plot(correlation ~ mean.of.class, data = cor, type = "o")

## Test for variance inflation in a model without any interaction terms
mod = bam(prop_forest ~ scaled_prop_forest_area + scaled_prop_land_area_deforested + scaled_disturbances + 
    scaled_dist_equator_1000km + scaled_geological_forest_time + 
    scaled_geological_forest_stability + scaled_alpha_plant_diversity + 
    s(taxa, bs = "re"),
    weights = n_spec, data = prop_forest_df, family = binomial(), discrete = TRUE, nthreads = 6)

col = plot(check_collinearity(mod)) + 
theme_classic() + 
theme(text = element_text(size = 15),
axis.text.x = element_text(angle = 45, hjust = 1)) +
labs(title = "Variance inflation factor for all variables in a full model\nfor all taxa with no interaction terms",
subtitle = NULL)

print(col)

dev.off()
