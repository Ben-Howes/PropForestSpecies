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
                                            scaled_disturbances = scale(disturbances)[,1])

################################# TESTING
library(mgcv)

mod_data = prop_forest_df %>% filter(taxa == "mammals") %>% slice_sample(prop = 1, replace = FALSE)
mod_data = mod_data %>% mutate(ECO_NAME = as.factor(ECO_NAME))

mod = bam(prop_forest ~ scaled_prop_forest_area*scaled_prop_land_area_deforested + scaled_disturbances + scaled_dist_equator_1000km + scaled_geological_forest_time + 
    scaled_geological_forest_stability + scaled_alpha_plant_diversity + s(x, y, bs = "gp") + s(ECO_NAME, bs = "re"),
    data = mod_data, family = binomial, weights = n_spec, discrete = TRUE, nthreads = 6)

testQuantiles(mod, plot = T)
testResiduals(mod, plot = T)
testZeroInflation(mod, plot = T)

mod_data = mutate(mod_data, res = residuals(mod))
correlog = with(mod_data, correlog(x, y, res, increment = 93694.99, resamp = 0, na.rm = TRUE))
cor = as.data.frame(correlog[1:3]) %>% filter(n > 100)
plot(correlation ~ mean.of.class, data = cor, type = "o")

##########################################################################################
## How much of the variation is explained by just distance to equator (latitude)?
##########################################################################################

## Run model predicting prop forest species
## with lattiude as the only predictor
## per taxa

latModels = list()

for(i in unique(prop_forest_df$taxa)) {
    dat = filter(prop_forest_df, taxa == i)
    mod = glmmTMB(prop_forest ~ scaled_dist_equator_1000km, weights = n_spec, data = dat, family = "binomial")
    latModels[[i]] = mod
}

save(latModels, file = paste0(gpath, "Results/simpleLatModels.RData"))

## Get model diagnostics

pdf(paste0(gpath, "Paper/Figures/Supplementary/latModelDiagnostics.pdf"), width = 10, height = 10)

for(i in names(latModels)) {

    mod = latModels[[i]]
    taxa = names(latModels[i])

    par(mfrow=c(2,3), oma = c(0, 0, 8, 0))
    testUniformity(mod, plot = T)
    testDispersion(mod, plot = T)
    plot(fitted(mod), residuals(mod, type = "pearson"))
    testZeroInflation(mod, plot = T)
    mtext(paste0("Diagnostics for latitude only model for ", i, " without accounting for spatial autocorrelation"), side = 3, line = 3, outer = TRUE, font = 2)

    ## Spatial autocorrelation using correlogram
    dat = filter(prop_forest_df, taxa == i)
    dat = mutate(dat, res = residuals(mod))
    correlog = with(dat, correlog(x, y, res, increment = 93694.99, resamp = 0, na.rm = TRUE))
    cor = as.data.frame(correlog[1:3]) %>% filter(n > 100 )
    plot(correlation ~ mean.of.class, data = cor, type = "o")

}

dev.off()

##########################################################################################
## Run model with all variables
##########################################################################################

fullModels = list()

for(i in unique(prop_forest_df$taxa)) {
    dat = filter(prop_forest_df, taxa == i)
    mod = glmmTMB(prop_forest ~ scaled_prop_forest_area*scaled_prop_land_area_deforested + scaled_disturbances + scaled_dist_equator_1000km + scaled_geological_forest_time + 
        scaled_geological_forest_stability + scaled_alpha_plant_diversity, weights = n_spec, data = dat, family = "binomial")
    fullModels[[i]] = mod
}

save(fullModels, file = paste0(gpath, "Results/simpleFullModels.RData"))

## Get model diagnostics

pdf(paste0(gpath, "Paper/Figures/Supplementary/fullModelDiagnostics.pdf"), width = 10, height = 10)

for(i in names(fullModels)) {

    mod = fullModels[[i]]
    taxa = names(fullModels[i])

    par(mfrow=c(2,3), oma = c(0, 0, 8, 0))
    testUniformity(mod, plot = T)
    testDispersion(mod, plot = T)
    plot(fitted(mod), residuals(mod, type = "pearson"))
    testZeroInflation(mod, plot = T)

    col = as.data.frame(check_collinearity(mod)) %>%
    mutate(term = gsub("_", " ", Term)) %>%
    mutate(term = gsub("scaled", "", term)) %>%
    mutate(term = gsub(":", " x \n", term))

    plot.new()
    plot.window(xlim = c(1, nrow(col)), ylim = c(floor(min(col$VIF_CI_low)), ceiling(max(col$VIF_CI_high))))

    segments(1:nrow(col), col$VIF_CI_low, 1:nrow(col), col$VIF_CI_high)
    points(1:nrow(col), col$VIF, cex = 1.5, col = "black", pch = 19)

    axis(1, at = 1:nrow(col), labels = FALSE)
    axis(2, las = 2)

    title("Collinearity Check", adj=0)
    title(xlab = "Model Terms", ylab = "Variance Inflation Factor (VIF)")

    mtext(paste0("Diagnostics for full model for ", i, " without accounting for spatial autocorrelation"), side = 3, line = 3, outer = TRUE, font = 2)

    ## Spatial autocorrelation using correlogram
    dat = filter(prop_forest_df, taxa == i)
    dat = mutate(dat, res = residuals(mod))
    correlog = with(dat, correlog(x, y, res, increment = 93694.99, resamp = 0, na.rm = TRUE))
    cor = as.data.frame(correlog[1:3]) %>% filter(n > 100 )
    plot(correlation ~ mean.of.class, data = cor, type = "o")

}

dev.off()

##########################################################################################
## Run latitude models accounting for spatial autocorrelation
##########################################################################################

randLatModels = list()

for(i in unique(prop_forest_df$taxa)) {
    dat = filter(prop_forest_df, taxa == i) %>% mutate(id = 1:nrow(.))
    mod = glmmTMB(prop_forest ~ scaled_dist_equator_1000km + (1|ECO_NAME) + (1|id), 
        zi = ~scaled_dist_equator_1000km, 
        weights = n_spec, data = dat, family = "binomial")
    randLatModels[[i]] = mod
}

save(randLatModels, file = paste0(gpath, "Results/randLatModels.RData"))

## Get model diagnostics

pdf(paste0(gpath, "Paper/Figures/Supplementary/randLatModelDiagnostics.pdf"), width = 10, height = 10)

for(i in names(randLatModels)) {

    mod = randLatModels[[i]]
    taxa = names(randLatModels[i])

    par(mfrow=c(2,3), oma = c(0, 0, 8, 0))
    testUniformity(mod, plot = T)
    testDispersion(mod, plot = T)
    plot(fitted(mod), residuals(mod, type = "pearson"))
    testZeroInflation(mod, plot = T)
    mtext(paste0("Diagnostics for latitude only model for ", i, " accounting for spatial autocorrelation"), side = 3, line = 3, outer = TRUE, font = 2)

    ## Spatial autocorrelation using correlogram
    dat = filter(prop_forest_df, taxa == i)
    dat = mutate(dat, res = residuals(mod))
    correlog = with(dat, correlog(x, y, res, increment = 93694.99, resamp = 0, na.rm = TRUE))
    cor = as.data.frame(correlog[1:3]) %>% filter(n > 100 )
    plot(correlation ~ mean.of.class, data = cor, type = "o")

}

dev.off()

##########################################################################################
## Run full models accounting for spatial autocorrelation
##########################################################################################

randFullModels = list()

for(i in unique(prop_forest_df$taxa)) {
    dat = filter(prop_forest_df, taxa == i)
    mod = glmmTMB(prop_forest ~ scaled_prop_forest_area*scaled_prop_land_area_deforested + scaled_disturbances + scaled_dist_equator_1000km + scaled_geological_forest_time + 
        scaled_geological_forest_stability + scaled_alpha_plant_diversity + (1|ECO_NAME),
    zi = ~scaled_prop_forest_area + scaled_prop_land_area_deforested,
    weights = n_spec, data = dat, family = binomial())
    randFullModels[[i]] = mod
}

save(randFullModels, file = paste0(gpath, "Results/randFullModels.RData"))

## Get model diagnostics

pdf(paste0(gpath, "Paper/Figures/Supplementary/randFullModelsDiagnostics.pdf"), width = 10, height = 10)

for(i in names(randFullModels)) {

    mod = randFullModels[[i]]
    taxa = names(randFullModels[i])

    par(mfrow=c(2,3), oma = c(0, 0, 8, 0))
    testUniformity(mod, plot = T)
    testDispersion(mod, plot = T)
    plot(fitted(mod), residuals(mod, type = "pearson"))
    testZeroInflation(mod, plot = T)

    col = as.data.frame(check_collinearity(mod)) %>%
    mutate(term = gsub("_", " ", Term)) %>%
    mutate(term = gsub("scaled", "", term)) %>%
    mutate(term = gsub(":", " x \n", term))

    plot.new()
    plot.window(xlim = c(1, nrow(col)), ylim = c(floor(min(col$VIF_CI_low)), ceiling(max(col$VIF_CI_high))))

    segments(1:nrow(col), col$VIF_CI_low, 1:nrow(col), col$VIF_CI_high)
    points(1:nrow(col), col$VIF, cex = 1.5, col = "black", pch = 19)

    axis(1, at = 1:nrow(col), labels = FALSE)
    axis(2, las = 2)

    title("Collinearity Check", adj=0)
    title(xlab = "Model Terms", ylab = "Variance Inflation Factor (VIF)")

    mtext(paste0("Diagnostics for full model for ", i, " accounting for spatial autocorrelation"), side = 3, line = 3, outer = TRUE, font = 2)

    ## Spatial autocorrelation using correlogram
    dat = filter(prop_forest_df, taxa == i)
    dat = mutate(dat, res = residuals(mod))
    correlog = with(dat, correlog(x, y, res, increment = 93694.99, resamp = 0, na.rm = TRUE))
    cor = as.data.frame(correlog[1:3]) %>% filter(n > 100 )
    plot(correlation ~ mean.of.class, data = cor, type = "o")

}

dev.off()

## Let's also run a model with all taxa (taxa as random intercept)

taxaLatModel = glmmTMB(prop_forest ~ scaled_dist_equator_1000km + (1|taxa) + (1|ECO_NAME),
    zi = ~scaled_dist_equator_1000km,
    weights = n_spec, data = prop_forest_df, family = "binomial")

save(taxaLatModel, file = paste0(gpath, "Results/taxaLatModel.RData"))

pdf(paste0(gpath, "Paper/Figures/Supplementary/taxaLatModel.pdf"), width = 10, height = 10)

    par(mfrow=c(2,3), oma = c(0, 0, 8, 0))
    testUniformity(taxaLatModel, plot = T)
    testDispersion(taxaLatModel, plot = T)
    plot(fitted(taxaLatModel), residuals(taxaLatModel, type = "pearson"))
    testZeroInflation(taxaLatModel, plot = T)

    col = as.data.frame(check_collinearity(taxaLatModel)) %>%
    mutate(term = gsub("_", " ", Term)) %>%
    mutate(term = gsub("scaled", "", term)) %>%
    mutate(term = gsub(":", " x \n", term))

    plot.new()
    plot.window(xlim = c(1, nrow(col)), ylim = c(floor(min(col$VIF_CI_low)), ceiling(max(col$VIF_CI_high))))

    segments(1:nrow(col), col$VIF_CI_low, 1:nrow(col), col$VIF_CI_high)
    points(1:nrow(col), col$VIF, cex = 1.5, col = "black", pch = 19)

    axis(1, at = 1:nrow(col), labels = FALSE)
    axis(2, las = 2)

    title("Collinearity Check", adj=0)
    title(xlab = "Model Terms", ylab = "Variance Inflation Factor (VIF)")

    mtext(paste0("Diagnostics for lat model for all taxa accounting for spatial autocorrelation"), side = 3, line = 3, outer = TRUE, font = 2)

    ## Spatial autocorrelation using correlogram
    prop_forest_df = mutate(prop_forest_df, res = residuals(taxaLatModel))
    correlog = with(prop_forest_df, correlog(x, y, res, increment = 93694.99, resamp = 0, na.rm = TRUE))
    cor = as.data.frame(correlog[1:3]) %>% filter(n > 100 )
    plot(correlation ~ mean.of.class, data = cor, type = "o")

dev.off()

taxaFullModel = glmmTMB(prop_forest ~ scaled_prop_forest_area*scaled_prop_land_area_deforested + scaled_disturbances + scaled_dist_equator_1000km + scaled_geological_forest_time + 
        scaled_geological_forest_stability + scaled_alpha_plant_diversity + (1|ECO_NAME) + (1|taxa),
    zi = ~scaled_prop_forest_area + scaled_prop_land_area_deforested,
    weights = n_spec, data = dat, family = binomial())

save(taxaFullModel, file = paste0(gpath, "Results/taxaFullModel.RData"))

pdf(paste0(gpath, "Paper/Figures/Supplementary/taxaFullModel.pdf"), width = 10, height = 10)

    par(mfrow=c(2,3), oma = c(0, 0, 8, 0))
    testUniformity(taxaFullModel, plot = T)
    testDispersion(taxaFullModel, plot = T)
    plot(fitted(taxaFullModel), residuals(taxaFullModel, type = "pearson"))
    testZeroInflation(taxaFullModel, plot = T)

    col = as.data.frame(check_collinearity(taxaFullModel)) %>%
    mutate(term = gsub("_", " ", Term)) %>%
    mutate(term = gsub("scaled", "", term)) %>%
    mutate(term = gsub(":", " x \n", term))

    plot.new()
    plot.window(xlim = c(1, nrow(col)), ylim = c(floor(min(col$VIF_CI_low)), ceiling(max(col$VIF_CI_high))))

    segments(1:nrow(col), col$VIF_CI_low, 1:nrow(col), col$VIF_CI_high)
    points(1:nrow(col), col$VIF, cex = 1.5, col = "black", pch = 19)

    axis(1, at = 1:nrow(col), labels = FALSE)
    axis(2, las = 2)

    title("Collinearity Check", adj=0)
    title(xlab = "Model Terms", ylab = "Variance Inflation Factor (VIF)")

    mtext(paste0("Diagnostics for full model for all taxa accounting for spatial autocorrelation"), side = 3, line = 3, outer = TRUE, font = 2)

    ## Spatial autocorrelation using correlogram
    prop_forest_df = mutate(prop_forest_df, res = residuals(taxaFullModel))
    correlog = with(prop_forest_df, correlog(x, y, res, increment = 93694.99, resamp = 0, na.rm = TRUE))
    cor = as.data.frame(correlog[1:3]) %>% filter(n > 100 )
    plot(correlation ~ mean.of.class, data = cor, type = "o")

dev.off()