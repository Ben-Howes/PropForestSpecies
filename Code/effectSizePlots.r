#############################################################
## Analyse the latmodels run on the HPC
## using glmmTMB to account for spatialautocorrelation
## Analyse the output of 4000 models run on the HPC
## using glmmTMB models accounting for spatial autocorrelation
## since models were so large, estimates and CI were output
## on HPC, and can then be analysed here
## 1000 models were run per taxa group
## only variable used is distance from equator, which is scaled
#############################################################

library(tidyverse)
library(ggdist) ## For halfeye density plots
library(parallel) ## mclapply etc
library(viridis) ## Colour scheme
library(ggstance) ## For specific positioning in ggplot
library(glmmTMB) ## extracting model output
library(ggeffects) ## Model predictions 
library(lemon) ## Facet_rep_wrap
library(broom)
library(broom.mixed) ## tidy
library(sf)

## Set path and wd
gpath = "/home/ben/Documents/PhD/PropForestSpecies/"
setwd(gpath)

## Load in models
load("Results/randLatModels.RData")
load("Results/randFullModels.RData")
load("Results/taxaLatModel.RData")
load("Results/taxaFullModel.RData")

## Combine the individual taxa models with the combined taxa models for plotting
randLatModels = c(randLatModels, list(all = taxaLatModel))
randFullModels = c(randFullModels, list(all = taxaFullModel))

latResults = lapply(names(randLatModels), function(x)
        tidy(randLatModels[[x]], parametric = TRUE, conf.int = TRUE) %>%
        mutate(taxa = x)) %>%
        bind_rows() %>%
        mutate(term = case_when(term == '(Intercept)' ~ "Intercept",
        term == "scaled_dist_equator_1000km" ~ "Latitude (Distance\nfrom Equator)"),
        taxa = factor(str_to_title(taxa), levels = c("All", "Amphibians", "Birds", "Mammals", "Reptiles")))

fullResults = lapply(names(randFullModels), function(x)
        tidy(randFullModels[[x]], parametric = TRUE, conf.int = TRUE) %>%
        mutate(taxa = x)) %>%
        bind_rows() %>%
        mutate(term = case_when(term == '(Intercept)' ~ "Intercept",
                        term == "scaled_prop_forest_area" ~ "Current Forest Cover",
                        term == "scaled_prop_land_area_deforested" ~ "Historical Deforestation",
                        term == "scaled_disturbances" ~ "Naturally Disturbed Area",
                        term == "scaled_dist_equator_1000km" ~ "Latitude (Distance\nfrom Equator)",
                        term == "scaled_geological_forest_time" ~ "Maximum Geological\nForest Time",
                        term == "scaled_geological_forest_stability" ~ "Geological Forest Stability",
                        term == "scaled_alpha_plant_diversity" ~ "Plant Alpha Diversity",
                        term == "scaled_altitude" ~ "Altitude",
                        term == "scaled_n_spec" ~ "Total Species Richness"),
        taxa = factor(str_to_title(taxa), levels = c("All", "Amphibians", "Birds", "Mammals", "Reptiles")))

## Add column to alpha non-significant results
fullResults = fullResults %>%
        mutate(alpha = ifelse(sign(conf.low) == sign(conf.high), "sig", "nonsig"))

############################
## Effect size plots
############################

## LATITUDE ONLY EFFECT PLOTS

latEffectPlot = ggplot(filter(latResults, term != "Intercept"), aes(estimate, term, col = taxa)) +
        geom_pointrange(aes(xmin = conf.low, xmax = conf.high), position = ggstance::position_dodgev(height = -1), size = 2, linewidth = 2) +
        geom_vline(xintercept = 0, linetype = "dotted", linewidth = 1) +
        theme_classic() +
        labs(x = "Standardised Effect Size", y = NULL, col = "Taxa") +
        scale_color_viridis_d(limits = c("All", "Amphibians", "Birds", "Mammals", "Reptiles")) +
        theme(text = element_text(size = 35)) +
        guides(linetype = "none", col = guide_legend(byrow = TRUE, nrow = 5, title.vjust = 1)) +
        theme(legend.position = c(0.85, 0.52),legend.box = "horizontal",
        legend.spacing.y = unit(0.15, 'cm'),
        legend.box.background = element_rect(colour = "black", 
        linewidth = 2),
        legend.margin = margin(t=10,r=20,b=10,l=15)) +
        xlim(-1.25, 0.5)

latEffectPlot

ggsave(paste0(gpath, "Paper/Figures/latEffectPlot.png"), plot = latEffectPlot, width = 20, height = 4)
ggsave(paste0(gpath, "Paper/Figures/latEffectPlot.pdf"), plot = latEffectPlot, width = 20, height = 4)

## FULL RESULTS PLOT

effectPlot = ggplot(filter(fullResults, term != "Intercept"), aes(estimate, term, col = taxa, alpha = alpha)) +
        geom_pointrange(aes(xmin = conf.low, xmax = conf.high), position = ggstance::position_dodgev(height = -0.75), size = 2, linewidth = 2) +
        geom_vline(xintercept = 0, linetype = "dotted", linewidth = 1) +
        theme_classic() +
        labs(x = "Standardised Effect Size", y = NULL, col = "Taxa") +
        scale_color_viridis_d(limits = c("All", "Amphibians", "Birds", "Mammals", "Reptiles")) +
        theme(text = element_text(size = 35)) +
        guides(linetype = "none", col = guide_legend(byrow = TRUE, nrow = 5, title.vjust = 1), 
        alpha = "none") +
        theme(legend.position = c(0.13, 0.86),legend.box = "horizontal",
        legend.spacing.y = unit(0.15, 'cm'),
        legend.box.background = element_rect(colour = "black", 
        linewidth = 2),
        legend.margin = margin(t=10,r=20,b=10,l=15)) +
        xlim(-0.85, 1.5) +
        scale_alpha_manual(values = c(0.3, 1)) +
        scale_y_discrete(limits = rev(c("Total Species Richness", "Current Forest Cover", "Historical Deforestation",
                                "Geological Forest Stability", "Altitude", "Naturally Disturbed Area",
                                "Plant Alpha Diversity", "Maximum Geological\nForest Time",
                                "Latitude (Distance\nfrom Equator)"))) 

effectPlot

ggsave(paste0(gpath, "Paper/Figures/effectPlot.png"), plot = effectPlot, width = 20, height = 12)
ggsave(paste0(gpath, "Paper/Figures/effectPlot.pdf"), plot = effectPlot, width = 20, height = 12)

#############################################
## Make predictions for examples in paper
#############################################

## Split data into each model for lapplying over
results_split = latResults %>% filter(taxa == "All" & model == "Binomial")

## Load in our actual data
dat = read_csv(paste0(gpath, "Data/proportion_forest_species_analysis_data.csv"))

## Just keep values we are interested in and
## Scale variables so they are comparable to our model
dat = dat %>% dplyr::select(dist_equator_1000km, taxa) %>%
                mutate(scaled_dist_equator_1000km = scale(dist_equator_1000km))

## Function to unscale values
unscale = function(x, term) {
        val = x * attr(dat[[term]], 'scaled:scale') + attr(dat[[term]], 'scaled:center')
        return(val)
}

newdat = data.frame(scaled_dist_equator_1000km = seq(min(dat$scaled_dist_equator_1000km), 
        max(dat$scaled_dist_equator_1000km), length.out = 100))

newdat = newdat %>% mutate(predicted = predict(taxaLatModel, newdata = newdat, type = "response", exclude= c("s(taxa)", "s(ECO_NAME)", "s(x,y)"),
        newdata.guaranteed = TRUE, discrete = FALSE))

## Unscale distance from equator
newdat = newdat %>% mutate(dist_equator_1000km = unscale(scaled_dist_equator_1000km, "scaled_dist_equator_1000km"))

## Convert distances to latitude
## So distance from equator need to be multipled by 1000000 to get to metres
## and an arbitrary x value of 0 added to the data frame
newdat = newdat %>% mutate(y = dist_equator_1000km*1000000, x = 0)
sfPredictions = st_as_sf(newdat, coords = c("x", "y"))
behr = "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
st_crs(sfPredictions) = behr

## Transform from behrmann to epsg 4326
sfPredictions = st_transform(sfPredictions, 4326)

## Split xy coordinates to get lat and long
sfPredictions = sfPredictions %>% mutate(x = st_coordinates(.)[,1], y = st_coordinates(.)[,2]) %>% as.data.frame()

## Get % forest species at different latitudes for use in text
sfPredictions %>% slice_min(y)
sfPredictions %>% filter(y > 29 & y < 31)
sfPredictions %>% filter(y > 49 & y < 51)
