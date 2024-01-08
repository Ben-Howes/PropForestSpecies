#######################################################
## Run models for proportion of forest species
## usiing quasi binomial models to account for
## overdispersion
## then plot comparisons with normal binomial models
## to show overdispersion does not effect results
#######################################################

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
                                            scaled_altitude = scale(altitude)[,1],
                                            scaled_n_spec = scale(n_spec)[,1],
                                            ECO_NAME = as.factor(ECO_NAME), 
                                            taxa = as.factor(taxa)) 

##########################################################################################
## Run latitude models accounting for spatial autocorrelation
##########################################################################################

quasiLatModels = list()

for(i in unique(prop_forest_df$taxa)) {
    print(i)
    dat = filter(prop_forest_df, taxa == i) %>% mutate(id = as.factor(1:nrow(.)))
    mod = tryCatch(bam(prop_forest ~ scaled_dist_equator_1000km + s(x, y, bs = "gp") + s(ECO_NAME, bs = "re"), 
    weights = n_spec, data = dat, family = quasibinomial(), discrete = TRUE, nthreads = 6), error = function(e) NULL,
    warning = function(w) NULL)
    quasiLatModels[[i]] = mod
}

save(quasiLatModels, file = paste0(gpath, "Results/quasiBinomialModels/quasiLatModels.RData"))

##########################################################################################
## Run full models accounting for spatial autocorrelation
##########################################################################################

quasiFullModels = list()

for(i in unique(prop_forest_df$taxa)) {
    dat = filter(prop_forest_df, taxa == i)
    mod = tryCatch(bam(prop_forest ~ scaled_prop_forest_area + scaled_prop_land_area_deforested + scaled_disturbances + scaled_dist_equator_1000km + scaled_geological_forest_time + 
    scaled_geological_forest_stability + scaled_alpha_plant_diversity + scaled_altitude +scaled_n_spec + s(x, y, bs = "gp") + s(ECO_NAME, bs = "re"),
    weights = n_spec, data = dat, family = quasibinomial(), discrete = TRUE, nthreads = 6), error = function(e) NULL,
    warning = function(w) NULL)
    quasiFullModels[[i]] = mod
}

save(quasiFullModels, file = paste0(gpath, "Results/quasiBinomialModels/quasiFullModels.RData"))

## Let's also run a model with all taxa (taxa as random intercept)

quasiTaxaLatModel = bam(prop_forest ~ scaled_dist_equator_1000km + s(x, y, bs = "gp") + 
    s(ECO_NAME, bs = "re") + s(taxa, bs = "re"),
    weights = n_spec, data = prop_forest_df, family = quasibinomial(), discrete = TRUE, nthreads = 6)

save(quasiTaxaLatModel, file = paste0(gpath, "Results/quasiBinomialModels/quasiTaxaLatModel.RData"))

quasiTaxaFullModel = bam(prop_forest ~ scaled_prop_forest_area + scaled_prop_land_area_deforested + scaled_disturbances + 
    scaled_dist_equator_1000km + scaled_geological_forest_time + 
    scaled_geological_forest_stability + scaled_alpha_plant_diversity + 
    scaled_altitude + scaled_n_spec +
    s(x, y, bs = "gp") + s(ECO_NAME, bs = "re") + s(taxa, bs = "re") + s(taxa, bs = "re"),
    weights = n_spec, data = prop_forest_df, family = quasibinomial(), discrete = TRUE, nthreads = 6)

save(quasiTaxaFullModel, file = paste0(gpath, "Results/quasiBinomialModels/quasiTaxaFullModel.RData"))

## Load in normal binomial models
load("Results/randLatModels.RData")
load("Results/randFullModels.RData")
load("Results/taxaLatModel.RData")
load("Results/taxaFullModel.RData")

## Convert all model outputs to tidy tables

## Combine the individual taxa models with the combined taxa models for plotting
randLatModels = c(randLatModels, list(all = taxaLatModel))
randFullModels = c(randFullModels, list(all = taxaFullModel))

quasiLatModels = c(quasiLatModels, list(all = quasiTaxaLatModel))
quasiFullModels = c(quasiFullModels, list(all = quasiTaxaFullModel))

latResults = lapply(names(randLatModels), function(x)
        tidy(randLatModels[[x]], parametric = TRUE, conf.int = TRUE) %>%
        mutate(taxa = x)) %>%
        bind_rows() %>%
        mutate(term = case_when(term == '(Intercept)' ~ "Intercept",
        term == "scaled_dist_equator_1000km" ~ "Latitude (Distance\nfrom Equator)"),
        taxa = factor(str_to_title(taxa), levels = c("All", "Amphibians", "Birds", "Mammals", "Reptiles"))) %>%
        mutate(model = "Binomial")

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
        taxa = factor(str_to_title(taxa), levels = c("All", "Amphibians", "Birds", "Mammals", "Reptiles"))) %>%
        mutate(model = "Binomial")

## Add column to alpha non-significant results
fullResults = fullResults %>%
        mutate(alpha = ifelse(sign(conf.low) == sign(conf.high), "sig", "nonsig"))

## Do the same for quasi models

quasiLatResults = lapply(names(quasiLatModels), function(x)
        tidy(quasiLatModels[[x]], parametric = TRUE, conf.int = TRUE) %>%
        mutate(taxa = x)) %>%
        bind_rows() %>%
        mutate(term = case_when(term == '(Intercept)' ~ "Intercept",
        term == "scaled_dist_equator_1000km" ~ "Latitude (Distance\nfrom Equator)"),
        taxa = factor(str_to_title(taxa), levels = c("All", "Amphibians", "Birds", "Mammals", "Reptiles"))) %>%
        mutate(model = "Quasi-Binomial")

quasiFullResults = lapply(names(quasiFullModels), function(x)
        tidy(quasiFullModels[[x]], parametric = TRUE, conf.int = TRUE) %>%
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
        taxa = factor(str_to_title(taxa), levels = c("All", "Amphibians", "Birds", "Mammals", "Reptiles"))) %>%
        mutate(model = "Quasi-Binomial")

## Add column to alpha non-significant results
quasiFullResults = quasiFullResults %>%
        mutate(alpha = ifelse(sign(conf.low) == sign(conf.high), "sig", "nonsig"))

## Bind rows of results together
latResults = bind_rows(latResults, quasiLatResults)
fullResults = bind_rows(fullResults, quasiFullResults)

latEffectPlot = ggplot(filter(latResults, term != "Intercept"), aes(estimate, term, col = taxa, shape = model)) +
        geom_pointrange(aes(xmin = conf.low, xmax = conf.high), position = ggstance::position_dodgev(height = -1), size = 2, linewidth = 2) +
        geom_vline(xintercept = 0, linetype = "dotted", linewidth = 1) +
        theme_classic() +
        labs(x = "Standardised Effect Size", y = NULL, col = "Taxa", shape = "Model Type") +
        scale_color_viridis_d(limits = c("All", "Amphibians", "Birds", "Mammals", "Reptiles")) +
        theme(text = element_text(size = 35)) +
        guides(linetype = "none", shape = "none",
        col = guide_legend(byrow = TRUE, nrow = 5)) +
        theme(legend.position = c(0.85, 0.53),legend.box = "horizontal",
        legend.spacing.y = unit(0.15, 'cm'),
        legend.box.background = element_rect(colour = "black", 
        linewidth = 2),
        legend.margin = margin(t=10,r=20,b=10,l=15)) +
        xlim(-1.25, 0.5)

latEffectPlot

ggsave(paste0(gpath, "Paper/Figures/Supplementary/quasiLatEffectPlot.png"), plot = latEffectPlot, width = 20, height = 4)

effectPlot = ggplot(filter(fullResults, term != "Intercept"), aes(estimate, term, col = taxa, alpha = alpha, shape = model)) +
        geom_pointrange(aes(xmin = conf.low, xmax = conf.high), position = ggstance::position_dodgev(height = -1), size = 2, linewidth = 2) +
        geom_vline(xintercept = 0, linetype = "dotted", linewidth = 1) +
        theme_classic() +
        labs(x = "Standardised Effect Size", y = NULL, col = "Taxa", shape = "Model Type") +
        scale_color_viridis_d(limits = c("All", "Amphibians", "Birds", "Mammals", "Reptiles")) +
        theme(text = element_text(size = 35)) +
        guides(linetype = "none", shape = "none", col = guide_legend(byrow = TRUE, nrow = 5), 
        alpha = "none") +
        theme(legend.position = c(0.13, 0.85),legend.box = "horizontal",
        legend.spacing.y = unit(0.15, 'cm'),
        legend.box.background = element_rect(colour = "black", 
        linewidth = 2),
        legend.margin = margin(t=10,r=20,b=10,l=15)) +
        xlim(-0.85, 0.5) +
        scale_alpha_manual(values = c(0.5, 1)) +
        scale_y_discrete(limits = rev(c("Total Species Richness", "Current Forest Cover", "Historical Deforestation",
                                "Geological Forest Stability", "Altitude", "Naturally Disturbed Area",
                                "Plant Alpha Diversity", "Maximum Geological\nForest Time",
                                "Latitude (Distance\nfrom Equator)"))) 

effectPlot

ggsave(paste0(gpath, "Paper/Figures/Supplementary/quasiREffectPlot.png"), plot = effectPlot, width = 20, height = 12)
