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

## Set path and wd
gpath = "/home/ben/Documents/PhD/PropForestSpecies/"
setwd(gpath)

## Get paths to models
result_paths = list.files(paste0(gpath, "Results/1000latResultsForest/"), full.names = T)

## Read all the csvs in and bind them together
results = mclapply(result_paths, read_csv, show_col_types = FALSE, mc.cores = 8) %>% bind_rows()

## Change names of predictors for the plot
results = results %>% 
           mutate(predictor = case_when(predictor == '(Intercept)' ~ "Intercept",
                                        predictor == "dist_equator_1000km" ~ "Latitude\n(Distance from Equator)"),
                taxa = factor(str_to_title(taxa), levels = c("Amphibians", "Birds", "Mammals", "Reptiles")))

############################
## Effect size plot
############################

effect_size_plot = ggplot(filter(results, predictor != "Intercept"), aes(est, y = predictor, col = taxa)) + 
                        stat_pointinterval(aes(linetype = after_stat(xmin < 0 & xmax > 0)), 
                                position = ggstance::position_dodgev(height = 0.4), linewidth = 15, 
                                size = 25, .width = 0.95) +
                        geom_vline(xintercept = 0, linetype = "dotted", linewidth = 1) +
                        theme_classic() +
                        labs(x = "Standardised Effect Size", y = NULL, col = NULL) +
                        scale_color_viridis_d(limits = c("Amphibians", "Birds", "Mammals", "Reptiles")) +
                        theme(text = element_text(size = 35)) +
                        guides(linetype = "none", col = guide_legend(byrow = TRUE, nrow = 1, title.vjust = -0.5)) +
                        theme(legend.position = "bottom",legend.box="horizontal",
                                legend.spacing.y = unit(0, 'cm'),
                                legend.box.background = element_rect(colour = "black", 
                                linewidth = 2),
                                legend.margin = margin(t=0,r=20,b=-5,l=15)) +
                        scale_y_discrete(limits = rev(c("Latitude\n(Distance from Equator)")))

effect_size_plot

ggsave(paste0(gpath, "Paper/Figures/lat_effect_size_plot.png"), plot = effect_size_plot, width = 20, height = 5)
save(effect_size_plot, file = paste0(gpath, "Results/lat_effect_size_plot_rfile.Rdata"))

#############################################
## Make predictions for examples in paper
#############################################

## Split data into each model for lapplying over
results_split = results %>% filter(taxa == "Birds") %>% group_split(run)

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

get_prediction = function(x) {

        dist_equator_1000km = x[2,2]
        intercept = x[1,2]

        make_predictions = function(x, y) {
                predicted =  (dist_equator_1000km*x) +
                        intercept
                return(predicted)
                }

        out = mapply(x = seq(min(dat$scaled_dist_equator_1000km ),
                max(dat$scaled_dist_equator_1000km), length.out = 100),
                make_predictions,
                SIMPLIFY = F) %>% bind_rows() %>%
                mutate(scaled_dist_equator_1000km = seq(min(dat$scaled_dist_equator_1000km ),
                max(dat$scaled_dist_equator_1000km), length.out = 100)) %>%
                mutate(taxa = x$taxa[[1]]) %>%
                mutate(run = x$run[[1]]) %>%
                mutate(est = exp(est)/1+exp(est))

}

predicted_results = mclapply(results_split, get_prediction, mc.cores = 8) %>% bind_rows()

## Unscale latitude
predicted_results = predicted_results %>% mutate(dist_equator_1000km = unscale(scaled_dist_equator_1000km, "scaled_dist_equator_1000km"))

## Find the mean at each dist equator to get a mean prediction
mean_predicted_results = predicted_results %>% 
        group_by(dist_equator_1000km) %>%
        summarise(propForest = mean(est))

## Get % forest species at different latitudes for use in text
mean_predicted_results %>% slice_min(dist_equator_1000km)
mean_predicted_results %>% filter(dist_equator_1000km > 1.9 & dist_equator_1000km < 2.1)
mean_predicted_results %>% filter(dist_equator_1000km > 4.9 & dist_equator_1000km < 5.1)
