#############################################################
## Analyse the spatmodels run on the HPC
## using glmmTMB to account for spatialautocorrelation
## Analyse the output of 4000 models run on the HPC
## using glmmTMB models accounting for spatial autocorrelation
## since models were so large, estimates and CI were output
## on HPC, and can then be analysed here
## 1000 models were run per taxa group
## and all variables are scaled so are comparable
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
result_paths = list.files(paste0(gpath, "Results/1000spatResults/"), full.names = T)

## Read all the csvs in and bind them together
results = mclapply(result_paths, read_csv, show_col_types = FALSE, mc.cores = 8) %>% bind_rows()

## Change names of predictors for the plot
results = results %>% 
           mutate(predictor = case_when(predictor == '(Intercept)' ~ "Intercept",
                                        predictor == "prop_forest_area" ~ "Proportion Forested Area",
                                        predictor == "historical_forest_loss" ~ "Proportion Historic Forest Loss",
                                        predictor == "disturbances" ~ "Naturally Disturbed Area",
                                        predictor == "dist_equator_1000km" ~ "Distance from Equator\n(1000km)",
                                        predictor == "geological_forest" ~ "Long Term Forest Status",
                                        predictor == "prop_forest_area:historical_forest_loss" ~ "Proportion Forested Area x\nProportion Historic Forest Loss"),
                taxa = factor(str_to_title(taxa), levels = c("Amphibians", "Birds", "Mammals", "Reptiles")))

############################
## Effect size plot
############################

effect_size_plot = ggplot(filter(results, predictor != "Intercept"), aes(est, y = predictor, col = taxa)) + 
                        stat_pointinterval(aes(linetype = after_stat(xmin < 0 & xmax > 0)), 
                                position = ggstance::position_dodgev(height = -0.5), linewidth = 15, 
                                size = 25, .width = 0.95) +
                        geom_vline(xintercept = 0, linetype = "dotted", linewidth = 1) +
                        theme_classic() +
                        labs(x = "Standardised Effect Size", y = NULL, col = "Taxa") +
                        scale_color_viridis_d(limits = c("Amphibians", "Birds", "Mammals", "Reptiles")) +
                        theme(text = element_text(size = 35)) +
                        guides(linetype = FALSE) +
                        theme(legend.position = c(0.15, 0.725),
                                legend.box.background = element_rect(colour = "black", 
                                linewidth = 1.5)) +
                        scale_y_discrete(limits = rev(c("Proportion Forested Area", "Proportion Historic Forest Loss", 
                                                        "Proportion Forested Area x\nProportion Historic Forest Loss",
                                                        "Naturally Disturbed Area", "Long Term Forest Status",
                                                        "Distance from Equator\n(1000km)")))

effect_size_plot

ggsave(paste0(gpath, "Paper/Figures/effect_size_plot.png"), plot = effect_size_plot, width = 20, height = 10)

############################
## Plot of interaction term between forest cover and historic forest loss
############################

## To do this we will use our bootstrapped samples and calculate a line for each run of our model
## The form of the interaction term equation is predicted = pred1x + pred2y + interactionxy + c
## So I need to do that for each of my models

## Split data into each model for lapplying over
results_split = results %>% filter(taxa == "Birds") %>% group_split(run)

## Load in our actual data
dat = read_csv(paste0(gpath, "Data/proportion_forest_species_analysis_data.csv"))

## Just keep values we are interested in and
## Scale variables so they are comparable to our model
dat = dat %>% dplyr::select(prop_forest_area, historical_forest_loss, taxa) %>%
                mutate(scaled_prop_forest_area = scale(prop_forest_area),
                scaled_historical_forest_loss = scale(historical_forest_loss))

## I think doing values of -1 SD, 0 SD, and +1SD make sense
## Though there is no -1 SD, so we will just use the lowest which is -0.571
## Which will be No, Average, and high

# Then we can get the maximum proporion_forest_area that accompany these
dat %>% slice_min(scaled_historical_forest_loss) %>% slice_max(scaled_prop_forest_area, with_ties = F) ## 2.68
dat %>% slice_min(scaled_historical_forest_loss) %>% slice_min(scaled_prop_forest_area, with_ties = F) ## -0.571

dat %>% filter(scaled_historical_forest_loss > 0) %>% slice_max(scaled_prop_forest_area, with_ties = F) ## 2.52
dat %>% filter(scaled_historical_forest_loss > 0) %>% slice_min(scaled_prop_forest_area, with_ties = F) ## -0.571

dat %>% filter(scaled_historical_forest_loss > 1) %>% slice_max(scaled_prop_forest_area, with_ties = F) ## 1.18
dat %>% filter(scaled_historical_forest_loss > 1) %>% slice_min(scaled_prop_forest_area, with_ties = F) ## -0.571

## So these can be our limits on this prediction grid
prediction_grid = expand.grid(x = seq(-0.571, 2.68, length.out = 100), y = -0.571) %>%
                  bind_rows(expand.grid(x = seq(-0.571, 2.52, length.out = 100), y = 0)) %>%
                  bind_rows(expand.grid(x = seq(-0.571, 1.18, length.out = 100), y = 1))

get_prediction = function(x) {

        prop_forest_area = x[2,2]
        historic_forest_loss = x[3,2]
        interaction = x[7,2]
        intercept = x[1,2]

make_predictions = function(x, y) {
        predicted =  (prop_forest_area*x) + 
                     (historic_forest_loss*y) + 
                     (x*y*interaction) + 
                     intercept
        return(predicted)
        }

out = mapply(x = prediction_grid$x, 
             y = prediction_grid$y, make_predictions,
             SIMPLIFY = F) %>% bind_rows() %>%
             mutate(prop_forest_area = prediction_grid$x) %>%
             mutate(historic_forest_loss = prediction_grid$y) %>%
             mutate(taxa = x$taxa[[1]]) %>%
             mutate(run = x$run[[1]]) %>%
             mutate(est = exp(est)/1+exp(est))

}

predicted_results = mclapply(results_split, get_prediction, mc.cores = 8) %>% bind_rows()

facet_labels =  c(
  "-0.571"="Low Amount of\nHistoric Forest Loss",
  "0"="Average Amount of\nHistoric Forest Loss",
  "1"="High Amount of\nHistoric Forest Loss"
)

interaction_term_plot = ggplot(predicted_results, aes(prop_forest_area, est)) + 
                        stat_lineribbon(alpha = 0.2, fill = "grey40") + 
                        facet_rep_wrap(. ~ historic_forest_loss, labeller = as_labeller(facet_labels)) +
                        theme_classic() +
                        labs(x = "Proportion Forested Area (Standardised)", y = 
                                "Probability of Forest Species") +
                        theme(text = element_text(size = 35))

interaction_term_plot

ggsave(paste0(gpath, "Paper/Figures/interaction_term_plot.png"), plot = interaction_term_plot, width = 20, height = 10)
