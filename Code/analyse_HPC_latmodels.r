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
result_paths = list.files(paste0(gpath, "Results/1000latResults/"), full.names = T)

## Read all the csvs in and bind them together
results = mclapply(result_paths, read_csv, show_col_types = FALSE, mc.cores = 8) %>% bind_rows()

## Change names of predictors for the plot
results = results %>% 
           mutate(predictor = case_when(predictor == '(Intercept)' ~ "Intercept",
                                        predictor == "dist_equator_1000km" ~ "Distance from Equator\n(1000km)"),
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
                        guides(linetype = "none", col = guide_legend(byrow = TRUE, title.vjust = 5)) +
                        theme(legend.position = c(0.15, 0.25),
                                legend.spacing.y = unit(-0.75, 'cm'),
                                legend.box.background = element_rect(colour = "black", 
                                linewidth = 2),
                                legend.margin = margin(t=50,r=20,b=1,l=10)) +
                        scale_y_discrete(limits = rev(c("Distance from Equator\n(1000km)")))

effect_size_plot

ggsave(paste0(gpath, "Paper/Figures/lat_effect_size_plot.png"), plot = effect_size_plot, width = 20, height = 10)
