####################################################
## See if proportion of matrix species
## changes with latitude
####################################################

library(tidtidyverse)
library(lme4)
library(brms)
library(janitor)
library(tidybayes)
library(sf)
library(ggeffects)

gpath = "/home/ben/Documents/PhD/matrix_Response/"
setwd(gpath)

spec = read_csv(paste0(gpath, "Data/species_Occurrence.csv")) %>% clean_names() %>% 
  filter(buffer == 200) %>% mutate(matrix_cover = matrix_cover/100)
matrix_Studies = spec %>% group_by(pid) %>% slice_max(matrix_cover, with_ties = F) %>% 
  filter(matrix_cover > 0.9) %>% distinct(pid)
duplicate_Studies = c("PID0100", "PID0101", "PID0102", "PID0103", "PID0104", "PID0105",
                      "PID0106", "PID0107", "PID0108", "PID0109", "PID0110", "PID0111", "PID0113", "PID0114", "PID0117", "PID0119", "PID0124",
                      "PID0138", "PID1002", "PID1005")

spec = spec %>% filter(pid %in% matrix_Studies$pid & !pid %in% duplicate_Studies)

taxa = read_csv(paste0(gpath, "Data/taxa_PID_Lookup.csv")) %>% clean_names()

geoms = load(paste0(gpath, "Data/study_Geometries.Rdata"))
study_Centroids = lapply(study_Plot_Points, function(x) 
  data.frame(PID = x$PID[[1]], Latitude = mean(x$Latitude), Longitude = mean(x$Longitude))) %>% bind_rows() %>% 
  clean_names() %>% mutate(small_abs_lat = abs(latitude)/100)

spec = spec %>% left_join(taxa) %>% left_join(study_Centroids)

#################################################################################
## Run mixed effect models to test for change in effect of matrix cover
## with latitude
#################################################################################

load(paste0(gpath, "Data/categorised_Species.Rdata"))
cat = cat %>% mutate(matrix = ifelse(slope > 0, 1, 0)) %>% left_join(study_Centroids) %>% 
  left_join(taxa)
spec = spec %>% left_join(study_Centroids)

## Frequentist approach with different random effects
fMod = glmer(occurrence ~ matrix_cover + matrix_cover:small_abs_lat + (matrix_cover|species) + (1|pid), data = spec, family = "binomial")
fMod1 = glmer(occurrence ~ matrix_cover + matrix_cover:small_abs_lat + (matrix_cover|pid/species), data = spec, family = "binomial")
fMod2 = glmer(occurrence ~ matrix_cover + matrix_cover:small_abs_lat + (matrix_cover|taxa/pid/species), data = spec, family = "binomial")

predicted_Incidence_Plot = plot(ggemmeans(fMod2, terms = c("matrix_cover", "small_abs_lat[0, 0.23, 0.46]"))) + 
  labs(x = "Proportion Matrix Cover 200m", y = "Probability of Incidence", col = "Latitude", title = NULL) + 
  theme_classic() + theme(text = element_text(size = 20))

# What can be inferred from this?
## The effect of matrix on incidence changes with latitude, becoming more positive at higher latitudes
## I think this could mean more matrix species, but it could also be stronger positive effects on the same number of matrix species?

## Compare species with positive sloeps in above random models
## to those from individual bayesian models

fixed = fixef(fMod1)
random = coef(fMod1)[[1]] %>% mutate(pid = str_split(rownames(.), ":", simplify = T)[,2], 
                                     species = str_split(rownames(.), ":", simplify = T)[,1]) %>% 
  left_join(study_Centroids) %>% mutate(small_abs_lat = abs(latitude)/100) %>% 
  mutate(slope = matrix_cover + (.[,3]*small_abs_lat)) %>% 
  mutate(matrix = ifelse(slope > 0, 1, 0))

## Compare to individual models
ind_Rand_Comp = random %>% left_join(select(cat, species, slope, pid), by = c("species", "pid"), suffix = c("_rand", "_ind")) %>%
  mutate(match = ifelse(sign(slope_rand) == sign(slope_ind), 1, 0))

summary(as.factor(ind_Rand_Comp$match))
# 4788 match of 6098 = 4788/6098 = ~78%

##########################################################################
## Models looking at how change in proportion of matrix species varies 
## with latitude
## This is a definitive increase/decrease in matrix species
## But does require us to define what a matrix species is - propegate error?
################################################################################

prop_Matrix = cat %>% group_by(pid) %>% summarise(n_Matrix = sum(matrix), prop_Matrix = n_Matrix/n(), 
                                                  avg_Slope = median(slope), n_Species = n()) %>% 
                                                  left_join(study_Centroids) %>% left_join(taxa)

prop_Matrix_Mod = glmer(prop_Matrix ~ small_abs_lat + (small_abs_lat|taxa), weights = n_Species, 
                        data = prop_Matrix, family = "binomial") ## Sig pos

plot(ggpredict(prop_Matrix_Mod, terms = c("small_abs_lat"))) + 
  scale_colour_viridis_d() + labs(x = "Absolute Latitude", 
                                  y = "Proportion Matrix Species", title = NULL) +
  theme_classic() + theme(text = element_text(size = 20)) + 
  geom_point(data = prop_Matrix, aes(x = small_abs_lat, y = prop_Matrix), 
             size = 3, inherit.aes = F)

a = ggpredict(prop_Matrix_Mod, terms = c("small_abs_lat"))

plot(ggpredict(prop_Matrix_Mod, terms = c("small_abs_lat", "taxa"), type = "random"), ci = F) + 
  scale_colour_viridis_d() + labs(x = "Absolute Latitude", 
                                  y = "Proportion Matrix Species", col = "Taxa", title = NULL) +
  theme_classic() + theme(text = element_text(size = 20)) + 
  geom_point(data = prop_Matrix, aes(x = small_abs_lat, y = prop_Matrix, col = taxa), 
             size = 3, inherit.aes = F)

prop_Matrix_Bird_Mod = glm(prop_Matrix ~ small_abs_lat, weights = n_Species, data = filter(prop_Matrix, taxa == "Birds"), family = "binomial") ## Sig pos
prop_Matrix_Arthropod_Mod = glm(prop_Matrix ~ small_abs_lat, weights = n_Species, data = filter(prop_Matrix, taxa == "Arthropods"), family = "binomial") ## Sig pos
prop_Matrix_Herp_Mod = glm(prop_Matrix ~ small_abs_lat, weights = n_Species, data = filter(prop_Matrix, taxa == "Herps"), family = "binomial") ## Sig pos
prop_Matrix_Mammal_Mod = glm(prop_Matrix ~ small_abs_lat, weights = n_Species, data = filter(prop_Matrix, taxa == "Mammals"), family = "binomial") # Sig pos

####################################################################################################################################################
## Another way to analyse - based on Matt's draft of nature paper
## Predict probability species will be defined as matrix vs non-matrix (1 vs 0) based on abs latitude
## Same issue with defining what a matrix species is a propegating error
## use logistic reression and study as random effect
####################################################################################################################################################

altMod = glmer(matrix ~ small_abs_lat + (small_abs_lat|taxa/pid), data = cat, family = "binomial")
## Strong increase in probability of a species being classified as matrix in higher latitudes

plot(ggpredict(altMod, terms = c("small_abs_lat"))) + 
  scale_colour_viridis_d() + labs(x = "Absolute Latitude", y = "Proportion Matrix Species", title = NULL) +
  theme_classic() + theme(text = element_text(size = 20))

altBirdMod = glmer(matrix ~ small_abs_lat + (1|pid), data = filter(cat, taxa == "Birds"), family = "binomial") ## Sig pos
altArthropodMod = glmer(matrix ~ small_abs_lat + (1|pid), data = filter(cat, taxa == "Arthropods"), family = "binomial") ## Non-sig pos
altHerpMod = glmer(matrix ~ small_abs_lat + (1|pid), data = filter(cat, taxa == "Herps"), family = "binomial") # Just sig pos
altMammalMod = glmer(matrix ~ small_abs_lat + (1|pid), data = filter(cat, taxa == "Mammals"), family = "binomial") ## Sig Pos

## Save all models including individual models
save(fMod, fMod1, fMod2, prop_Matrix, prop_Matrix_Mod,prop_Matrix_Bird_Mod, prop_Matrix_Arthropod_Mod, 
     prop_Matrix_Herp_Mod, prop_Matrix_Mammal_Mod, altMod, altBirdMod, altArthropodMod, altHerpMod, altMammalMod,
     file = paste0(gpath, "Data/freq_Models_and_Results.Rdata"))

## Load
load(paste0(gpath, "Data/freq_Models_and_Results.Rdata"))

###################
## Bayesian
###################

bayesMod = brm(occurrence ~ matrix_cover + matrix_cover:small_abs_lat + (matrix_cover|taxa/pid/species), data = spec, family = "bernoulli",
               chains = 2, cores = 2, iter = 50)

bayesMammal = brm(occurrence ~ matrix_cover + matrix_cover:small_abs_lat + (matrix_cover|pid/species), data = filter(spec, taxa == "Mammals"), 
                  family = "bernoulli", chains = 2, cores = 2)
bayesBirds = brm(occurrence ~ matrix_cover + matrix_cover:small_abs_lat + (matrix_cover|pid/species), data = filter(spec, taxa == "Birds"), 
                 family = "bernoulli", chains = 2, cores = 2, iter = 50)
bayesHerps = brm(occurrence ~ matrix_cover + matrix_cover:small_abs_lat + (matrix_cover|pid/species), data = filter(spec, taxa == "Herps"), 
                 family = "bernoulli", chains = 2, cores = 2, iter = 200)
bayesArths = brm(occurrence ~ matrix_cover + matrix_cover:small_abs_lat + (matrix_cover|pid/species), data = filter(spec, taxa == "Arthropods"), 
                 family = "bernoulli", chains = 2, cores = 2, iter = 50)

# load(paste0(gpath, "Data/species_Bayes_Models.Rdata"))
# save(bayesMod, bayesMammal, bayesBirds, bayesHerps, bayesArths, file = paste0(gpath, "Data/bayes_Models.Rdata"))
