#############################################
## Test whether HWI is significantly different
## between forest and non-forest species
## based on our categorisation
#############################################

library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(DHARMa)

gpath = "/home/ben/Documents/PhD/PropForestSpecies/"
setwd(gpath)

## Load in bird habitat data
bird_habitats = read.csv(paste0(gpath, "Data/Ranges/bird_habitats.csv")) %>% mutate(species = tolower(gsub(" ", "_", species)))

## Categorise as forest or non-forest species
birds = bird_habitats %>% group_by(species) %>%
  filter(suitability == "Suitable") %>% mutate(species = tolower(gsub(" ", "_", species))) %>%
  summarise(forest_habitats = length(grep("forest", habitat, ignore.case = T) %>% 
                                       grep("degraded", ., ignore.case = T, invert = T)), 
            nonforest_habitats = length(grep("forest", habitat, ignore.case = T, invert = T)), 
            total_habitats = n(), taxa = taxa[[1]]) %>% mutate(forest = ifelse(forest_habitats == total_habitats, 1, 0))

## Identify marine bird species
marine_bird_habitats = bird_habitats %>% group_by(species) %>% 
  summarise(marine_habitats = length(grep("marine", habitat, ignore.case = T))) %>% 
  mutate(marine = ifelse(marine_habitats > 0, "true", "false")) %>% dplyr::select(-marine_habitats)

## Remove marine bird species
birds = birds %>% left_join(marine_bird_habitats, by = c("species" = "species")) %>% filter(marine != "true") %>%
    dplyr::select(species, forest)

## Load in avonet for HWI data
avonet = read.csv(paste0(gpath, "Data/avonet_traits.csv")) %>% 
    dplyr::select(Species1, Hand.Wing.Index) %>%
    rename("species" = 1, "HWI" = 2) %>%
    mutate(species = tolower(gsub(" ", "_", species)))

## Join bird habitat data with HWI data
birds = birds %>% left_join(avonet)

## Make "forest" a factor for modelling
birds = birds %>% mutate(forest = as.factor(forest))

## Model if HWI of forest(1) is larger than non-forest(0)
HWImod = glmmTMB(HWI ~ forest, data = birds)
testSimulatedResiduals(HWImod) ## Looks a bit dodgy
testQuantiles(HWImod) ## Pretty close so not too bad

summary(HWImod)
confint(HWImod)

## Get data for plotting
modDat = as.data.frame(ggpredict(HWImod)) %>%
    rename("forest" = 1, "est" = 2, "std" = 3, "CIlow" = 4,
    "CIhigh" = 5, "group" = 6) %>%
    dplyr::select(forest, est, CIlow, CIhigh)

HWIfig = ggplot(modDat, aes(forest, est)) +
    geom_pointrange(aes(ymin = CIlow, ymax = CIhigh), size = 2.5, linewidth = 2.5) +
    theme_classic() +
    labs(x = NULL, y = "Hand Wing Index") +
    theme(text = element_text(size = 30)) +
    scale_x_discrete(labels = c("0" = "Non-Forest", "1" = "Forest"))

ggsave(paste0(gpath, "Paper/Figures/HWIfig.png"), plot = HWIfig, width = 10, height = 10)
