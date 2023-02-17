################################################################
## Categorise species based on relative habitat abundance
## Data taken from ESA 2020 land cover dataset
## and biofrag for species data
################################################################

library(tidyverse)
library(lme4)
library(brms)
library(janitor)
library(tidybayes)
library(sf)
library(ggeffects)

gpath = "/home/ben/Documents/PhD/matrix_Response/"
setwd(gpath)

landcover = read_csv(paste0(gpath, "Data/esa_landuse.csv")) %>% dplyr::select(-...1)

## Take the highest proportion landcover for each site and then categorise that site
## by that land cover. Since it is 200m most will be close to 100%
site_landcover = landcover %>% group_by(pid, plot) %>% slice_max(percentage, with_ties = F)

## How many non-tree sites do we have?
summary(as.factor(site_landcover$landuse)) ## 1332 of 11748 sites

## How many studies have at least one non-tree site?
non_tree = site_landcover %>% filter(landuse != "trees") %>% ungroup() %>% 
  distinct(pid) ## 43 studies of 114 total. Not great.

spec = read_csv(paste0(gpath, "Data/species_Occurrence.csv")) %>% clean_names() %>% 
  filter(buffer == 200) %>% dplyr::select(-buffer, -matrix_cover)

duplicate_Studies = c("PID0100", "PID0101", "PID0102", "PID0103", "PID0104", "PID0105",
                      "PID0106", "PID0107", "PID0108", "PID0109", "PID0110", "PID0111", "PID0113", "PID0114", "PID0117", "PID0119", "PID0124",
                      "PID0138", "PID1002", "PID1005")

spec = spec %>% left_join(site_landcover) %>% filter(pid %in% non_tree$pid & !pid %in% duplicate_Studies) %>%
  rename(cover = percentage) %>% dplyr::select(-area)

taxa = read_csv(paste0(gpath, "Data/taxa_PID_Lookup.csv")) %>% clean_names()

geoms = load(paste0(gpath, "Data/study_Geometries.Rdata"))
study_centroids = lapply(study_Plot_Points, function(x) 
  data.frame(PID = x$PID[[1]], Latitude = mean(x$Latitude), Longitude = mean(x$Longitude))) %>% bind_rows() %>% 
  clean_names() %>% mutate(small_abs_lat = abs(latitude)/100)

spec = spec %>% left_join(taxa) %>% left_join(study_centroids)

## Relative Habitat Usaeg taken from O-Reilly et al 2022

RHU = function(data, study) {
  
  data = filter(data, pid == study)
  
  habitats = data %>% distinct(landuse) %>% first()
  species = data %>% distinct(species) %>% first()
  sites = data %>% distinct(plot) %>% nrow()
  species_totals = data %>% group_by(species) %>% summarise(n = sum(occurrence))
  

  habitat_RHU = function(data, habitat) {
    
    data = data %>% filter(landuse == habitat)
    n_habitat = data %>% distinct(plot) %>% nrow()
    
    RHU_species = function(data, s) {
      
      data = filter(data, species == s)
      
      n = data %>% dplyr::select(occurrence) %>% sum() # number of individuals in this habitat
      p = n_habitat # number of sitse of habitat in study
      N = species_totals %>% filter(species == s) %>% select(n) %>% first() # total number of individuals across all habitats
      P = sites # total number of sites in the study
    
      ## RHU returns inf if the species only occurs in one habitat type
      ## So check this and if it only occurs in one habitat type set RHU to 5 as per O'Reilly et al 2022
      if(n != N) { 
        RHU = (n/p)/((N-n)/(P-p))
      } else {
        RHU = 5
      }
      
      out = data.frame(pid = study, species = s, habitat = habitat, RHU = RHU, habitat_occ = n, all_occ = N,
                       n_habitat = p, all_habitat = P)
      
    }
    
    out = lapply(species, RHU_species, data = data) %>% bind_rows()
    return(out)
  }
  
  res = lapply(habitats, habitat_RHU, data = data) %>% bind_rows()

return(res)
  
}

PIDs = spec %>% ungroup() %>% distinct(pid) %>% first()
all_RHU = lapply(PIDs, RHU, data = spec) %>% bind_rows()  

## RHU is sensitive to p and P so need to implement minimum occurrence for species
## Since we have lots of rare species we'll use total occurrence of >= 10
## and occurrence at particular habitat of >= 5

all_RHU = all_RHU %>% filter(all_occ >= 10 & habitat_occ >= 5)

all_RHU %>% arrange(pid, species)
all_RHU %>% distinct(species) %>% nrow() ## 1029 species left

## In general RHU < 1 is not associated, 1-2 is associated, and > 2 is strongly associated
## So I think I can categorise a species as a specialist if it has RHU > 1 and everything else < 1

specialists = all_RHU %>% group_by(pid, species) %>% slice_max(RHU, with_ties = F) %>% filter(RHU > 1)

prop_matrix = specialists %>% mutate(matrix = ifelse(habitat != "trees", 1, 0)) %>% 
  group_by(pid) %>% summarise(n_spec = n(), n_matrix = sum(matrix), prop_matrix = n_matrix/n_spec) %>%
  left_join(taxa) %>% left_join(study_centroids)

ggplot(prop_matrix, aes(small_abs_lat, prop_matrix)) + geom_point() + 
  theme_classic() + geom_smooth(method = "lm")

prop_matrix_mod = glmer(prop_matrix ~ small_abs_lat + (small_abs_lat|taxa), data = prop_matrix, 
                        weights = n_spec, family = "binomial")

plot(ggemmeans(prop_matrix_mod, terms = c("small_abs_lat"))) + theme_classic() +
  labs(x = "Absolute Latitude", y = "Proportion Matrix Species", title = NULL) + 
  theme(text = element_text(size = 20)) + 
  geom_point(data = prop_matrix, aes(x = small_abs_lat, y = prop_matrix), 
             size = 3, inherit.aes = F)
