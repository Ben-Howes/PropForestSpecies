####################################################
## Assess HWI changes between
## 1) Forest vs Matrix
## 2) Fragmentation sensitive vs non-sensitive
## 3) Latitude
####################################################

library(tidyverse)
library(lme4)
library(janitor)
library(sf)
library(tidybayes)
library(ggeffects)
library(viridis)

gpath = "/home/ben/Documents/PhD/matrix_Response/"
setwd(gpath)

## Load hand wing index data
HWI = read.csv(paste0(gpath, "Data/AVONET_BirdLife_Habitat.csv")) %>% clean_names() %>% 
  dplyr::select(species1, hand_wing_index) %>% rename(species = species1, hwi = hand_wing_index) %>%
  mutate(species = gsub(" ", "_", tolower(species)))

## Load Geom and Taxa data

taxa = read_csv(paste0(gpath, "Data/taxa_PID_Lookup.csv")) %>% clean_names()

geoms = load(paste0(gpath, "Data/study_Geometries.Rdata"))
study_Centroids = lapply(study_Plot_Points, function(x) 
  data.frame(PID = x$PID[[1]], Latitude = mean(x$Latitude), Longitude = mean(x$Longitude))) %>% bind_rows() %>% 
  clean_names() %>% mutate(small_abs_lat = abs(latitude)/100)

########################
### Matrix vs Forest
########################

load(paste0(gpath, "Data/categorised_Species.Rdata"))
cat = cat %>% mutate(matrix = ifelse(slope > 0, 1, 0)) %>% left_join(study_Centroids) %>% 
  left_join(taxa) %>% filter(taxa == "Birds") %>% left_join(HWI)

ggplot(cat, aes(x = as.factor(matrix), y = hwi, fill = as.factor(matrix))) + geom_boxplot() + theme_classic() +
  labs(x = NULL, y = "Hand Wing Index") +
  theme(text = element_text(size = 20), legend.position = "none") +
  scale_fill_viridis_d(begin = 0.5) + 
  scale_x_discrete(labels = c('Forest Species','Matrix Species'))

hwi_mod = glmer(matrix ~ hwi + (hwi|pid), data = cat, family = "binomial")
## Positive significant slope, only p = 0.03 though

########################
### Fragmentation sensitive vs non-sensitive
########################

load(paste0(gpath, "Data/fragCat.Rdata"))
fragCat = fragCat %>% left_join(cat) %>% left_join(HWI) %>% filter(taxa == "Birds")

ggplot(filter(fragCat, !is.na(matrix)& matrix == 1), aes(as.factor(sen), hwi, fill = as.factor(sen))) + 
  geom_boxplot() + theme_classic() +
  labs(x = NULL, y = "Hand Wing Index") +
  theme(text = element_text(size = 20), legend.position = "none") +
  scale_fill_viridis_d(begin = 0.5) +
  scale_x_discrete(labels = c('Not Sensitive to Matrix \nFragmentation','Sensitive to Matrix \nFragmentation'))

hwi_sen_mod = lmer(hwi ~ sen + (sen|pid), data = filter(fragCat, matrix == 1))

########################
### Latitude
########################

ggplot(fragCat, aes(small_abs_lat, hwi, col = as.factor(matrix))) + geom_point() + 
  theme_classic() + geom_smooth(method = "lm") + 
  scale_colour_viridis_d(begin = 0.5, labels = c("Forest", "Matrix")) +
  labs(x = "Absolute Latitude", y = "Hand Wing Index", col = "Species Type") + 
  theme(text = element_text(size = 20))

hwi_lat_mod = lmer(hwi ~ small_abs_lat*matrix + (small_abs_lat|pid), data = fragCat)

