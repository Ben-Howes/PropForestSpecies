################################################################
## Extract species matrix from biofrag dataset
## and combine with landscape metrics (forest cover %)

## Ben Howes, 2nd Year PhD Student, Imperial College London
################################################################

library(tidyverse)
library(sf)
library(raster)
library(geojsonio)
library(landscapemetrics)
library(janitor)
library(gtools)
library(pbapply)
library(units)

gpath = "/home/ben/Documents/PhD/matrix_Response/"
setwd(gpath)

## List all biofrag study folders

biofrag_Studies = list.files("Data/BIOFRAG_Jan2020Status", pattern = "PID*")

#########################
## Load in Species Data
#########################

## Make a list of all study PIDs and just keep percentage cover of matrix (0) and forest (1) for now from landscape metrics

study_PIDs = sapply(1:length(biofrag_Studies), function (x) {str_split(biofrag_Studies, pattern = "_")[[x]][1]})

## Load in Species Matrix

species_Matrix = mapply(function(x,y) {
    read_csv(paste0(gpath, "Data/BIOFRAG_Jan2020Status/", x, "/", y, "_species_matrix.csv"), locale = locale(encoding = "windows-1252")) %>% clean_names()
}, biofrag_Studies, study_PIDs)

## Replace NaNs and NAs with 0

for (i in 1:length(species_Matrix)) {
    species_Matrix[[i]][is.na(species_Matrix[[i]])] = 0
}

## Pivot data to long format

species_Matrix = species_Matrix %>% purrr:::map(~ pivot_longer(data = .x, !plot, names_to = "Species", values_to = "Abundance")) %>%
                        purrr:::map(mutate, plot = as.character(plot))

## Add PID to each list object
species_Matrix = mapply(function(x, y) add_column(x, PID = as.character(y)), species_Matrix, study_PIDs, SIMPLIFY = F)

## Bind together and add column for occurrence, which is just 0 or 1
species_Matrix = species_Matrix %>% bind_rows() %>% mutate(Occurrence = ifelse(Abundance > 0, 1, 0)) %>% relocate(PID, plot, Species, Abundance, Occurrence)

## Add column stating whether a study has abundance data or just pres/abs
species_Matrix = species_Matrix %>% group_by(PID) %>% mutate(Data = ifelse(max(Abundance) > 1 & sum(Abundance)%%1 == 0, "Abund", "Occ"))

## Get study PIDs for each broad taxa type (birds, mammals, bats, herps) so we can filter before running models

bird_Study_PIDs = str_split(grep(pattern = "Bird", x = biofrag_Studies, ignore.case = T, value = T), pattern = "_", simplify = T)[,1]
mammal_Study_PIDs = str_split(grep(pattern = "Mammal", x = biofrag_Studies, ignore.case = T, value = T), pattern = "_", simplify = T)[,1]
bat_Study_PIDs = str_split(grep(pattern = "Bat", x = biofrag_Studies, ignore.case = T, value = T), pattern = "_", simplify = T)[,1]
herps_Study_PIDs = str_split(grep(pattern = "Herp|Rept|Lizard|Amphib|Frog", x = biofrag_Studies, ignore.case = T, value = T), pattern = "_", simplify = T)[,1]
arthropods_Study_PIDs = str_split(grep(pattern = "Arthro|Beetle|_Ant|Butterfl|Bee|Moth|Arach|Fly|Flies|Hopper|Spider|insect", x = biofrag_Studies, ignore.case = T, value = T), pattern = "_", simplify = T)[,1]
other_Study_PIDs = study_PIDs %>% .[!. %in% c(bird_Study_PIDs, mammal_Study_PIDs, bat_Study_PIDs, herps_Study_PIDs, arthropods_Study_PIDs)]

## Filter out PIDs for studies that use taxa i'm not interested in - e.g Trees, Gastropods, Annelids

species_Matrix = species_Matrix %>% filter(!PID %in% other_Study_PIDs)

## Load in forest cover data

pland_Metrics = read_csv(paste0(gpath, "Data/percentage_Forest_Cover.csv"), locale = locale(encoding = "windows-1252"))

## Add forest cover for each plot to species matrix
## Remove unwanted columns (e.g percentange_inside)
## Filter to only keep matrix remove class 1 which is forest

species_Occurrence = left_join(species_Matrix, pland_Metrics, by = c("plot" = "plot_id", "PID" = "PID")) %>% 
                            filter(class == 0) %>% dplyr:::select(-class, -metric, -percentage_inside) %>%
                            relocate(PID, plot, buffer, Species, value, Abundance, Occurrence) %>% rename("Matrix_Cover" = "value", "Plot" = "plot", "Buffer" = "buffer")

if(!file.exists(paste0(gpath, "Data/species_Occurrence.csv"))) {
    write.csv(species_Occurrence, file = paste0(gpath, "Data/species_Occurrence.csv"), row.names = F)
}


