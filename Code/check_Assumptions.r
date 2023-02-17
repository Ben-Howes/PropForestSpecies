#################################################################
## Check whether the % matrix sites differs with latitude
## and whether we have more occurrences in higher latitudes
#################################################################

library(tidyverse)
library(janitor)
library(sf)
library(lemon)

gpath = "/home/ben/Documents/PhD/matrix_Response/"
setwd(gpath)

## Load data
spec = read_csv(paste0(gpath, "Data/species_Occurrence.csv")) %>% 
    clean_names() %>%
    mutate(matrix_cover = matrix_cover/100) %>% 
        filter(buffer == 200)

## Get data for just sites to see how they chnange with latitude
sites = spec %>% distinct(pid, plot, matrix_cover) %>%
    mutate(matrix = ifelse(matrix_cover >= 0.9, 1, 0))

matrix_studies = sites %>% filter(matrix == 1) %>% distinct(pid)

## Load latitudinal data
geoms = load(paste0(gpath, "Data/study_Geometries.Rdata"))
study_centroids = lapply(study_Plot_Points, function(x)
    data.frame(PID = x$PID[[1]], Latitude = mean(x$Latitude), 
        Longitude = mean(x$Longitude))) %>% bind_rows() %>%
            clean_names()

duplicate_studies = c("PID0100", "PID0101", "PID0102", "PID0103", "PID0104", "PID0105",
                      "PID0106", "PID0107", "PID0108", "PID0109", "PID0110", "PID0111", "PID0113", "PID0114", "PID0117", "PID0119", "PID0124",
                      "PID0138", "PID1002", "PID1005")

sites = sites %>% left_join(study_centroids) %>% mutate(abs_lat = abs(latitude)) %>% 
        filter(!pid %in% duplicate_studies)

ggplot(sites %>% filter(pid %in% matrix_studies$pid) %>%
    group_by(pid, latitude, abs_lat) %>%
        summarise(prop_matrix = sum(matrix)/n()),
            aes(abs_lat, prop_matrix)) + geom_point(size = 2) + theme_classic() +
                labs(x = "Absolute Latitude", y = "Proportion Matrix Sites") +
                    geom_smooth(method = "lm") +
                        theme(text = element_text(size = 20))

## What do the occurrences of species look like across latitudes
## Worried we have more occurrences, and therefore higher chance of finding matrix species in higher latitudes # nolint

taxa = read.csv(paste0(gpath, "Data/taxa_PID_Lookup.csv")) %>% clean_names()

spec_freq = spec %>% filter(pid %in% matrix_studies$pid & !pid %in% duplicate_studies) %>%
    left_join(study_centroids) %>% 
        mutate(abs_lat = abs(latitude)) %>% 
            group_by(pid, species, latitude, abs_lat) %>% 
                summarise(count = sum(occurrence)) %>% 
                    left_join(taxa)

study_freq = spec_freq %>% group_by(pid, latitude, abs_lat) %>% 
    mutate(over5 = ifelse(count > 5, 1, 0), 
        over10 = ifelse(count > 10, 1, 0)) %>%
            summarise(prop5 = sum(over5)/n(), prop10 = sum(over10)/n())

ggplot(study_freq, aes(x = abs_lat, y = prop5)) + 
    geom_point(size = 2) + theme_classic() +
        labs(x = "Absolute Latitude", 
            y = "Proportion of species occuring > 5 times") +
                theme(text = element_text(size = 20))

ggplot(study_freq, aes(x = abs_lat, y = prop10)) + 
    geom_point(size = 2) + theme_classic() + 
        labs(x = "Absolute Latitude", 
            y ="Proportion of species occuring > 10 times") +
                theme(text = element_text(size = 20))

ggplot(filter(spec_freq, taxa == "Birds"), aes(count)) + geom_density() +
    facet_rep_wrap(. ~ round(abs_lat)) +
        theme_classic() + labs(x = "Count", y = "Density") +
            geom_vline(xintercept = 10, linetype = 2) +
                geom_vline(xintercept = 5, linetype = 3) +
                        xlim(0,50) +
                            theme(text = element_text(size = 20))
