##################################
## Find number of studies in biofrag with abundance data
## Find nummber of studies with matrix dominant sites
###################################

library(tidyverse)
library(janitor)
library(lemon)
library(sf)

gpath = "/home/ben/Documents/PhD/matrix_Response/"
setwd(gpath)

## List all biofrag study folders
biofrag_Studies = list.files(paste0(gpath, "Data/BIOFRAG_Jan2020Status"), pattern = "PID*")
## 114 studies

study_PIDs = sapply(1:length(biofrag_Studies), function (x) {str_split(biofrag_Studies, pattern = "_")[[x]][1]})

## Check if data is abundance or pres/abs
data_Type = mapply(function(x,y) {
    read.csv(paste0(gpath, "Data/BIOFRAG_Jan2020Status/", x, "/", y, "_species_matrix.csv")) %>% clean_names() %>% 
        select(-plot) %>% pivot_longer(everything(), names_to = "species", values_to = "occ") %>% 
            slice_max(occ, with_ties = F) %>% mutate(type = ifelse(occ == 1, "occ", "abund"), PID = y)
}, biofrag_Studies, study_PIDs, SIMPLIFY = F) %>% bind_rows()

## Add taxa type
taxa = read_csv(paste0(gpath, "Data/taxa_PID_Lookup.csv"))
data_Type = left_join(data_Type, taxa) %>% clean_names() %>% select(type, pid, taxa)
plot_Type = data_Type %>% group_by(taxa, type) %>% summarise(n = n())

data_Type_Plot = ggplot(plot_Type, aes(taxa, n, fill = type)) + geom_col(position = "dodge") + theme_classic() + 
    labs(x = "Taxa", y = "Count", fill = "Data Type") + theme(text = element_text(size = 20)) + theme(text = element_text(size = 20)) +
        scale_y_continuous(expand = expansion(mult = c(0, .1)))

ggsave(filename = "Results/Exploratory/data_Type_Plot.png", plot = data_Type_Plot, device = "png", width = 20, height = 10, units = "in")

## 99/104 have abund data

### How many studies have matrix dominant sites?
## We want them to be in matrix at local (50m) and slightly larger (200m) scale
## So always above a large % matrix at local, say 100% for now and another % at 200m
##  for now use 3 percentages that could be "matrix dominant", over 50%, over 75% and over 90%

## Load in site data
matrix_Cov = read_csv(paste0(gpath, "Data/percentage_Forest_Cover.csv")) %>% left_join(taxa) %>% 
    clean_names() %>% filter(buffer %in% c("50", "200") & class == 0) %>% select(-metric, -percentage_inside)

### Remove studies which do not have abundance data
abund_PID = data_Type %>% filter(type == "abund")
matrix_Cov = matrix_Cov %>% filter(pid %in% abund_PID$pid)

## Join 50 and 200 so they are in one single row for filtering
matrix_Cov50 = matrix_Cov %>% filter(buffer == 50)
matrix_Cov200 = matrix_Cov %>% filter(buffer == 200)

matrix_Wider = left_join(matrix_Cov50, matrix_Cov200, by = c("class", "plot_id", "pid", "taxa"), suffix = c("_50", "_200")) %>% 
    select(-buffer_50, -buffer_200)

## Filter so only sites with 100% matrix at the local scale remain
matrix_Local = matrix_Wider %>% filter(value_50 == 100)

matrix_Local %>% distinct(pid) %>% nrow()
# 68 of 99 abund studies have sites in 100% local (50m) matrix#

## How many studies have matrix sites at the > 50% level at 200m?
matrix_Local %>% filter(value_200 > 50) %>% distinct(pid) %>% nrow()
# 68 of 99 abund studies have sites in 100% local (50m) matrix and > 50% 200m scale

## How many studies have matrix sites at the > 75% level at 200m?
matrix_Local %>% filter(value_200 > 75) %>% distinct(pid) %>% nrow()
# 66 of 99 abund studies have sites in 100% local (50m) matrix and > 75% 200m scale

## How many studies have matrix sites at the > 90% level at 200m?
matrix_Local %>% filter(value_200 > 90) %>% distinct(pid) %>% nrow()
# 61 of 99 abund studies have sites in 100% local (50m) matrix and > 90% 200m scale

## Looking deeper at the 90% level
matrix_90 = matrix_Local %>% filter(value_200 > 90)
matrix_90_Plot = matrix_90 %>% distinct(pid, taxa) %>% group_by(taxa) %>% summarise(n = n())

matrix_Studies_Plot = ggplot(matrix_90_Plot, aes(taxa, n)) + geom_col() + theme_classic() + 
    labs(x = "Taxa", y = "Count", title = "Studies with at least one \n predominantly matrix site \n 
        (100% local, >90% )") + theme(text = element_text(size = 20)) +
            scale_y_continuous(expand = expansion(mult = c(0, .1)))

ggsave(filename = "Results/Exploratory/matrix_Studies_Plot.png", plot = matrix_Studies_Plot, device = "png", width = 20, height = 10, units = "in")

## How many sites does each study have?
matrix_Sites = matrix_90 %>% group_by(taxa, pid) %>% summarise(sites = n())

matrix_Sites_Plot = ggplot(matrix_Sites, aes(sites, fill = taxa)) + geom_histogram(binwidth = 5) + theme_classic() +
    labs(y = "Frequency", x = "Number of Matrix Sites", fill = "Taxa") + theme(text = element_text(size = 20)) +
        facet_rep_wrap(. ~ taxa) + scale_y_continuous(expand = expansion(mult = c(0, .1)))

ggsave(filename = "Results/Exploratory/matrix_Sites_Plot.png", plot = matrix_Sites_Plot, device = "png", width = 20, height = 10, units = "in")

## Load in spatial data so we can see if we have a good latitudinal spread of studies
geoms = load(paste0(gpath, "Data/study_Geometries.Rdata"))
study_Centroids = lapply(study_Plot_Points, function(x) data.frame(PID = x$PID[[1]], Latitude = mean(x$Latitude), Longitude = mean(x$Longitude))) %>% bind_rows()

### Join taxa, filter to only include abundance studies and those with at least one predominantly matrix site
location_Data = study_Centroids %>% left_join(taxa) %>% filter(PID %in% matrix_90$pid) %>% st_as_sf(coords = c("Longitude", "Latitude"))

world = map_data("world")
study_Map = ggplot() + geom_map(data = world, map = world, aes(long, lat, map_id = region))  + geom_sf(data = location_Data, aes(col = Taxa), size = 5) +
    theme_classic() + labs(x = "Longitude", y = "Latitude") + theme(text = element_text(size = 20))

ggsave(filename = "Results/Exploratory/study_Map.png", plot = study_Map, device = "png", width = 20, height = 10, units = "in")meter