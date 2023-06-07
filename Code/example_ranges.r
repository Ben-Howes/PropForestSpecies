###################################
## Make example ranges
## and habitat classifications
## for presentations or SI
###################################

library(tidyverse)
library(sf)
library(terra)
library(tidyterra)

gpath = "/home/ben/Documents/PhD/PropForestSpecies/"
setwd(gpath)

mammal_path = "/home/ben/Documents/PhD/PropForestSpecies/Data/Ranges/MAMMALS/MAMMALS.shp"

## Load in mammal ranges and habitats
mammals = st_read(mammal_path)
if(!file.exists(paste0(gpath, "Data/Ranges/mammal_habitats.csv"))) {
  mammals_species = mammals %>% distinct(species) %>% first()
  mammals_habitats = lapply(mammals_species, get_habitat, taxa = "mammals") %>% bind_rows()
  write.csv(mammals_habitats, paste0(gpath, "Data/Ranges/mammal_habitats.csv"))
} else {
  mammal_habitats = read.csv(paste0(gpath, "Data/Ranges/mammal_habitats.csv"))
}

mammal_habitats = mammal_habitats %>% dplyr::select(species, taxa, habitat, suitability)

## Categorise species as forest or non-forest specialist, or generalist ()uses forest and non-forest
mammal_habitats = mammal_habitats %>% group_by(species) %>%
  filter(suitability == "Suitable") %>% mutate(species = tolower(gsub(" ", "_", species))) %>%
  summarise(forest_habitats = length(grep("forest", habitat, ignore.case = T) %>% 
                                       grep("degraded", ., ignore.case = T, invert = T)), 
            nonforest_habitats = length(grep("forest", habitat, ignore.case = T, invert = T)), 
            total_habitats = n(), taxa = taxa[[1]]) %>% mutate(forest = ifelse(forest_habitats == total_habitats, 1, 0))

mammals = mammals %>% mutate(species = tolower(gsub(" ", "_", binomial))) %>% 
  left_join(mammal_habitats) %>% filter(marine != "true" & !is.na(forest)) %>% 
  dplyr::select(binomial, species, genus, family, order_, taxa, presence, origin, seasonal, marine, SHAPE_Leng, SHAPE_Area, 
                forest_habitats, nonforest_habitats, total_habitats, forest, geometry)

## Keep ranges in which the species are definitely extant (presence code 1)
## Native (origin code 1) and Resident (season code 1)
mammals = mammals %>% filter(presence == 1 & origin == 1 & seasonal %in% c(1,2))

## Filter by taxa of interest, let's do a Lemur family
lemur = mammals %>% filter(family == toupper("Lemuridae"))

## Transform shapefiles to behrmann equal area projection
behr = "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
lemur = st_transform(lemur, behr)

## Join ranges of the same species together
lemur = lemur %>% group_by(species, forest) %>% summarise(range = st_union(geometry))

## Make template raster using terra and behrmann equal area projection, each cell is now 93694.99m (93.7km x 93.7km)
template_raster = rast() %>% project(behr)

## Load in coast
coast = st_read("../Raw_Data/gshhg-shp-2.3.7/GSHHS_shp/c/GSHHS_c_L1.shp")
coast = st_transform(coast, behr)
coast = st_crop(coast, st_bbox(lemur))
coast_outline = st_cast(coast, "MULTILINESTRING")

## Plot ranges
ggplot() + 
    geom_sf(data = coast, fill = "grey95") +
    geom_sf(data = coast_outline, linewidth = 0.25) +
    geom_sf(data = lemur, aes(fill = species)) +
    theme_void() +
    scale_fill_viridis_d() +
    facet_wrap(. ~ species, ncol = 7) +
    theme(legend.position = "none") +
    theme(strip.background = element_blank(),
    strip.text.x = element_blank())

ggsave(paste0(gpath, "Data/Ranges/example_lemur_ranges.png"))

ggplot() + 
    geom_sf(data = coast, fill = "grey95") +
    geom_sf(data = coast_outline, linewidth = 0.25) +
    geom_sf(data = lemur, aes(fill = as.factor(forest))) +
    theme_void() +
    scale_fill_viridis_d() +
    facet_wrap(. ~ species, ncol = 7) +
    theme(legend.position = "none") +
    theme(strip.background = element_blank(),
    strip.text.x = element_blank())
    
ggsave(paste0(gpath, "Data/Ranges/example_lemur_specialisations.png"))

## Make example raster maps
species_richness_latlong = function(dat) {

  vect = as(dat, "SpatVector")
  rasts = lapply(1:nrow(vect), function(x) terra::rasterize(vect[x,], template_raster, touches = TRUE))
  list = sprc(rasts)
  map = terra::mosaic(list, fun = "sum")

  out = list(forest = dat$forest[[1]], map = map)
  return(out)

}

richness = lemur %>% ungroup() %>% group_split(forest) %>% map(., ~species_richness_latlong(.))

total_lemur = terra::mosaic(richness[[1]]$map, richness[[2]]$map, fun = "sum")
prop_lemur = richness[[2]]$map/total_lemur

ggplot() +
    geom_spatraster(data = crop(total_lemur, coast)) +
    theme_void() +
    scale_fill_viridis_c(na.value = "transparent") +
    labs(fill = "Species\nRichness") +
    theme(text = element_text(size = 25))

ggsave(paste0(gpath, "Data/Ranges/example_total_species.png"))

ggplot() +
    geom_spatraster(data = crop(prop_lemur, coast)) +
    theme_void() +
    scale_fill_viridis_c(na.value = "transparent") +
    labs(fill = "Proportion Forest\nSpecies") +
    theme(text = element_text(size = 25))

ggsave(paste0(gpath, "Data/Ranges/example_prop_forest_raster.png"))
