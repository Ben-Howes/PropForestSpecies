######################################################################
## Get range values for amphibians, mammals, and repitles
## Already have bird range values from avonet
######################################################################

library(tidyverse)
library(janitor)
library(sf)
library(taxize)
library(rredlist)
library(raster)
library(terra)
library(parallel)

gpath = "/home/ben/Documents/PhD/matrix_Response/"
setwd(gpath)

## Get paths to IUCN shapefiles
amphib_path = "/home/ben/Documents/PhD/matrix_Response/Data/Ranges/AMPHIBIANS/AMPHIBIANS.shp"
mammal_path = "/home/ben/Documents/PhD/matrix_Response/Data/Ranges/MAMMALS/MAMMALS.shp"
reptile_path = "/home/ben/Documents/PhD/matrix_Response/Data/Ranges/REPTILES/REPTILES.shp"

## Get habitat data for each taxa group

get_habitat = function(species, taxa) {
  print(species)
  hab = rl_habitats(gsub("_", " ", species))$result
  if(length(hab) > 0) {
    out = hab %>% mutate(species = species, taxa = taxa, .before = 1) %>% dplyr::select(-code) 
    return(out)
  }
  Sys.sleep(1)
}

amphibians = st_read(amphib_path)
if(!file.exists(paste0(gpath, "Data/Ranges/amphibian_habitats.csv"))) {
  amphibian_species = unique(amphibians$binomial)
  amphibian_habitats = lapply(amphibian_species, get_habitat, taxa = "amphibians") %>% bind_rows()
  write.csv(amphibian_habitats, paste0(gpath, "Data/Ranges/amphibian_habitats.csv"))
} else {
  amphibian_habitats = read.csv(paste0(gpath, "Data/Ranges/amphibian_habitats.csv"))
}

mammals = st_read(mammal_path)
if(!file.exists(paste0(gpath, "Data/Ranges/mammal_habitats.csv"))) {
  mammals_species = mammals %>% distinct(species) %>% first()
  mammals_habitats = lapply(mammals_species, get_habitat, taxa = "mammals") %>% bind_rows()
  write.csv(mammals_habitats, paste0(gpath, "Data/Ranges/mammal_habitats.csv"))
} else {
  mammal_habitats = read.csv(paste0(gpath, "Data/Ranges/mammal_habitats.csv"))
}

reptiles = st_read(reptile_path)
if(!file.exists(paste0(gpath, "Data/Ranges/reptile_habitats.csv"))) {
  reptiles_species = reptiles %>% distinct(species) %>% first()
  reptiles_habitats = lapply(reptiles_species, get_habitat, taxa = "reptiles") %>% bind_rows()
  write.csv(reptiles_habitats, paste0(gpath, "Data/Ranges/reptile_habitats.csv"))
} else {
  reptile_habitats = read.csv(paste0(gpath, "Data/Ranges/reptile_habitats.csv"))
}

birds = st_read(paste0(gpath, "Data/Ranges/BIRDS/BOTW.gdb"))
if(!file.exists(paste0(gpath, "Data/Ranges/bird_habitats.csv"))) {
  bird_species = birds$sci_name %>% unique()
  bird_habitats = lapply(bird_species, get_habitat, taxa = "birds") %>% bind_rows()
  write.csv(bird_habitats, paste0(gpath, "Data/Ranges/bird_habitats.csv"))
} else {
  bird_habitats = read.csv(paste0(gpath, "Data/Ranges/bird_habitats.csv")) %>% mutate(species = tolower(gsub(" ", "_", species)))
}

## 162 ranges for bird data are stored as multisurface objects which
## cannot be used in later functions e.g rasterize
## So we need to change them to multipolygons using the gdal method at:
## https://github.com/r-spatial/sf/issues/1573 

birds = st_cast(birds, to = "MULTIPOLYGON")

# fix_multisurfaces = function(x) {

#   type = st_geometry_type(x)

#   if(type == "MULTISURFACE") {
#     st_write(x, file.path(tempdir(), "temp.gpkg"), append = FALSE)
#     #> Writing layer `object3' to data source `/tmp/RtmpnKprFn/object3.gpkg' using driver `GPKG'
#     #> Writing 1 features with 0 fields and geometry type Curve Polygon.

#     gdal_utils("vectortranslate", 
#               file.path(tempdir(), "temp.gpkg"), 
#               file.path(tempdir(), "object3_1.gpkg"), 
#               options=c("-nlt", "CONVERT_TO_LINEAR"))

#     out = read_sf(file.path(tempdir(), "object3_1.gpkg"))

#     file.remove(file.path(tempdir(), "object3_1.gpkg"))

#     return(out$geom)
#   } else {
#     return(x)
#   }
# }

## In addition to multisurface errors, some bird polygons are invalid because they have polygons that overlap
## the meridian. This doesn't cause issues now but does once we transform to equal area
## So we need to turn S2 off, then identify these invalid polygons and crop their extent by 0.1 degree
## to remove the overlapping part of the polygon

sf_use_s2(FALSE)

## First remove species which are extinct, or have non-bredding seasonal ranges
## to minimise the time the next functions take to run
birds = birds %>% filter(presence == 1 & origin == 1 & seasonal %in% c(1,2))

## Also remove all marine species, these won't be included in the analysis and will again
## reduce the time taken for the next functions
## Marine species removed by categorising any species that just uses habitats including the word marine as marine
marine_bird_habitats = bird_habitats %>% group_by(species) %>% 
  summarise(marine_habitats = length(grep("marine", habitat, ignore.case = T))) %>% 
  mutate(marine = ifelse(marine_habitats > 0, "true", "false")) %>% dplyr::select(-marine_habitats)

birds = birds %>% left_join(marine_bird_habitats, by = c("sci_name" = "species")) %>% filter(marine != "true")

## Make geometries valid and then crop out the edges of xmax and xmin so they don't overlap and circle back around
## This takes hours to run
birds = birds %>% st_make_valid() %>% st_crop(xmin = -179.9, ymin = -90, xmax = 179.9, ymax = 90)

## Combine all habitat data together
taxa_habitats = bind_rows(amphibian_habitats, mammal_habitats, reptile_habitats, bird_habitats) %>% dplyr::select(species, taxa, habitat, suitability)

## Categorise species as forest or non-forest specialist, or generalist ()uses forest and non-forest
taxa_habitats = taxa_habitats %>% group_by(species) %>%
  filter(suitability == "Suitable") %>% mutate(species = tolower(gsub(" ", "_", species))) %>%
  summarise(forest_habitats = length(grep("forest", habitat, ignore.case = T) %>% 
                                       grep("degraded", ., ignore.case = T, invert = T)), 
            nonforest_habitats = length(grep("forest", habitat, ignore.case = T, invert = T)), 
            total_habitats = n(), taxa = taxa[[1]]) %>% mutate(forest = ifelse(forest_habitats == total_habitats, 1, 
                                      ifelse(nonforest_habitats == total_habitats, 0, "generalist")))

## Since the birdlife dataset has a few different column names and some missing column names we need to add them
## Firstly change species names to binomial and shape to geometry to match IUCN
## Then need to add marine column to remove marine species
## this will be done by categorising any species that just uses habitats including the word marine as marine
marine_bird_habitats = bird_habitats %>% group_by(species) %>% 
  summarise(marine_habitats = length(grep("marine", habitat, ignore.case = T))) %>% 
  mutate(marine = ifelse(marine_habitats > 0, "true", "false")) %>% dplyr::select(-marine_habitats)

birds = birds %>% rename(binomial = sci_name, geometry = Shape, 
  SHAPE_Leng = Shape_Length, SHAPE_Area = Shape_Area) %>% 
  left_join(marine_bird_habitats, by = c("binomial" = "species"))

taxa = bind_rows(amphibians, mammals, reptiles, birds) %>% mutate(species = tolower(gsub(" ", "_", binomial))) %>% 
  left_join(taxa_habitats) %>% filter(marine != "true" & !is.na(forest)) %>% 
  dplyr::select(binomial, species, taxa, presence, origin, seasonal, marine, SHAPE_Leng, SHAPE_Area, 
                forest_habitats, nonforest_habitats, total_habitats, forest, geometry)

## Remove large groups of shapefiles that are no longer needed
## to free up RAM
rm(reptiles, amphibians, mammals, birds)
gc()

## Keep ranges in which the species are definitely extant (presence code 1)
## Native (origin code 1) and Resident (season code 1)
taxa = taxa %>% filter(presence == 1 & origin == 1 & seasonal %in% c(1,2))

## Transform shapefiles to behrmann equal area projection
behr = "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
taxa = st_transform(taxa, behr)

## Make template raster using terra and behrmann equal area projection, each cell is now 93694.99m (93.7km x 93.7km)
template_raster = rast() %>% project(behr)

## Calculate the species richness of forest vs nonforest species for each taxa
## Do this in a function in which we have to convert sf objects to SpatVectors
## Which can then be rasterized using the template raster before using mosaic
## to sum together the cells
## Touches = TRUE is used in rasterization whcih means a species is counted as being in a cell
## if just one part of the polygon overlaps the cell (doesn't have to cross the center of the cell)

taxa = taxa %>% mutate(type = ifelse(forest == 1, 1, 0))

species_richness_latlong = function(dat) {

  print(paste(dat$taxa[[1]], dat$type[[1]]))

  vect = as(dat, "SpatVector")
  rasts = lapply(1:nrow(vect), function(x) terra::rasterize(vect[x,], template_raster, touches = TRUE))
  list = sprc(rasts)
  map = terra::mosaic(list, fun = "sum")

  out = list(taxa = dat$taxa[[1]], type = dat$type[[1]], map = map)
  return(out)

}

richness = taxa %>% group_split(taxa, type) %>% map(., ~species_richness_latlong(.))

total_amphibians = terra::mosaic(richness[[1]]$map, richness[[2]]$map, fun = "sum")
total_birds = terra::mosaic(richness[[3]]$map, richness[[4]]$map, fun = "sum")
total_mammals = terra::mosaic(richness[[5]]$map, richness[[6]]$map, fun = "sum")
total_reptiles = terra::mosaic(richness[[7]]$map, richness[[8]]$map, fun = "sum")
prop_amphibians = richness[[2]]$map/total_amphibians
prop_birds = richness[[4]]$map/total_birds
prop_mammals = richness[[6]]$map/total_mammals
prop_reptiles = richness[[8]]$map/total_reptiles
all_list = list(total_amphibians, total_mammals, total_reptiles, total_birds, prop_amphibians, prop_mammals, prop_reptiles, prop_birds)
names = c("total_amphibians", "total_mammals", "total_reptiles", "total_birds", "prop_amphibians", "prop_mammals", "prop_reptiles", "prop_birds")

richness_rasts = rast(all_list)
names(richness_rasts) = names

writeRaster(richness_rasts, paste0(gpath, "Data/Ranges/richness_rasts_resident_breeding.tif"), overwrite = TRUE)
