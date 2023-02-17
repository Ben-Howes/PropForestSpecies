################################################################
## Calculate Forest Cover for each site at different buffers
## Using 70% Forest cover as binary forest/non-forest cut off

## Ben Howes, 2nd Year PhD, Imperial College London, 2021
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
library(parallel)
library(terra)

gpath = "/home/ben/Documents/PhD/matrix_Response/"
setwd(gpath)

## List all biofrag study folders

biofrag_Studies = list.files("Data/BIOFRAG_Jan2020Status", pattern = "PID*")

## All biofrag folders have the same folders:
## Landscape meta, plots, species matrix, plot meta, and species meta
## Plot meta contains the point for each plot/site
## Landscape meta contains the central points for the study

####################################################
## Load in points for each plot for each study
####################################################

study_Plot_Points = lapply(biofrag_Studies, function(x) {
    study_Code = str_split(x, pattern = "_")[[1]][1]
    df = read.csv(paste0(gpath, "Data/BIOFRAG_Jan2020Status/", x, "/", study_Code, "_Plot.csv"))
    if(ncol(df) < 4) {df$PID = study_Code}
    colnames(df) = c("Plot", "Latitude", "Longitude", "PID")
    df = df %>% dplyr:::select(Plot, Latitude, Longitude, PID) %>% mutate(Plot = as.character(Plot), Latitude = as.numeric(Latitude), 
                                                                            Longitude = as.numeric(Longitude), PID = as.character(PID))
    return(df)})

## Remove rows with NAs

study_Plot_Points = lapply(study_Plot_Points, function(x) x[complete.cases(x),])

## Convert to geometry objects

study_Plot_Points_sf = lapply(study_Plot_Points, function(x) {
    df = st_as_sf(x, coords = c("Longitude", "Latitude"), crs = 4326)
    return(df)})

############################################################
## Create shapefiles for the different study sites
############################################################

## Make each study into a simple four sided shape file
## each corner is the maximum longitude and latitude in a direction
## and then an additional 0.1 degrees is added, which is about 10km
## to ensure that we have all data when buffers are then added to sites

Create_Study_Shapefile = function(study) {

    min_Long = min(study$Longitude) - 1; max_Long = max(study$Longitude) + 1;
    min_Lat = min(study$Latitude) - 1; max_Lat = max(study$Latitude) + 1;

    study_Coords = rbind(c(min_Long, min_Lat), c(max_Long, min_Lat), c(max_Long, max_Lat),
        c(min_Long, max_Lat), c(min_Long, min_Lat))

    study_Shapefile = st_geometry(st_polygon(list(study_Coords)))
    st_crs(study_Shapefile) = 4326
    
    return(study_Shapefile)
}

## Create shapefiles encompassing all plots of every study

study_Shapefiles = lapply(study_Plot_Points, Create_Study_Shapefile)

############################################################
## Calculate forest cover
############################################################

## Load in rasters

raster_Paths = list.files(paste0(gpath, "/Data/BIOFRAG_Jan2020Status"), pattern = "hansen_image_", recursive = T, full.names = T)

## Calculate Landscape Metrics for each plot in each study

calculate_Frag_Metrics = function(site, raster_File, study_Shapefile, buffer_Size) {

    site_Raster = stack(raster_File)[[1]] ## Load raster and select treecover2000 layer
    study_PID = site$PID[1]

    ## Set up AEQD projection format 

    azimuthal_Proj_Fmt = "+proj=aeqd +lat_0=%0.8f +lon_0=%0.8f +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m no_defs"

    if (!file.exists(paste0(str_split(raster_File, pattern = "hansen_")[[1]][1], "landscape_metrics_", study_PID, "_", buffer_Size, ".csv"))) {

        print(c(study_PID, buffer_Size))

        ## Only keep pixels which have forest cover equal or greather than 70% 
        ## this is a general rule but we could change it to see how this threshold affects results

        tree = stack(raster_File)[[1]]
        mask = stack(raster_File)[[5]]

        ## Change no data and water cells to NA in the mask
        mask[mask == 2 | mask == 0] = NA

        ## Update values in tree to be NA if they are in NA of the mask
        tree = tree*mask

        ## Classify cells 70 or above to 1, and 69 or below as 0, leaving sea and water as NA values
        tree[tree < 70 & !is.na(tree)] = 0
        tree[tree >= 70] = 1

        forest_Raster = tree

        # behrmann_Proj = "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m"

        ## Find centre of study shapefile to then set AEQD to specific location

        centre_Coords = study_Shapefile %>% st_centroid() %>% st_coordinates()

        ## Set AEQD specific for this site

        azimuthal_Proj = sprintf(azimuthal_Proj_Fmt, centre_Coords[,'Y'], centre_Coords[,'X'])

        ## Transform sites and Raster to AEQD projection

        site = st_transform(site, azimuthal_Proj)

        # forest_Raster = projectRaster(forest_Raster, crs = crs(site), method = 'ngb', res = 30)
        forest_Raster1 = rast(forest_Raster)
        forest_Raster = terra:::project(x = forest_Raster1, y = crs(site), method = 'near', res = 30)

        # Calculate landscape metrics
    
        frag_Metrics = sample_lsm(forest_Raster, site, plot_id = site$Plot, shape = "circle", size = buffer_Size, all_classes = T, 
            what = "lsm_c_pland", progress = T)

        ## Remove NA classes, then turn NAs in value column to 0

        frag_Metrics = frag_Metrics[!is.na(frag_Metrics$class), ]
        frag_Metrics[is.na(frag_Metrics$value), ]$value = 0

        write.csv(frag_Metrics, file = paste0(str_split(raster_File, pattern = "hansen_")[[1]][1],
                                            "landscape_metrics_", study_PID, "_", buffer_Size, ".csv"), row.names = FALSE)
    }
}

mcmapply(calculate_Frag_Metrics, study_Plot_Points_sf, raster_Paths, study_Shapefiles, rep(c(50, 200, 2000), each = length(raster_Paths)), mc.cores = 1)

## Load in plot landscape metrics
## 0 = matrix, 1 = forest
## WILL HAVE TO ADD BUFFER SIZE INTO THIS SECTION AND GET IT FROM THE FILE NAME, ACTUALLY CAN JUST ONLY SEARCH FILES WITH 50, 200 AND 2000 IN THEM
lsm_Paths = mixedsort(list.files(paste0(gpath, "/Data/BIOFRAG_Jan2020Status"), pattern = "landscape_metrics_PID", recursive = T, full.names = T))
lsm_Paths = grep(pattern = "_50.c|_200.c|_2000.c", lsm_Paths, value = T)

plot_Landscape_Metrics = pblapply(lsm_Paths, function(x) {read.csv(x) %>% dplyr:::select(class, metric, value, plot_id, percentage_inside) %>% 
                                                            add_column(PID = str_split(str_split(x[[1]], pattern = "/")[[1]][[10]], pattern = "_")[[1]][1])})

## Add buffer column

buffer_List = c(50,200,2000)

plot_Landscape_Metrics = mapply(function(x,y) add_column(x, buffer = y), plot_Landscape_Metrics, buffer_List, SIMPLIFY = F)

pland_Metrics = plot_Landscape_Metrics %>% purrr:::map(mutate, plot_id = as.character(plot_id)) %>% bind_rows() %>% filter(metric == "pland")

if(!file.exists(paste0(gpath, "Data/percentage_Forest_Cover.csv"))) {
    write.csv(pland_Metrics, file = paste0(gpath, "Data/percentage_Forest_Cover.csv"), row.names = F)
}

