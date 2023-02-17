##########################################################
## Plot proportion of forest vs matrix species
## from latitudes and habitat usage taken from 
## IUCN for mammals, amphibs, and reptiles
## and from AVONET for birds
##########################################################

library(tidyverse)
library(janitor)
library(sf)
library(lemon)
library(mgcv)
library(lme4)
library(ggeffects)
library(MuMIn)
library(terra)
library(maptools)
library(colorspace)
library(tidyterra)
library(glmmTMB)
library(DHARMa)

gpath = "/home/ben/Documents/PhD/matrix_Response/"
setwd(gpath)

## Make behrmann template raster for use in transforming in future
behr = "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
template_raster = rast() %>% project(behr)

## Load basic world outline from GSHHS to use in cropping/masking
coast = st_read("../Raw_Data/gshhg-shp-2.3.7/GSHHS_shp/c/GSHHS_c_L1.shp")
coast = st_transform(coast, behr)
coast_outline = st_cast(coast, "MULTILINESTRING")

## Load all latitude and habitat data
forest_prop_maps = rast(paste0(gpath, "Data/Ranges/richness_rasts_resident_breeding.tif"))

ggplot() + geom_spatraster(data = forest_prop_maps[[5:8]]) + geom_sf(data = coast_outline) + theme_classic() +
  scale_fill_viridis_c(na.value = "white") + labs(fill = "Proportion Forest Species") +
  facet_wrap(~lyr, ncol = 2) 

################################################
## Are these trends explained by forest area?
################################################

## Load forest area data calculated from Hansen 2013
## Forest categorised at 30, 50, and 70% tree cover cut-offs
## For now just use 70% cut off values and later use others for sensitivity analysis

forest_rasts = rast(paste0(gpath, "Data/Ranges/forest_rasts.tif"))

forest_rasts = mask(forest_rasts, coast)

ggplot() + geom_spatraster(data = forest_rasts[[3]]) + geom_sf(data = coast_outline) + theme_classic() +
  scale_fill_viridis_c(na.value = "white") + labs(fill = "Proportion Forest Area")

## Add this raster to the stack of proportion forest species rasters
forest_prop_maps = c(forest_prop_maps, forest_rasts)

######################################################
## Could differences be due to historic forest loss?
######################################################
## Load in data from Betts et al tranformed into behr projection and % calculated at 1 cell
## area 96x96km
loss_rast = rast(paste0(gpath, "Data/Ranges/loss_rast.tif"))
loss_rast = mask(loss_rast, coast)

ggplot() + geom_spatraster(data = loss_rast) + geom_sf(data = coast_outline) + scale_fill_viridis_c(na.value = "white") +
  theme_classic() + labs(fill = "Proportion of Forest Loss")

## Add this raster to the stack of proportion forest species rasters
forest_prop_maps = c(forest_prop_maps, loss_rast)

######################################################
## Could differences be due to historic disturbances?
######################################################
## Load in data from Betts et al transformed into behrmann and categorised as 1 (disturbed) or 0 (undisturbed)
dist_rasts = rast(paste0(gpath, "Data/Ranges/dist_rasts.tif"))
dist_rasts = mask(dist_rasts, coast)

## Combine all disturbances into one raster
dist_rast = terra::mosaic(dist_rasts[[1]], dist_rasts[[2]], dist_rasts[[3]], fun = "max")
names(dist_rast) = "disturbances"

ggplot() + geom_spatraster(data = dist_rast) + geom_sf(data = coast_outline) + scale_fill_viridis_c(na.value = "white") +
  theme_classic() + labs(fill = "Disturbances")

## Add this raster to the stack of proportion forest species rasters
forest_prop_maps = c(forest_prop_maps, dist_rast)

######################################################
## How about how long an area has had that same type 
## of forest?
######################################################

eocene = st_read(paste0(gpath, "Data/Ranges/geological_forest_type/eocene.shp")) %>% st_transform("EPSG:4326")
oligocene = st_read(paste0(gpath, "Data/Ranges/geological_forest_type/oligocene.shp")) %>% st_transform("EPSG:4326")
miocene = st_read(paste0(gpath, "Data/Ranges/geological_forest_type/miocene.shp")) %>% st_transform("EPSG:4326")
pliocene = st_read(paste0(gpath, "Data/Ranges/geological_forest_type/pliocene.shp")) %>% st_transform("EPSG:4326")
holocene = st_read(paste0(gpath, "Data/Ranges/geological_forest_type/holocene.shp"))

eocene = project(as(eocene, 'SpatVector'), behr)
oligocene = project(as(oligocene, 'SpatVector'), behr)
miocene = project(as(miocene, 'SpatVector'), behr)
pliocene = project(as(pliocene, 'SpatVector'), behr)
holocene = project(as(holocene, 'SpatVector'), behr)

eocene = terra::rasterize(eocene, template_raster, field = "type")
oligocene = terra::rasterize(oligocene, template_raster, field = "type")
miocene = terra::rasterize(miocene, template_raster, field = "type")
pliocene = terra::rasterize(pliocene, template_raster, field = "type")
holocene = terra::rasterize(holocene, template_raster, field = "type")

forest_type = c(eocene, oligocene, miocene, pliocene, holocene)
forest_type = as.data.frame(forest_type, xy = T) %>% rename("eocene" = 3, "oligocene" = 4, "miocene" = 5, "pliocene" = 6, "holocene" = 7)
forest_type = forest_type %>% mutate(geological_forest = case_when(eocene == oligocene & oligocene == miocene & miocene == pliocene & pliocene == holocene ~ 5, 
oligocene == miocene & miocene == pliocene & pliocene == holocene ~ 4 , miocene == pliocene & pliocene == holocene ~ 3, 
pliocene == holocene ~ 2)) %>% mutate(geological_forest = ifelse(is.na(geological_forest), 1, geological_forest))

## Select just x y and geological forest, where the larger the value of geological forest the longer
## the forest in that area has continued to be the same forest type
forest_type_rast = rast(dplyr::select(forest_type, x,y,geological_forest), crs = behr)
forest_type_rast = mask(forest_type_rast, coast) 

ggplot() + geom_spatraster(data = forest_type_rast) +
  geom_sf(data = coast_outline) +
  theme_classic() + scale_fill_viridis_c(na.value = "white") + 
  labs(fill = "Geological Forest Time")

## Add this raster to the stack of proportion forest species rasters
forest_prop_maps = c(forest_prop_maps, forest_type_rast)

######################################################
## Analysis
######################################################

## Turn these rasters into a data frame for use in models
total_df = as.data.frame(forest_prop_maps, xy = T) %>% dplyr::select(-prop_amphibians, -prop_mammals, -prop_reptiles, -prop_birds) %>% 
  pivot_longer(cols = c(total_amphibians, total_mammals, total_reptiles, total_birds), names_to = "taxa", values_to = "n_spec") %>%
  mutate(taxa = sub("^[^_]*_", "", taxa))

prop_forest_df = as.data.frame(forest_prop_maps, xy = T) %>% dplyr::select(-total_amphibians, -total_mammals, -total_reptiles, -total_birds) %>%
  pivot_longer(cols = c(prop_amphibians, prop_mammals, prop_reptiles, prop_birds), names_to = "taxa", values_to = "prop_forest") %>%
  mutate(taxa = sub("^[^_]*_", "", taxa))

prop_forest_df = prop_forest_df %>% left_join(total_df) %>% relocate(x,y,prop_forest,n_spec,prop_forest_area,taxa)

######################################################################################
## Plot proportion forest species of each taxa
## Small land area cells (<100km) have been removed to aid visualisation
## But no minimum number of species used
######################################################################################

figure1_data = prop_forest_df %>% filter(!is.na(n_spec) & n_spec >= 10 & !is.na(land) & land > 0.1) %>%
  mutate(prop_forest = ifelse(is.na(prop_forest), 0, prop_forest))

figure1_data = dplyr::select(figure1_data, x, y, prop_forest, taxa) %>% group_split(taxa) %>% 
  map(., ~rast(dplyr::select(., x,y,prop_forest), crs = behr) %>% 
  project(template_raster)) %>% rast()

names(figure1_data) = c("amphibians", "birds", "mammals", "reptiles")

ggplot() + geom_spatraster(data = figure1_data) + facet_wrap(.~lyr, ncol = 2) + 
  geom_sf(data = coast_outline) +
  theme_classic() + scale_fill_viridis_c(na.value = "white") + 
  labs(fill = "Proportion Forest Species")

## Remove rows that have no species and no land area
## Also use a cut off for what cells are included based on total species from which proportions are calculated
## and the total land area
## In this case a minimum of 10 total species, and 100km2 land area
prop_forest_df = prop_forest_df %>% filter(!is.na(n_spec) & n_spec >= 10 & !is.na(land) & land > 0.1) %>%
  mutate(prop_forest = ifelse(is.na(prop_forest), 0, prop_forest), 
  prop_forest_area = ifelse(is.na(prop_forest_area), 0, prop_forest_area), dist_equator_1000km = abs(y/1000000),
  disturbances = as.factor(disturbances), geological_forest = case_when(geological_forest == 5 ~ 55, geological_forest == 4 ~ 33, geological_forest == 3 ~ 23,
                                                                        geological_forest == 2 ~ 5, geological_forest == 1 ~ 0))

## Scale all continuous variables so their effect sizes are comparable
prop_forest_df = prop_forest_df %>% mutate(scaled_prop_forest_area = scale(prop_forest_area),
                                            scaled_historical_forest_loss = scale(historical_forest_loss),
                                            scaled_dist_equator_1000km = scale(dist_equator_1000km),
                                            scaled_geological_forest = scale(geological_forest))

## In line with Doug Bates thoughts, we will go for the most complex random effect structure which does not create convergence or singlarity issues
## I have already run the model with all variables as random slopes, and this gave singular warnings
## So let's try and just include the most important variables as random slopes
## The most complex I can get without changing the glmmTMB control is including prop_forest_area and historical_forest_loss as additive random slopes
dat = prop_forest_df %>% filter(taxa == "amphibians")
dat$pos <- numFactor(scale(dat$x), scale(dat$y))
# then create a dummy group factor to be used as a random term
dat$ID <- factor(rep(1, nrow(dat)))

model_runtime = function(sample) {

  dat = prop_forest_df %>% filter(taxa == "amphibians") %>% slice_sample(n = sample)
  dat$pos <- numFactor(scale(dat$x), scale(dat$y))
  # then create a dummy group factor to be used as a random term
  dat$ID <- factor(rep(1, nrow(dat)))

  func = function() {
    full_mod = glmmTMB(prop_forest ~ scaled_prop_forest_area*scaled_historical_forest_loss + disturbances + scaled_dist_equator_1000km + scaled_geological_forest + mat(pos+0|ID), 
                    data = dat, 
                    family = "binomial", 
                    weights = n_spec)
  }

  start_time = Sys.time()
  mod = tryCatch(func(), error = function(e) NULL)
  end_time = Sys.time()

  time = end_time - start_time
  out = list(data = dat, model = mod, run = data.frame(sample, time))
  return(out)
}

timings1 = mclapply(seq(100, 600, by = 100), model_runtime, mc.cores = 6)
time_mod = lm(time ~ sample, data = timings)
predict(time_mod, data.frame(sample = nrow(dat)))

testResiduals(a)
simulateResiduals(a, plot = T)
testZeroInflation(a)
testSpatialAutocorrelation(a, x = a$frame$x, y = a$frame$y)
summary(full_mod)

# y_partial = remef(full_mod, fix = c("prop_forest_area", "disturbances", "dist_equator_1000km"), ran = "all")
# test = prop_forest_df %>% add_column(partial_historic_res = y_partial)

## Interaction between forest area and historia forest loss: 
## If you used to have lots of forest then the effect of losing forest is greater
## E.g if you have 0.25 forest now, and you used to have 0.5 (+0.25) then the difference will be smaller than
## if you have 0.5 forest now and you used to have 0.75 (+0.25)
## Interaction between forest area and disturbance:
## Disturbance has a greater effect at low levels of forest area, and little effect at high levels
## Interaction between forest area and distance from equator:
## Forest area has a greater effect far from the equator as opposed to close to the equator when it has a smaller effect

# ggplot(test) + geom_raster(aes(x, y, fill = partial_historic_res)) + theme_classic() +
#   scale_fill_continuous_diverging(mid = 0, palette = "Purple-Green", guide = "colourbar", na.value = "black") + 
#   labs(x = NULL, y = NULL) +
#   theme(axis.title = element_blank(),axis.text = element_blank(),axis.ticks.x = element_blank(),
#         axis.line=element_blank()) + facet_rep_wrap(. ~ taxa, ncol = 2)

plot(ggpredict(full_mod, terms = c("prop_forest_area", "historical_forest_loss[0, 0.2, 0.4]", "disturbances", "dist_equator_1000km"), type = "random")) + theme_classic() +
labs(x = "Proportion of Forested Area", y = "Predicted Proportion of Forest Species", col = "1000 Km \nfrom Equator", title = NULL) + 
theme(text = element_text(size = 20))

model_residuals = prop_forest_df %>% add_column(res = resid(full_mod, type = "response"))

ggplot(model_residuals) + geom_raster(aes(x, y, fill = res)) + theme_classic() +
  scale_fill_continuous_diverging(mid = 0, palette = "Purple-Green", guide = "colourbar", na.value = "black") + 
  labs(x = NULL, y = NULL) + geom_sf(data = coast_outline) +
  theme(axis.title = element_blank(),axis.text = element_blank(),axis.ticks.x = element_blank(),
        axis.line=element_blank()) + facet_rep_wrap(. ~ taxa, ncol = 2)

## Model comparisons
lat_mod = glmer(prop_forest ~ dist_equator_1000km + (dist_equator_1000km|taxa), 
                          data = prop_forest_df, family = "binomial", weights = n_spec)
forest_mod = glmer(prop_forest ~ prop_forest_area + (prop_forest_area|taxa), 
                          data = prop_forest_df, family = "binomial", weights = n_spec)
forest_loss_mod = glmer(prop_forest ~ historical_forest_loss + (historical_forest_loss|taxa), 
                          data = prop_forest_df, family = "binomial", weights = n_spec)
disturb_mod = glmer(prop_forest ~ disturbances + (disturbances|taxa), 
                          data = prop_forest_df, family = "binomial", weights = n_spec)

AIC(lat_mod, forest_mod, forest_loss_mod, disturb_mod)

# Run model accounting for spatial autocorrelation
library(spaMM)
mod_data = prop_forest_df %>% filter(taxa == "amphibians") %>% slice_sample(n = 600)
mod_data = mod_data %>% mutate(n_forest = round(n_spec*prop_forest), n_nonforest = n_spec-n_forest)
spat_forest_amount_mod = fitme(cbind(n_forest, n_nonforest) ~ prop_forest_area*historical_forest_loss + disturbances + dist_equator_1000km + geological_forest + Matern(1|x+y), 
  data = mod_data, family = "binomial")

summary(spat_forest_amount_mod)

testResiduals(spat_forest_amount_mod)
simulateResiduals(spat_forest_amount_mod, plot = T)
testZeroInflation(spat_forest_amount_mod)
testSpatialAutocorrelation(spat_forest_amount_mod, x = mod_data$x, y = mod_data$y)

dd = dist(mod_data[,c("x","y")])
mm = MaternCorr(dd, nu = 1.81, rho = 9.8e-07)

plot(as.numeric(dd), as.numeric(mm), 
  xlab = "Distance between pairs of location [in m]", 
  ylab = "Estimated correlation")