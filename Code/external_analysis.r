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
library(raster)f
library(maptools)
library(sf)
library(colorspace)

gpath = "/home/ben/Documents/PhD/matrix_Response/"
setwd(gpath)

## Load all latitude and habitat data
mammal_hab = read.csv(paste0(gpath, "Data/Ranges/mammal_habitats.csv"))
mammal_lat = read.csv(paste0(gpath, "Data/Ranges/mammals_quick_ranges.csv"))
amphib_hab = read.csv(paste0(gpath, "Data/Ranges/amphibian_habitats.csv"))
amphib_lat = read.csv(paste0(gpath, "Data/Ranges/amphibians_quick_ranges.csv"))
reptile_hab = read.csv(paste0(gpath, "Data/Ranges/reptile_habitats.csv"))
reptile_lat = read.csv(paste0(gpath, "Data/Ranges/reptiles_quick_ranges.csv"))

## Remove species that live in marine environments. even if they also live terrestrially
## This removes all things like seals as well as obvious ones like whales
## Also remove ranges where the species is not definitely extant
## or where the species is just a migrant or vagrant
species_lat = bind_rows(mammal_lat, amphib_lat, reptile_lat) %>% filter(marine != "true" & category == "Extant (resident)")

## Categorise species as forest or matrix
## Forest species are categorised as any species with a SUITABLE 
## habitat containing the word FOREST
## We also want to remove all marine species (this dataset includes Whales etc)
## so we will remove anything with a SUITABLE habitat that contains the word MARINE

habitat = bind_rows(mammal_hab, amphib_hab, reptile_hab) %>% group_by(taxa, species) %>%
  filter(suitability == "Suitable") %>% 
  summarise(forest = ifelse(sum(grep("forest", habitat, ignore.case = T)) > 0 , 1, 0))

habitat = bind_rows(mammal_hab, amphib_hab, reptile_hab) %>% group_by(taxa, species) %>%
  filter(suitability == "Suitable") %>% 
  summarise(forest_habitats = length(grep("forest", habitat, ignore.case = T) %>% 
                                       grep("degraded", ., ignore.case = T, invert = T)), 
            nonforest_habitats = length(grep("forest", habitat, ignore.case = T, invert = T)), 
            total_habitats = n()) %>% 
  mutate(forest = ifelse(forest_habitats == total_habitats, 1, 
                         ifelse(nonforest_habitats == total_habitats, 0, "generalist")))

## Join the latitude and habitat data together
## and remove species with NA matrix/habitat values as these are either those that have no habitat data (minority)
## or those that did not fulfill our criteria (ma rine species)
## We lost 530 species of mammals that are either marine or with unknown habitat types
## We lose 588 amphibian species which are either marine or with unknown habitat types
hab_lat = species_lat %>% left_join(habitat) %>% filter(!is.na(forest))

## Bird data is in a different format where habitat and max/min lat are all in one csv file
## Remove non-terrestrial species, rename species, and select only relevant columns

bird_data =  read.csv(paste0(gpath, "Data/AVONET_BirdLife_Habitat.csv")) %>% clean_names() %>% 
dplyr::select(species1, habitat, min_latitude, max_latitude, centroid_longitude) %>% 
mutate(species = gsub(" ", "_", tolower(species1)), .before = species1) %>% 
dplyr::select(-species1) %>% filter(!habitat %in% c("Coastal", "Marine") &
                                      !is.na(habitat))

## Bin into latitudes of 1, where a species is present as long as its range crosses the latitude
## e.g a species with min range -5 and max -2 would be shown at -5,-4,-3 and -2
make_lat_rows = function(x) {
  if(!is.na(x$min_latitude)) {
    min = floor(x$min_latitude)
    max = ceiling(x$max_latitude)
    
    min_latitudes = seq(min, max-1, by = 1)
    max_latitudes = seq(min+1, max, by = 1)
    
    out = data.frame(species = x$species, habitat = x$habitat, min_lat = min_latitudes,
                     max_lat = max_latitudes, min_long = x$centroid_longitude)
  }
}

bird_lat = bird_data %>% group_split(species) %>% lapply(., make_lat_rows) %>% bind_rows() %>%
  mutate(forest = ifelse(habitat == "Forest", 1, 0), taxa = "birds") %>% dplyr::select(-habitat) %>%
  mutate(forest = as.character(forest))

## Join the bird data to the other taxa
hab_lat = bind_rows(hab_lat, bird_lat)

## Make new column saying whether a species is in Americas (-180/-30, Europe, Africa and ME (-30/60), or
## Asia and Australia (60/180))
## If species have min and max long we can see if they are in multiple regions and will have to duplicate the rows

find_long_region = function(x) {
  
  min = x$min_long[[1]]
  max = x$max_long[[1]]
  
    if(is.na(max)) {
      region = case_when(min >= -180 & min < -30 ~ "Americas",
                         min>= -30 & min < 60 ~ "Europe_Africa_ME",
                         min >= 60 & min <= 180 ~ "Asia_Australia")
      x = x %>% add_column(region)
      return(x)
    } else {
      min_region = case_when(min >= -180 & min < -30 ~ "Americas",
                         min>= -30 & min < 60 ~ "Europe_Africa_ME",
                         min >= 60 & min <= 180 ~ "Asia_Australia")
      max_region = case_when(max >= -180 & max < -30 ~ "Americas",
                             max >= -30 & max < 60 ~ "Europe_Africa_ME",
                             max >= 60 & max <= 180 ~ "Asia_Australia")
      if(min_region == max_region) {
        x = x %>% add_column(region = min_region)
        return(x)
      } else if(min_region == "Europe_Africa_ME" | max_region == "Europe_Africa_ME") {
        y = x
        x = x %>% add_column(region = min_region)
        y = y %>% add_column(region = max_region)
        out = bind_rows(x,y)
        return(out)
      } else {
        y = x
        z = x
        x = x %>% add_column(region = "Americas")
        y = y %>% add_column(region = "Europe_Africa_ME")
        z = z %>% add_column(region = "Asia_Australia")
        out = bind_rows(x,y,z)
        return(out)
      }
    }
}

hab_lat = hab_lat %>% group_split(species) %>% map(., ~find_long_region(.)) %>% bind_rows()
hab_lat = hab_lat %>% mutate(region = factor(region, levels=c('Americas','Europe_Africa_ME','Asia_Australia')))

## Plot of regions
data("wrld_simpl")
map = st_as_sf(wrld_simpl)

ggplot() + geom_sf(data = map) + theme_classic() + geom_rect(aes(xmin = -180, xmax = -30, ymin = -90, ymax = 90), fill = "red", alpha = 0.2) +
  geom_rect(aes(xmin = -30, xmax = 60, ymin = -90, ymax = 90), fill = "green", alpha = 0.2) +
  geom_rect(aes(xmin = 60, xmax = 180, ymin = -90, ymax = 90), fill = "blue", alpha = 0.2)

## Make data frame with species richness of matrix species and proportion of matrix species at each degree of lat
prop_lat = hab_lat %>% group_by(taxa, min_lat, region) %>% 
  summarise(n_forest = sum(forest == 1), n_spec = n(), prop_forest = n_forest/n_spec)

## Remove the 3 NA bird species - should be fixed with birdliufe data
prop_lat = prop_lat %>% filter(!is.na(region))

ggplot(filter(hab_lat, !is.na(region)), aes(min_lat, col = as.factor(forest))) + geom_density(size = 2) + theme_classic() +
  labs(x = "Latitude", y = "Density", col = "Species Type") +
  scale_colour_viridis_d(end = 0.5, labels = c("Matrix", "Forest")) + theme(text = element_text(size = 20)) +
  geom_vline(xintercept = 0, linetype = 2) + facet_rep_grid((taxa ~ region))

ggplot(filter(prop_lat, !is.na(prop_forest) & n_spec > 10), aes(min_lat, prop_forest)) + geom_point(aes(size = n_spec)) + 
  theme_classic() + labs(x = "Latitude", y = "Proportion Forest Species", title = "Proportion of Inland Forest Species",
                         size = "Number of \nSpecies") +
  scale_size(breaks = c(25, 100, 1000), range = c(0.1, 7.5)) +
  theme(text = element_text(size = 20)) + geom_smooth(method = "loess") + facet_rep_grid(taxa ~ region) +
  scale_x_continuous(breaks = c(-80, -60, -40,-20,0,20,40,60,80), limits = c(-75, 75)) +
  geom_vline(xintercept = 0, linetype = 2)

################################################
## Are these trends explained by forest area?
################################################

## Load forest area data calculated from Hansen 2013
## at 1 degree latitude intervals with 1000m scale
## Forest categorised at 30, 50, and 70% tree cover cut-offs

forest_area = read.csv(paste0(gpath, "Data/Ranges/forest_by_latitude.csv")) %>% 
  mutate(forest_area_km = forest_area/1000000, prop_forest_area = forest_area/land_area) %>%
  dplyr::select(-X) %>% mutate(region = case_when(longitude_min == -180 ~ "Americas",
                                                  longitude_min == -30 ~ "Europe_Africa_ME",
                                                  longitude_min == 60 ~ "Asia_Australia"))

ggplot(forest_area, aes(latitude_min, prop_forest_area, col = as.factor(forest_def))) + geom_line(size = 2) + 
  theme_classic() + labs(x = "Latitude", y = "Proportion Forested Area", 
                         col = "Tree Cover Forest\n Cut-Off Value") +
  theme(text = element_text(size = 20)) + scale_colour_viridis_d() +
  scale_x_continuous(breaks = c(-90, -75, -60, -45, -30,-15,0,15,30,45,60, 75, 90)) +
  facet_rep_wrap(. ~ region)

## Join forest area data to proportional matrix data and plot to visualise relationship
prop_lat_forest = prop_lat %>% left_join(forest_area, by = c("min_lat" = "latitude_min", "region")) %>% 
  mutate(small_lat = min_lat/100)

## Hansen dataset including the land area mask does not go all the way from -90 to 90
## As it does not model areas with 0 tree cover
## As such we remove those latitudes
## Secondly, we are currently creating proportions from extremely variable sample sizes
## for example at very low latitudes there often less than 3 amphibian species
## for now we will use 10 as an arbitrary cut off in minimum number of species to 
## calculate a proportion

prop_lat_forest = prop_lat_forest %>% filter(land_area > 0 & n_spec > 10) %>% 
  mutate(region = factor(region, levels=c('Americas','Europe_Africa_ME','Asia_Australia')))

ggplot(filter(prop_lat_forest, forest_def == 70), aes(min_lat, prop_forest, col = prop_forest_area)) + geom_point(size = 4) +
  theme_classic() + facet_rep_grid((taxa ~ region)) + scale_colour_viridis_c() +
  labs(x = "Latitude", y = "Proportion Forest Species", col = "Proportion Forested\nArea") +
  theme(text = element_text(size = 20)) + ggtitle("Tree Cover Forest Cut-Off = Cols\n Taxa = Rows")

## Plot and model prop_marix ~ proportion forested area

ggplot(filter(prop_lat_forest, forest_def == 70), aes(prop_forest_area, prop_forest, col = taxa)) +
  geom_point() + theme_classic() + scale_colour_viridis_d(end = 0.8) +
  labs(x = "Proportion Forested Area", y = "Proportion Forest Species", col = "Taxa") +
  geom_smooth() + theme(text = element_text(size = 20)) + facet_rep_grid(taxa ~ region)

forest_amount_mod = glmer(prop_forest ~ prop_forest_area + (prop_forest_area|region/taxa), 
                          data = filter(prop_lat_forest, forest_def == 70), family = "binomial", weight = n_spec)

forest_amount_mod %>% summary()
r.squaredGLMM(forest_amount_mod)

plot(ggpredict(forest_amount_mod, terms = c("prop_forest_area", "taxa", "region"), type = "random"), add.data = T) + theme_classic() +
labs(x = "Proportion of Forested Area", y = "Predicted Proportion of Forest Species", col = "Taxa", title = NULL) + 
theme(text = element_text(size = 20))

model_residuals = filter(prop_lat_forest, forest_def == 70) %>% add_column(res = resid(forest_amount_mod, type = "response"))

ggplot(model_residuals, aes(min_lat, res, col = taxa)) + geom_point() + theme_classic() +
scale_colour_viridis_d() + labs(x = "Latitude", y = "Residuals", col = "Taxa") +
theme(text = element_text(size = 20)) + facet_rep_grid(taxa ~ region) +
scale_x_continuous(breaks = c(-90, -75, -60, -45, -30,-15,0,15,30,45,60, 75, 90), limits = c(-90,90)) +
geom_hline(yintercept = 0, linetype = "dashed")

## Model residuals against latitude to determine lattidues with higher matrix specis than expected
## We could also just add latitude to the above model, but this way is more visual
model_residuals = model_residuals %>% mutate(region = as.factor(region), region_int = as.factor(case_when(
  region == "Americas" ~ 1, region == "Europe_Africa_ME" ~ 2, region == "Asia_Australia" ~ 3)))

mammal_gam = gam(res ~ s(min_lat) + s(min_lat, by = region_int), data = filter(model_residuals, taxa == "mammals" & forest_def == 70))
summary(mammal_gam)

## Plot output on map of areas with more and less matrix species than expected
## based on proportion of forested area at that latitude

## Load in simple map of world
data("wrld_simpl")
map = st_as_sf(wrld_simpl)

# Make base raster with one pixel per degree latitude
base = raster(nrows = 141, ncols = 361, xmn = -180, xmx = 180, ymn = -60, ymx = 80,
            vals = rep(seq(80, -60, -1), each = 361))
base1 = setValues(base, values = as.factor(rep(c(rep(1, times = 150),
                                       rep(2, times = 90),
                                       rep(3, times = 121)), times = 141)))
full_base = stack(base, base1)

## Predict models onto this raster
names(full_base) = c("min_lat", "region_int")
mammal_predictions = predict(full_base, mammal_gam)

## Mask new predicted rasters with world map to get world shaped map
mammal_predictions = mask(mammal_predictions, map)

ggplot() + geom_raster(data = as.data.frame(mammal_predictions,xy=T), 
                     aes(x,y,fill=layer)) + theme_classic() +
geom_hline(yintercept = 0, linetype = "dashed") +
ylim(-60, 90) + scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0,
                                                 na.value = "white") +
theme(legend.position = "bottom", axis.line = element_blank(), axis.text.x = element_blank(), 
      axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), 
      axis.title.y = element_blank()) + 
guides(fill = guide_colourbar(
        ticks = FALSE,
        title.vjust = 1.5, 
        frame.linewidth = 1.5,
        barwidth = 25,
        barheight = 1.5,
        frame.colour = "black", 
        ticks.colour = "black",
        ticks.linewidth = 0.3,
        title = "More Matrix Species                   More Forest Species",
        title.position = "top")) + theme(text = element_text(size = 15))

## Could the discrepancies be caused by historic forest loss?
## Load and join historic forest loss to data frame

forest_loss = read_csv(paste0(gpath, "Data/forest_loss.csv"))

model_residuals = model_residuals %>% left_join(forest_loss, by = c("min_lat" = "lat_min", "region")) %>% 
  mutate(land_area_km = land_area/1000000, prop_loss = forest_loss_area_km/land_area_km,
         total_prop_forest = (forest_loss_area_km+forest_area_km)/land_area_km)

## Betts et al did not model forest loss in the boreal regions so this part of the analyssi only looks at latitude between -55 and 60
## degrees of lattiude
## You can see this lack of modelling in these regions in the figure below
ggplot(model_residuals, aes(min_lat, prop_loss, col = as.factor(region))) + geom_line(size = 2) + theme_classic() +
  scale_x_continuous(breaks = c(seq(-90,90, by = 10))) + scale_colour_viridis_d() +
  labs(x = "Latitude",y = "Proportion of Forest Loss from Total Land Area", col = "Region") +
  theme(text = element_text(size = 20))

## This actually looks quite promising, especially for amphibians
ggplot(filter(model_residuals), aes(prop_loss, res, col = min_lat)) + 
  geom_point() + theme_classic() +
  scale_colour_viridis_c(end = 0.8) + labs(x = "Proportion Forest Loss", y = "Residuals", col = "Latitude") +
  theme(text = element_text(size = 20)) + facet_rep_grid(taxa ~ region) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm", se = F)

forest_loss_mod = lmer(res ~ prop_loss + (prop_loss|region/taxa), 
                       data = filter(model_residuals))

summary(forest_loss_mod)
r.squaredGLMM(forest_loss_mod)

plot(ggpredict(forest_loss_mod, terms = c("prop_loss", "taxa", "region"), type = "random"), add.data = T) + 
  theme_classic() + labs(x = "Proportion Forest Loss of Total Land Area", y = "Residuals", col = "Taxa") +
  theme(text = element_text(size = 20)) + ggtitle(NULL)

historic_forest_mod = glmer(prop_forest ~ total_prop_forest + (total_prop_forest|region/taxa), 
             data = filter(model_residuals), family = "binomial",
             weights = n_spec)

summary(historic_forest_mod)
r.squaredGLMM(historic_forest_mod)

current_forest_comp_mod = glmer(prop_forest ~ prop_forest_area + (prop_forest_area|region/taxa), 
                                data = filter(model_residuals), 
                                family = "binomial", weights = n_spec)

summary(current_forest_comp_mod)
r.squaredGLMM(current_forest_comp_mod) ## Marginal R2 of 0.09, conditional of 0.242
r.squaredGLMM(historic_forest_mod) ## Marginal R2 of 0.129, conditional of 0.247
AIC(current_forest_comp_mod, historic_forest_mod)

## Could disturbance explain the differences? Extinction filter hypothesis

disturbed = read.csv(paste0(gpath, "Data/disturbed_area.csv")) %>% dplyr::select(-X)
model_residuals = model_residuals %>% left_join(disturbed, by = c("min_lat" = "lat_min", "region")) %>% 
  mutate(prop_dist = disturbed_area/land_area_km)

ggplot(model_residuals, aes(x = min_lat, prop_dist, col = region)) + geom_line(size = 2) + theme_classic() +
  scale_colour_viridis_d(end = 0.8) + labs(x = "Latitude", y = "Proportion Disturbed Land Area", col = "Region") +
  theme(text = element_text(size = 20))

ggplot(model_residuals, aes(prop_dist, prop_forest, col = min_lat)) + geom_point() + geom_smooth(se = F) +
  facet_rep_grid(taxa ~ region) + theme_classic() + scale_colour_viridis_c(end = 0.8) +
  labs(x = "Proportion Disturbed Land Area", y = "Proportion Forest Species", col = "Latitude") + theme(text = element_text(size = 20))

ggplot(model_residuals, aes(prop_dist, res, col = min_lat)) + geom_point() + geom_smooth(method = "lm", se = F) +
  facet_rep_grid(taxa ~ region) + theme_classic() + scale_colour_viridis_c(end = 0.8) +
  labs(x = "Proportion Disturbed Land Area", y = "Residuals", col = "Latitude") + theme(text = element_text(size = 20)) +
  geom_hline(yintercept = 0, linetype = "dashed")

## It doesn't look like historic natural disturbance explains either the main difference in proportion of forest species
## or the proportion of forest species once we account for total forested areas

###### PLant diversity as a proxy for strucutral complexity across latitudes
## We expect the structural complexity of forests to vary more steeply across latitudes
## than those of nonforests
## WHich may drive the semi-latitudinal gradient seen

plant_diversity = read.csv(paste0(gpath, "Data/Ranges/plant_diversity.csv"))

model_residuals = model_residuals %>% left_join(plant_diversity, by = c("min_lat" = "lat_min", "region")) %>% 
  mutate(prop_plant_rich = for_rich/(for_rich + nonfor_rich), prop_plant_mean = for_mean/(for_mean + nonfor_mean))

ggplot(model_residuals, aes(min_lat, prop_plant_rich)) + geom_point() + theme_classic() + 
  facet_rep_grid(. ~ region) + geom_smooth()
