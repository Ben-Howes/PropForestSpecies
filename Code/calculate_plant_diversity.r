#######################################################################
## Calculate plant diversity ebtween forest and non-forest at each latitude
## across teh lonigutde regions
#######################################################################

library(tidyverse)
library(raster)

forest_div = raster(paste0(gpath, "Data/Ranges/PLANT_DIVERSITY/w3_tile_sr1000_for.tif"))
nonforest_div = raster(paste0(gpath, "Data/Ranges/PLANT_DIVERSITY/w3_tile_sr1000_nonfor.tif"))

lat_values = data.frame(lat_min = rep(seq(round(ymin(forest_div)), round(ymax(forest_div)-1), 1), each = 3), 
                        lat_max = rep(seq(round(ymin(forest_div)+1), round(ymax(forest_div)), 1), each = 3))

long_values = data.frame(long_min = c(-180, -30, 60), long_max = c(-30, 60, 180))

band_values = cbind(lat_values, long_values)

calculate_plant_diversity = function(lat_min, lat_max, long_min, long_max) {
  band = raster(xmn = long_min, xmx = long_max, ymn = lat_min, ymx = lat_max, crs = crs(forest_div), resolution = res(forest_div))
  band_cropped_for = crop(forest_div, band)
  band_cropped_nonfor = crop(nonforest_div, band)
  for_rich = cellStats(band_cropped_for, "max", na.rm = T)
  for_mean = cellStats(band_cropped_for, "mean", na.rm = T)
  nonfor_rich = cellStats(band_cropped_nonfor, "max", na.rm = T)
  nonfor_mean = cellStats(band_cropped_nonfor, "mean", na.rm = T)
  region = case_when(long_min == -180 ~ "Americas", long_min == -30 ~ "Europe_Africa_ME", long_min == 60 ~ "Asia_Australia")
  out = data.frame(lat_min, lat_max, for_rich, for_mean, nonfor_rich, nonfor_mean, region)
  return(out)
}

plant_diversity = mapply(calculate_plant_diversity, band_values$lat_min, band_values$lat_max, band_values$long_min, 
                     band_values$long_max, SIMPLIFY = F) %>% bind_rows()

## Latitudes/Longitues with no cells of forest/non-forest will show either Inf or NaN depending on
## whether it was max or mean, so change those to 0
plant_diversity = plant_diversity %>% replace(is.na(.), 0) %>% mutate(for_rich = ifelse(is.infinite(for_rich), 0, for_rich),
                                                                       nonfor_rich = ifelse(is.infinite(nonfor_rich), 0, nonfor_rich))

write.csv(plant_diveristy, paste0(gpath, "Data/Ranges/plant_diversity.csv"), row.names = F)