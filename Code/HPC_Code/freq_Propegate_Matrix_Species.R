library(tidyverse)
library(janitor)
library(lme4)
library(coda)
library(tidybayes)
library(brms)

gpath = "/rds/general/user/bh719/home/matrix_Models/"
setwd(gpath)

## Load latitude data
geoms = load(paste0(gpath, "Data/study_Geometries.Rdata"))
study_centroids = lapply(study_Plot_Points, function(x) 
  data.frame(PID = x$PID[[1]], Latitude = mean(x$Latitude), Longitude = mean(x$Longitude))) %>% bind_rows() %>% 
  clean_names() %>% mutate(small_abs_lat = abs(latitude)/100)

## Load taxa data
taxa = read.csv(paste0(gpath, "Data/taxa_PID_Lookup.csv")) %>% clean_names()

load(paste0(gpath, "Data/category_Models.Rdata"))

draw_sample = function(x) {
  
  dat = x[[1]]
  mod = x[[2]]
  
  species = dat$species
  pid = dat$pid
  
  draws = spread_draws(mod, b_matrix_cover, ndraws = 1000) %>% 
    mutate(matrix = ifelse(b_matrix_cover > 0, 1, 0)) %>% 
    select(b_matrix_cover, matrix) %>% add_column(.before = "b_matrix_cover", pid, species)
  
  return(draws)
  
}

run_model = function(dat) {
  
  dat = lapply(dat, draw_sample) %>% bind_rows() %>% 
    left_join(taxa) %>% left_join(study_centroids)
  mod = glmer(matrix ~ small_abs_lat + (small_abs_lat|taxa/pid/species), data = dat,
          family = "binomial")
  
  return(mod)
  
}

model = run_model(res)

save(model, file = paste0(gpath, "Results/freq_Prop_Model.Rdata"))
