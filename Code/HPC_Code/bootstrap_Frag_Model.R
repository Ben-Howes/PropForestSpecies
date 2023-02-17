library(tidyverse)
library(janitor)
library(brms)
library(coda)
library(tidybayes)
library(lme4)

seed = commandArgs(trailingOnly = TRUE)
set.seed(seed)

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
load(paste0(gpath, "Data/fragmented_Species_Models.Rdata"))

draw_sample = function(x, y) {
  
  dat = x[[1]]
  mod = x[[2]]
  species = dat$species
  pid = dat$pid
  
  frag_dat = y[[1]]
  frag_mod = y[[2]]
  frag_species = frag_dat$species
  frag_pid = frag_dat$pid
  
  
  matrix_draws = spread_draws(mod, b_matrix_cover, ndraws = 1) %>% 
    mutate(matrix = ifelse(b_matrix_cover > 0, 1, 0)) %>% 
    select(b_matrix_cover, matrix) %>% add_column(.before = "b_matrix_cover", pid, species)
  
  frag_draws = spread_draws(frag_mod, b_fragmentation, ndraws = 1) %>% 
    mutate(frag = ifelse(b_fragmentation < 0, 1, 0)) %>% 
    select(b_fragmentation, frag) %>% add_column(.before = "b_fragmentation", frag_pid, frag_species) %>% 
    rename(species = frag_species)
  
  out = left_join(matrix_draws, frag_draws, by = c("species"))
  
  return(out)
  
}

run_model = function(dat, dat1) {
  
    dat2 = mapply(draw_sample, dat, dat1, SIMPLIFY = F) %>% bind_rows() %>% 
        left_join(taxa) %>% left_join(study_centroids) %>% 
          filter(matrix == 1)
    rm(dat,dat1) ## Remove so the objects don't get saved with the model!
    mod = glmer(frag ~ small_abs_lat + (small_abs_lat|taxa/pid), data = dat2, 
                 family = "binomial")

    return(mod)
  
}

result = run_model(res, fragmented_Models)

save(result, file = paste0(gpath, "Results/boot_Frag_Mods/model_", seed, ".Rdata"))