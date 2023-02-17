##########################################################################################
## Propegate error in matrix species model
## Sample from parameter space for each speceis, categorising as matrix/non-matrix
## Based on positive/negative effect of increasing matrix cover
## Run 1000 models sampling from each species each time
## Can then generate 95% CI from these results
## Run on imperial HPC
##########################################################################################

library(tidyverse)
library(lme4)
library(janitor)
library(sf)
library(tidybayes)
library(ggeffects)

gpath = "/home/ben/Documents/PhD/matrix_Response/"
setwd(gpath)

## Load models
mod_paths = list.files("Results/boot_Mods", pattern = "model_*", full.names = T)

mod_list = list()
prop_list = list()

for(i in 1:length(mod_paths)) {
  print(i)
  load(mod_paths[[i]])
  fix = fixef(result)
  mod_list[[i]] = data.frame(intercept = fix[[1]], small_abs_lat = fix[[2]],
                            sig = ifelse(coef(summary(result))[2,4] < 0.05, 1, 0))
  prop_list[[i]] = getData(result) %>% group_by(pid) %>% 
      summarise(matrix_spec = sum(matrix), total_spec = n(), prop_matrix = matrix_spec/total_spec) %>% 
        mutate(run = i)
  rm(result)
}

mod_param = mod_list %>% bind_rows()
prop_matrix = prop_list %>% bind_rows()

## Plot parameter densities
ggplot(mod_param, aes(small_abs_lat)) + stat_halfeye(outline_bars = TRUE) + 
    theme_classic() + labs(x = "Effect of Latitude on \n Probability of Matrix Species", 
                           y = "Density") + 
    xlim(-0.5, 3) + theme(text = element_text(size = 20)) +
      geom_vline(xintercept = 0, linetype = 2)

summary(as.factor(mod_param$sig))
## 956 of 1000 models are significant positive
## so 95.6%, and all positive

## Plot individual regressions and median line
min_lat = 0
max_lat = 0.5 ## (50 degrees)
x_plot = seq(min_lat, max_lat, 0.01)
y_plot = sapply(1:nrow(mod_param), function(x)
  plogis(mod_param[x,1] + mod_param[x,2] * x_plot) %>% 
    as_tibble() %>% add_column(mod = x), simplify = F) %>% bind_rows()
plot_data = data.frame(x_plot, y_plot) %>% rename("y_plot" = "value")

avg_line = plot_data %>% group_by(x_plot) %>% summarise(y_plot = median(y_plot))

ggplot(plot_data) + geom_line(aes(x_plot, y_plot, group = mod), alpha = 0.01) + 
  labs(x = "Absolute Latitude", y = "Probability of Matrix Species") +
  scale_y_continuous(limits = c(0, 1)) + theme_classic() +
    theme(text = element_text(size = 20)) +
      geom_line(data = avg_line, aes(x_plot, y_plot), linetype = 2, size = 1)

##########################################################################################
## Frequentist model with 1000 samples drawn and species as random effect
##########################################################################################

freq_mod = load(paste0(gpath, "Results/freq_Prop_Model.Rdata"))

plot(ggpredict(model)) + theme_classic() + 
  labs(x = "Absolute Latitude", y = "Probability of Matrix Species") +
      theme(text = element_text(size = 20)) 

##########################################################################################
## Proportion model and figure
##########################################################################################

geoms = load(paste0(gpath, "Data/study_Geometries.Rdata"))
study_centroids = lapply(study_Plot_Points, function(x) 
  data.frame(PID = x$PID[[1]], Latitude = mean(x$Latitude), Longitude = mean(x$Longitude))) %>% bind_rows() %>% 
  clean_names() %>% mutate(small_abs_lat = abs(latitude)/100)

prop_matrix_sd = prop_matrix %>% group_by(pid) %>% 
  summarise(prop_matrix_mean = mean(prop_matrix), sd = sd(prop_matrix)) %>% 
    left_join(study_centroids)

ggplot(prop_matrix_sd, aes(small_abs_lat, prop_matrix_mean)) + 
  geom_pointrange(aes(ymin = prop_matrix_mean-sd, ymax = prop_matrix_mean+sd)) + 
    theme_classic()
