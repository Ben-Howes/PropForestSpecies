##########################################################################################
## Propegate error in fragmented species models
## Sample from parameter space for each species, categorising as matrix/non-matrix
## Based on positive/negative effect of increasing matrix cover
## Then categorise matrix species as sensitive or non-sensitive to
## fragmentaion which is 1 - matrix_cover at 2000m
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
mod_paths = list.files("Results/boot_Frag_Mods", pattern = "model_*", full.names = T)

mod_list = list()
prop_list = list()

for(i in 1:length(mod_paths)) {
  print(i)
  load(mod_paths[[i]])
  fix = fixef(result)
  mod_list[[i]] = data.frame(intercept = fix[[1]], small_abs_lat = fix[[2]],
                             sig = ifelse(coef(summary(result))[2,4] < 0.05, 1, 0))
  prop_list[[i]] = getData(result) %>% group_by(pid) %>% 
    summarise(sen_spec = sum(frag), total_spec = n(), prop_sen = sen_spec/total_spec) %>% 
    mutate(run = i)
  rm(result)
}

mod_param = mod_list %>% bind_rows()
prop_sen = prop_list %>% bind_rows()

mod_param$sig %>% sum() ## 41, so 4.1% of results are significant, all positive

## Plot parameter densities
ggplot(mod_param, aes(small_abs_lat)) + stat_halfeye(outline_bars = TRUE) + 
  theme_classic() + labs(x = "Effect of Latitude on Probability of Matrix Species \nbeing Categorised as Fragmentation Sensitive", 
                         y = "Density") + 
  xlim(-1, 2.5) + theme(text = element_text(size = 20)) +
  geom_vline(xintercept = 0, linetype = 2)


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
  labs(x = "Absolute Latitude", y = "Probability of Matrix Species\n being Categorised as Fragmentation Sensitive") +
  scale_y_continuous(limits = c(0, 1)) + theme_classic() +
  theme(text = element_text(size = 20)) +
  geom_line(data = avg_line, aes(x_plot, y_plot), linetype = 2, size = 1)
