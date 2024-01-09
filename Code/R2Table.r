
library(tidyverse)
library(ggstance) ## For specific positioning in ggplot
library(glmmTMB) ## extracting model output
library(lemon) ## Facet_rep_wrap
library(broom)
library(broom.mixed) ## tidy
library(kableExtra) ## kable_styling 
library(magick)
library(webshot2)

## Set path and wd
gpath = "/home/ben/Documents/PhD/PropForestSpecies/"
setwd(gpath)

load("Results/randLatModels.RData")
load("Results/randFullModels.RData")
load("Results/taxaLatModel.RData")
load("Results/taxaFullModel.RData")

## Combine the individual taxa models with the combined taxa models for plotting
randLatModels = c(randLatModels, list(all = taxaLatModel))
randFullModels = c(randFullModels, list(all = taxaFullModel))

latResults = lapply(names(randLatModels), function(x) {
    mod = randLatModels[[x]]
    rsq = summary(mod)$r.sq
    out = data.frame(taxa = x, rsq = rsq)
})

latResults = bind_rows(latResults) %>%
    mutate(taxa = str_to_title(taxa), rsq = round(rsq, 2)) %>%
    rename("Taxa" = 1, "Latitude-Only R2" = 2) %>%
    arrange(-desc(Taxa))

fullResults = lapply(names(randFullModels), function(x) {
    mod = randFullModels[[x]]
    rsq = summary(mod)$r.sq
    out = data.frame(taxa = x, rsq = rsq)
})

fullResults = bind_rows(fullResults) %>%
    mutate(taxa = str_to_title(taxa), rsq = round(rsq, 2)) %>%
    rename("Taxa" = 1, "Full Model R2" = 2) %>% 
    arrange(-desc(Taxa)) %>%
    left_join(latResults, by = "Taxa") %>%
    relocate(1,3,2)

kable(fullResults, align = c("l", "c", "c")) %>%
kable_classic() %>%
kable_styling(latex_options = c("striped", "scale_down"), full_width = FALSE) %>%
save_kable(file = "Paper/Tables/R2Table.png", zoom = 2)
