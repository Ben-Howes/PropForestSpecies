#############################################################
## Analyse the spatmodels run on the HPC
## using glmmTMB to account for spatialautocorrelation
## Models were run separately per taxa
## with 100 models run per taxa
## and run with 10% of data for each of them in each model
#############################################################

library(tidyverse)
library(DHARMa) ## For model testing
library(glmmTMB) ## Model data extraction
library(ggdist) ## For halfeye density plots

## Set path and wd
gpath = "/home/ben/Documents/PhD/PropForestSpecies/"
setwd(gpath)

## Get paths to models
model_paths = list.files(paste0(gpath, "Results/spatResults/Results"), full.names = T)

## Write a function to open each model
## get the relevant information
## then clear it from the system
get_data = function(path) {

    load(path) ## Load model
    mod = out$model ## get model
    taxa = out$taxa ## get taxa name

    ## Get coefficients
    coef = summary(mod)$coef$cond %>%
            as.data.frame() %>%
            rename("est" = 1, "std_err" = 2, "z_val" = 3, "p_val" = 4)

    CI = confint(mod, parm = rownames(coef)) %>% 
            as.data.frame() %>% 
            rename("low_2.5" = 1,"up_97.5" = 2, "est" = 3) %>%
            dplyr::select(-est)

    final = data.frame(predictor = rownames(coef), coef, CI, taxa = taxa, run = out$PBS)

    return(final)
}


test = mclapply(model_paths[1:10], get_data, mc.cores = 8)

