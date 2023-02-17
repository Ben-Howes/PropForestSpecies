####################################################
## Categorize species into matrix or forest
## using bayes theorem
####################################################

library(tidyverse)
library(lme4)
library(brms)
library(janitor)
library(tidybayes)
library(sf)
library(ggeffects)

gpath = "/home/ben/Documents/PhD/matrix_Response/"
setwd(gpath)

spec = read_csv(paste0(gpath, "Data/species_Occurrence.csv")) %>% clean_names() %>% 
    filter(buffer == 200) %>% mutate(matrix_cover = matrix_cover/100)
matrix_Studies = spec %>% group_by(pid) %>% slice_max(matrix_cover, with_ties = F) %>% 
    filter(matrix_cover > 0.9) %>% distinct(pid)
duplicate_Studies = c("PID0100", "PID0101", "PID0102", "PID0103", "PID0104", "PID0105",
                      "PID0106", "PID0107", "PID0108", "PID0109", "PID0110", "PID0111", "PID0113", "PID0114", "PID0117", "PID0119", "PID0124",
                      "PID0138", "PID1002", "PID1005")

spec = spec %>% filter(pid %in% matrix_Studies$pid & !pid %in% duplicate_Studies)

taxa = read_csv(paste0(gpath, "Data/taxa_PID_Lookup.csv")) %>% clean_names()
    
geoms = load(paste0(gpath, "Data/study_Geometries.Rdata"))
study_Centroids = lapply(study_Plot_Points, function(x) 
    data.frame(PID = x$PID[[1]], Latitude = mean(x$Latitude), Longitude = mean(x$Longitude))) %>% bind_rows() %>% 
    clean_names() %>% mutate(small_abs_lat = abs(latitude)/100)

spec = spec %>% left_join(taxa) %>% left_join(study_Centroids)

## Load in avonet data to join with spec data to remove non-terrestrial species
avonet = read.csv(paste0(gpath, "Data/AVONET_BirdLife_Habitat.csv")) %>% clean_names() %>% 
    dplyr::select(species1, habitat) %>% rename(species = species1) %>%
    mutate(species = gsub(" ", "_", tolower(species)))

########################################################################################################
## Run bayesian models, one for each species/study combo 
## predicting occurrence based on matrix cover at 200m
########################################################################################################

x = spec %>% filter(species == unique(spec$species[[1]]))

## Set generic normal dist priors
p1 = set_prior("normal(0,5)", class = "b")
p2 =  set_prior("normal(0,5)", class = "Intercept")

## Make model template
ref_Mod = brm(formula =  occurrence ~ matrix_cover,  
            data = x, 
            family = bernoulli(link = "logit"),
            warmup = 500, 
            iter = 2000, 
            chains = 2, 
            init = "0", 
            cores=2,
            seed = 123,
            prior = c(p1, p2))

## Run models for all species in all studies
loop_Species_Mod = function(x) {

    print(x$species[[1]])

    dat = data.frame(pid = x$pid[[1]], species = x$species[[1]], 
        occurrences = sum(x$occurrence))

    mod = update(ref_Mod, newdata = x)

    out = list(dat, mod)

    return(out)

}

res = spec %>% group_split(pid, species) %>% map(., ~loop_Species_Mod(.))

# save(res, file = paste0(gpath, "Data/category_Models.Rdata"))
load(paste0(gpath, "Data/category_Models.Rdata"))

## Then need to check convergence - Rhat is less than 1.1 (for all parameters)## 2) Other things in coda package
## Next step is to categorise species based on the models that converged correctly
## Forest, Generalist, Matrix

categorise_Species = function(x) {

    dat = x[[1]]
    mod = x[[2]]
    rhats = gelman.diag(as.mcmc(mod))[[1]][,2]

    if(sum(rhats > 1.1) == 0) {
        slope = fixef(mod)[2,]
        low = slope[[3]] 
        up = slope[[4]]
        sig = if(sign(low) == -1 & sign(low) == sign(up)) {
            type = "forest"
        } else if(sign(low) == 1 & sign(low) == sign(up)) {
            type = "matrix"
        } else {
            type = "generalist"
        }
    } else {
        type = "uncertain"
    }

    out = data.frame(pid = dat$pid, species = dat$species, occurrences = dat$occurrences, 
        slope = slope[[1]], low = low, up = up)

}

cat = res %>% map(., ~categorise_Species(.)) %>% bind_rows()
save(cat, file = paste0(gpath, "Data/categorised_Species.Rdata"))
load(paste0(gpath, "Data/categorised_Species.Rdata"))
cat = cat %>% mutate(matrix = ifelse(slope > 0, 1, 0)) %>% left_join(study_Centroids) %>% 
        left_join(taxa)

#################################################################################
## Run mixed effect models to test for change in effect of matrix cover
## with latitude
#################################################################################

## Frequentist approach with different random effects

fMod = glmer(occurrence ~ matrix_cover + matrix_cover:small_abs_lat + (matrix_cover|species) + (1|pid), data = spec, family = "binomial")
fMod1 = glmer(occurrence ~ matrix_cover + matrix_cover:small_abs_lat + (matrix_cover|pid/species), data = spec, family = "binomial")
fMod2 = glmer(occurrence ~ matrix_cover + matrix_cover:small_abs_lat + (matrix_cover|taxa/pid/species), data = spec, family = "binomial")

predicted_Incidence_Plot = plot(ggemmeans(fMod2, terms = c("matrix_cover", "small_abs_lat[0, 0.23, 0.46]"))) + 
    labs(x = "Proportion Matrix Cover 200m", y = "Probability of Incidence", col = "Latitude", title = NULL) + 
        theme_classic() + theme(text = element_text(size = 20))

# What can be inferred from this?
## The effect of matrix on incidence changes with latitude, becoming more positive at higher latitudes
## I think this could mean more matrix species, but it could also be stronger positive effects on the same number of matrix species?

## Compare species with positive sloeps in above random models
## to those from individual bayesian models

fixed = fixef(fMod1)
random = coef(fMod1)[[1]] %>% mutate(pid = str_split(rownames(.), ":", simplify = T)[,2], 
    species = str_split(rownames(.), ":", simplify = T)[,1]) %>% 
        left_join(study_Centroids) %>% mutate(small_abs_lat = abs(latitude)/100) %>% 
            mutate(slope = matrix_cover + (.[,3]*small_abs_lat)) %>% 
                mutate(matrix = ifelse(slope > 0, 1, 0))

## Compare to individual models
ind_Rand_Comp = random %>% left_join(select(cat, species, slope, pid), by = c("species", "pid"), suffix = c("_rand", "_ind")) %>%
    mutate(match = ifelse(sign(slope_rand) == sign(slope_ind), 1, 0))

summary(as.factor(ind_Rand_Comp$match))
# 4788 match of 6098 = 4788/6098 = ~78%

##########################################################################
## Models looking at how change in proportion of matrix species varies 
## with latitude
## This is a definitive increase/decrease in matrix species
## But does require us to define what a matrix species is - propegate error?
################################################################################

prop_Matrix = cat %>% group_by(pid) %>% summarise(n_Matrix = sum(matrix), prop_Matrix = n_Matrix/n(), 
    avg_Slope = median(slope), n_Species = n()) %>% 
        left_join(study_Centroids) %>% left_join(taxa)

prop_Matrix_Mod = glmer(prop_Matrix ~ abs_lat + (abs_lat|taxa), weights = n_Species, data = prop_Matrix, family = "binomial") ## Sig pos

plot(ggpredict(prop_Matrix_Mod, terms = c("abs_lat"))) + 
    scale_colour_viridis_d() + labs(x = "Absolute Latitude", 
                                    y = "Proportion Matrix Species", title = NULL) +
    theme_classic() + theme(text = element_text(size = 20)) + 
    geom_point(data = prop_Matrix, aes(x = abs_lat, y = prop_Matrix), 
               size = 3, inherit.aes = F)

a = ggpredict(prop_Matrix_Mod, terms = c("abs_lat"))

plot(ggpredict(prop_Matrix_Mod, terms = c("abs_lat", "taxa"), type = "random"), ci = F) + 
    scale_colour_viridis_d() + labs(x = "Absolute Latitude", 
                                    y = "Proportion Matrix Species", col = "Taxa", title = NULL) +
        theme_classic() + theme(text = element_text(size = 20)) + 
            geom_point(data = prop_Matrix, aes(x = abs_lat, y = prop_Matrix, col = taxa), 
                       size = 3, inherit.aes = F)

prop_Matrix_Bird_Mod = glm(prop_Matrix ~ abs_lat, weights = n_Species, data = filter(prop_Matrix, taxa == "Birds"), family = "binomial") ## Sig pos
prop_Matrix_Arthropod_Mod = glm(prop_Matrix ~ abs_lat, weights = n_Species, data = filter(prop_Matrix, taxa == "Arthropods"), family = "binomial") ## Sig pos
prop_Matrix_Herp_Mod = glm(prop_Matrix ~ abs_lat, weights = n_Species, data = filter(prop_Matrix, taxa == "Herps"), family = "binomial") ## Sig pos
prop_Matrix_Mammal_Mod = glm(prop_Matrix ~ abs_lat, weights = n_Species, data = filter(prop_Matrix, taxa == "Mammals"), family = "binomial") # Sig pos

####################################################################################################################################################
## Another way to analyse - based on Matt's draft of nature paper
## Predict probability species will be defined as matrix vs non-matrix (1 vs 0) based on abs latitude
## Same issue with defining what a matrix species is a propegating error
## use logistic reression and study as random effect
####################################################################################################################################################

altMod = glmer(matrix ~ small_abs_lat + (small_abs_lat|taxa/pid), data = cat, family = "binomial")
## Strong increase in probability of a species being classified as matrix in higher latitudes

plot(ggpredict(altMod, terms = c("small_abs_lat"))) + 
    scale_colour_viridis_d() + labs(x = "Absolute Latitude", y = "Proportion Matrix Species", title = NULL) +
    theme_classic() + theme(text = element_text(size = 20))

altBirdMod = glmer(matrix ~ small_abs_lat + (1|pid), data = filter(cat, taxa == "Birds"), family = "binomial") ## Sig pos
altArthropodMod = glmer(matrix ~ small_abs_lat + (1|pid), data = filter(cat, taxa == "Arthropods"), family = "binomial") ## Non-sig pos
altHerpMod = glmer(matrix ~ small_abs_lat + (1|pid), data = filter(cat, taxa == "Herps"), family = "binomial") # Just sig pos
altMammalMod = glmer(matrix ~ small_abs_lat + (1|pid), data = filter(cat, taxa == "Mammals"), family = "binomial") ## Sig Pos

## Save all models including individual models
save(fMod, fMod1, fMod2, prop_Matrix, prop_Matrix_Mod,prop_Matrix_Bird_Mod, prop_Matrix_Arthropod_Mod, 
    prop_Matrix_Herp_Mod, prop_Matrix_Mammal_Mod, altMod, altBirdMod, altArthropodMod, altHerpMod, altMammalMod,
        file = paste0(gpath, "Data/freq_Models_and_Results.Rdata"))

## Load
load(paste0(gpath, "Data/freq_Models_and_Results.Rdata"))

###################
## Bayesian
###################

bayesMod = brm(occurrence ~ matrix_cover + matrix_cover:small_abs_lat + (matrix_cover|taxa/pid/species), data = spec, family = "bernoulli",
    chains = 2, cores = 2, iter = 50)

bayesMammal = brm(occurrence ~ matrix_cover + matrix_cover:small_abs_lat + (matrix_cover|pid/species), data = filter(spec, taxa == "Mammals"), 
    family = "bernoulli", chains = 2, cores = 2)
bayesBirds = brm(occurrence ~ matrix_cover + matrix_cover:small_abs_lat + (matrix_cover|pid/species), data = filter(spec, taxa == "Birds"), 
    family = "bernoulli", chains = 2, cores = 2, iter = 50)
bayesHerps = brm(occurrence ~ matrix_cover + matrix_cover:small_abs_lat + (matrix_cover|pid/species), data = filter(spec, taxa == "Herps"), 
    family = "bernoulli", chains = 2, cores = 2, iter = 200)
bayesArths = brm(occurrence ~ matrix_cover + matrix_cover:small_abs_lat + (matrix_cover|pid/species), data = filter(spec, taxa == "Arthropods"), 
    family = "bernoulli", chains = 2, cores = 2, iter = 50)

# load(paste0(gpath, "Data/species_Bayes_Models.Rdata"))
# save(bayesMod, bayesMammal, bayesBirds, bayesHerps, bayesArths, file = paste0(gpath, "Data/bayes_Models.Rdata"))

###########################################################
## Test effect of fragmentation (matrix % at 2000m)
###########################################################

frag_Dat = read_csv(paste0(gpath, "Data/species_Occurrence.csv")) %>% clean_names() %>% mutate(matrix_cover2000 = matrix_cover/100) %>% 
    filter(buffer == 2000) %>% filter(!pid %in% duplicate_Studies) %>% select(-matrix_cover, -abundance, -data, -buffer)

frag_Dat = spec %>% left_join(frag_Dat, by = c("pid", "plot", "species", "occurrence")) %>% mutate(fragmentation = 1 - matrix_cover2000)

#### Find species sensitive to fragmentation

x = frag_Dat %>% group_split(pid, species)

## Set generic normal dist priors
p1 = set_prior("normal(0,5)", class = "b")
p2 =  set_prior("normal(0,5)", class = "Intercept")

## Make model template
ref_Mod1 = brm(formula =  occurrence ~ matrix_cover + fragmentation,  
            data = x[[1]], 
            family = binomial(link = "logit"),
            warmup = 500, 
            iter = 2000, 
            chains = 2, 
            init = "0", 
            cores=2,
            seed = 123,
            prior = c(p1, p2))

## Individual models that categorise species as matrix and also the effect of the fragmentation of the matrix
## Then can use this output to find the proportion of fragmentation sensitive matrix species by latitude

loop_Frag_Species_Mod = function(x) {

    print(x$species[[1]])

    dat = data.frame(pid = x$pid[[1]], species = x$species[[1]], 
        occurrences = sum(x$occurrence))

    mod = update(ref_Mod1, newdata = x)

    out = list(dat, mod)

    return(out)

}

fragmented_Models = x %>% map(., ~loop_Frag_Species_Mod(.))

save(fragmented_Models, file = paste0(gpath, "Data/fragmented_Species_Models.Rdata"))

categorise_Frag_Species = function(x) {

    dat = x[[1]]
    mod = x[[2]]
    rhats = gelman.diag(as.mcmc(mod))[[1]][,2]

    if(sum(rhats > 1.1) == 0) {
        matrix_slope = fixef(mod)[2,1]
        frag_slope = fixef(mod)[3,1]
        if(sign(matrix_slope) > 0) {
            type = "matrix"
        } else if(sign(matrix_slope) <= 0) {
            type = "forest"
        }
        if(sign(frag_slope) < 0) {
            frag = "sensitive"
        } else if(sign(frag_slope) >= 0) {
            frag = "non-sensitive"
        }
    
    } else {
        type = "uncertain"
        frag = "uncertain"
    }

    out = data.frame(pid = dat$pid, species = dat$species, occurrences = dat$occurrences, 
        type = type, frag_sensitive = frag, matrix_slope = matrix_slope[[1]], frag_slope = frag_slope[[1]],
        frag_slope_low = fixef(mod)[3,3], frag_slope_up = fixef(mod)[3,4])

    return(out)

}

fragCat = fragmented_Models %>% map(., ~categorise_Frag_Species(.)) %>% bind_rows() %>% mutate(sen = ifelse(frag_sensitive == "sensitive", 1, 0)) %>% 
    left_join(taxa) %>% left_join(study_Centroids)

save(fragCat, file = paste0(gpath, "Data/fragCat.Rdata"))
load(paste0(gpath, "Data/fragCat.Rdata"))

fragCat = fragCat %>% left_join(dplyr::select(cat, pid, species, matrix, low, up)) %>% 
    mutate(sig_matrix = ifelse(low > 0 & up > 0, 1, 0)) %>% 
    filter(sig_matrix == 1)

fragSenMod = glmer(sen ~ small_abs_lat + (small_abs_lat|taxa/pid), data = fragCat, family = "binomial")

prop_Sen = fragCat %>% group_by(pid) %>%
    summarise(prop_Sen = sum(sen)/n(), num_Sen = sum(sen), num_Spec= n()) %>% left_join(taxa) %>% 
        left_join(study_Centroids)

ggplot(prop_Sen, aes(small_abs_lat, prop_Sen, col = taxa)) + geom_point() + theme_classic() +
        geom_smooth(method = "lm") + facet_wrap(. ~ taxa)

prop_Frag_Sen_Mod = glmer(prop_Sen ~ small_abs_lat + (small_abs_lat|taxa), data = prop_Sen, family = "binomial", weights = num_Spec)

plot(ggemmeans(prop_Frag_Sen_Mod, terms = c("small_abs_lat"))) +
        labs(x = "Absolute Latitude", y = "Proportion of Species \nSensitive to Matrix Fragmentation",
                title = NULL) + theme_classic() + theme(text = element_text(size = 20)) +
                    geom_point(data = prop_Sen, aes(small_abs_lat, y = prop_Sen), size = 2)

plot(ggpredict(prop_Frag_Sen_Mod, terms = c("small_abs_lat", "taxa"), type = "random"), ci = F) +
    labs(x = "Absolute Latitude", y = "Probability of Matrix \n Fragmentation Sensitive",
         title = NULL) + theme_classic() + theme(text = element_text(size = 20)) +
    geom_point(data = prop_Sen, aes(small_abs_lat, y = prop_Sen, col = taxa), 
               size = 2, inherit.aes = F)

## Instead run a bayesian approach to see if that helps matters
prop_Sen_Bayes_Mod = brm(num_Sen | trials(num_Spec) ~ small_abs_lat + (small_abs_lat|taxa), 
                         data = prop_Sen, family = "binomial", init = "0", iter = 10000, chains = 2, cores = 2)

## So from a relatively basic level as you increase latitude the number of matrix fragmentation
## Sensitive species increases
## Interestingly if I run this analysis with only species that respond positively to local matrix cover then
## I do not see this pattern
## So the pattern is driven by forest species, and what I find is that as we increase the latitude
## the proportion of forest species that are negatively affect by fragmentation of the matrix (or loss of frgmentation of forest)
## increases
## I think this backs up what Matt found, basically we have high proportions of forest species
## that like the edge of patches at high latitudes and more core forest species at low lattiudes
## But it seems from this basic analysis that this is not the case for matrix species - do error propegation methods though to
## check this
