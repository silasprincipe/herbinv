#### Modelling of Herbivorous invertebrates of coral reefs ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2021

### Evaluating results of the Virtual Species testing ###

# Load libraries and set folder ----
library(tidyverse)
setwd("data/vsp/modeling")

# Load data and separate for each VSP ----
metrics <- read.csv("ens_metrics.csv")
head(metrics)

# Separate VSP
vsp1 <- metrics[grep("vsp1", metrics$model),]
vsp1.ln <- vsp1[grep("ln", vsp1$model),]
vsp1.hn <- vsp1[grep("hn", vsp1$model),]
rm(vsp1)

vsp2 <- metrics[grep("vsp2", metrics$model),]

# Evaluate models ----
##### Evaluate VSP 1 ----

vsp1.ln <- vsp1.ln %>%
        separate(model, c("vsp", "strat", "sampn", "patype", "pan", "rep"))


ggplot(vsp1.ln) +
        geom_boxplot(aes(x = patype, y = tss, fill = pan)) +
        facet_wrap(~strat)

ggplot(vsp1.ln) +
        geom_boxplot(aes(x = patype, y = jac, fill = pan)) +
        facet_wrap(~strat)


##### Evaluate VSP1 HN ----
vsp1.hn <- vsp1.hn %>%
        separate(model, c("vsp", "strat", "sampn", "patype", "pan", "rep"))


ggplot(vsp1.hn) +
        geom_boxplot(aes(x = patype, y = tss, fill = pan)) +
        facet_wrap(~strat)

ggplot(vsp1.hn) +
        geom_boxplot(aes(x = patype, y = jac, fill = pan)) +
        facet_wrap(~strat)


##### Evaluate VSP2 ----
vsp2 <- vsp2 %>%
        separate(model, c("vsp", "strat", "patype", "pan", "rep"))


ggplot(vsp2) +
        geom_boxplot(aes(x = patype, y = tss, fill = pan)) +
        facet_wrap(~strat)

ggplot(vsp2) +
        geom_boxplot(aes(x = patype, y = jac, fill = pan)) +
        facet_wrap(~strat)


# Load models and project to future ----

setwd("../../..")
# Load needed functions ----
source("functions/vsp_modeling.r")

# Load enviromental data ----
source("functions/Varload.r")
env <- var.load(folder = "crop_layers/", "env_layers.txt", "2_100")
r12 <- var.load(folder = "proj_layers/ssp126/", "env_layers.txt", "2_100")
r24 <- var.load(folder = "proj_layers/ssp245/", "env_layers.txt", "2_100")
r37 <- var.load(folder = "proj_layers/ssp370/", "env_layers.txt", "2_100")
r85 <- var.load(folder = "proj_layers/ssp585/", "env_layers.txt", "2_100")

proj.rast <- list(env, r12, r24, r37, r85)
names(proj.rast) <- c("current", "r12", "r24", "r37", "r85")


setwd("data/vsp")


# Load original rasters
vsp1 <- raster("vsp1_pa.tif")
vsp2 <- raster("vsp2_pa.tif")

# Future projections
vsp1.fut <- stack(list.files(pattern = "vsp1_r"))
vsp2.fut <- stack(list.files(pattern = "vsp2_r"))

setwd("modeling")

models <- list.files("models")
models <- models[grep("2t", models[grep("mapa", models)])]

ens.jacc <- data.frame(matrix(nrow = length(models), ncol = length(proj.rast)))
names(ens.jacc) <- names(proj.rast)

for (i in 1:length(models)) {
        
        model <- readRDS(paste0("models/", models[i]))
        
        names(model) <- c("gbm", "rf", "glm", "gam",
                          "gbmt", "rft", "glmt", "gamt")
        
        ### Make projections ----
        pred.gbm <- pred.sdm(proj.rast, model$gbm, model$gbmt)
        pred.rf <- pred.sdm(proj.rast, model$rf, model$rft)
        pred.glm <- pred.sdm(proj.rast, model$glm, model$glmt)
        pred.gam <- pred.sdm(proj.rast, model$gam, model$gamt)
        
        par(mfrow = c(2,2))
        for (z in 1:5) {
                plot(pred.gbm[[z]], legend = F,
                     main = names(proj.rast)[z])
                plot(pred.rf[[z]], legend = F, main = models[i])
                plot(pred.glm[[z]])
                plot(pred.gam[[z]])
        }
        
        ### Creates an ensemble ----
        ens.model <- pred.gbm + pred.rf + pred.glm + pred.gam
        ens.model <- ens.model >= 3
        
        ### Get jaccard values for the current and future projections
        if (grepl("vsp1", models[i])) {
                ens.jacc[i,1] <- jaccard.eval(ens.model[[1]], vsp1)
                for (z in 1:4) {
                        ens.jacc[i, (z + 1)] <-
                                jaccard.eval(ens.model[[(z + 1)]], vsp1.fut[[z]])
                }
        } else{
                ens.jacc[i,1] <- jaccard.eval(ens.model[[1]], vsp2)
                for (z in 1:4) {
                        ens.jacc[i, (z + 1)] <-
                                jaccard.eval(ens.model[[(z + 1)]], vsp2.fut[[z]])
                }
        }

}

ens.jacc$models <- models

ens.jacc <- ens.jacc %>% 
        separate(models, c("vsp", "strat", "sampn", "patype", "pan", "rep")) %>%
        select(-vsp, -pan)

ens.jacc.m <- ens.jacc %>%
        group_by(strat, sampn, patype) %>%
        summarise_if(is.numeric, mean)

ens.jacc.m <- ens.jacc %>%
        pivot_longer(1:5,
                     names_to = "scenario", values_to = "jacc")

ggplot(ens.jacc.m) +
        geom_boxplot(aes(x = patype, y = jacc, fill = scenario))+
        facet_wrap(~strat+sampn)+
        theme_bw()

# Plot some to have an overview
