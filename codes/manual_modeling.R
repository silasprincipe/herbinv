#### Modelling of Herbivorous invertebrates of coral reefs ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2021

### Modeling - Sea-urchins ###

# Load packages ----
library(dismo)
library(randomForest)

# Set seed for replicability
set.seed(2932)

# Set species
species <- "lyva"

# Load environmental data ----
source("functions/Varload.r")
env <- var.load(folder = "crop_layers/", layers = "env_layers.txt")

r12 <- var.load(folder = "proj_layers/ssp126/", layers = "env_layers.txt")

r24 <- var.load(folder = "proj_layers/ssp245/", layers = "env_layers.txt")

r37 <- var.load(folder = "proj_layers/ssp370/", layers = "env_layers.txt")

r85 <- var.load(folder = "proj_layers/ssp585/", layers = "env_layers.txt")

# Load species data ----
# Load species data
pts <- read.csv(paste("data/", species, "/", species, "_pts.csv", sep = ""))

dsp <- read.csv(paste("data/", species, "/", species, "_dsp.csv", sep = ""))
dsp[is.na(dsp$dsp),] <- 0

# Load split table (block cross validation)
blocks <- read.csv(paste0("data/", species, "/", species, "_foldindex.csv"))
blocks <- blocks[,1]

blocks.lat <- read.csv(paste0("data/", species, "/", species, "_foldindex_lat.csv"))
blocks.lat <- blocks.lat[,1]

data <- raster::extract(env, pts)

data <- cbind(dsp, data)

# normal blocks / unweighted
gbm1 <- gbm.step(data, gbm.x = 2:8, gbm.y = 1, fold.vector = blocks, n.folds = 5)

# normal blocks / weigthed
gbm2 <- gbm.step(data, gbm.x = 2:8, gbm.y = 1, fold.vector = blocks,
                 n.folds = 5, site.weights = wgt)

# lat blocks / unweighted
gbm3 <- gbm.step(data, gbm.x = 2:8, gbm.y = 1, fold.vector = blocks.lat, n.folds = 5)

# lat blocks / weigthed
gbm4 <- gbm.step(data, gbm.x = 2:8, gbm.y = 1, fold.vector = blocks.lat,
                 n.folds = 5, site.weights = wgt)


pred.gbm1 <- predict(env, gbm1)
pred.gbm2 <- predict(env, gbm2)
pred.gbm3 <- predict(env, gbm3)
pred.gbm4 <- predict(env, gbm4)

par(mfrow = c(1,2))

thr1 <- threshold(eval <- dismo::evaluate(data[dsp$dsp == 1,2:8], data[dsp$dsp == 0,2:8], gbm1))
thr2 <- threshold(eval <- dismo::evaluate(data[dsp$dsp == 1,2:8], data[dsp$dsp == 0,2:8], gbm2))
thr3 <- threshold(eval <- dismo::evaluate(data[dsp$dsp == 1,2:8], data[dsp$dsp == 0,2:8], gbm3))
thr4 <- threshold(eval <- dismo::evaluate(data[dsp$dsp == 1,2:8], data[dsp$dsp == 0,2:8], gbm4))

thresholded <- lapply(list(pred.gbm1, pred.gbm2, pred.gbm3, pred.gbm4),
                      function(x){
                              x[x >= thr$spec_sens] <- 1
                              x[x < thr$spec_sens] <- 0
                              x
                      })

lapply(thresholded, plot)



thr <- function(x, thres){
        x[x >= thres] <- 1
        x[x < thres] <- 0
        x
}


pred.gbm1 <- thr(pred.gbm1, thr1$spec_sens)
pred.gbm2 <- thr(pred.gbm2, thr2$spec_sens)
pred.rf1 <- thr(pred.rf1, thr.rf$spec_sens)



# Random forest
library(randomForest)

rf1 <- randomForest(dsp ~ ., data = data, strata = factor(c(0,1)))

pred.rf1 <- predict(env, rf1)

thr.rf <- threshold(eval <- dismo::evaluate(data[dsp$dsp == 1,2:8], data[dsp$dsp == 0,2:8], rf1))

varimp <- var.importance(data, model = gbm1)

varimp <- apply(varimp, 2, mean)


varimp2 <- var.importance(data, model = rf1)

varimp2 <- apply(varimp2, 2, mean)

###
par(mfrow = c(1, 3))
plot(pred.gbm1)
points(pts[dsp == 1,], col = "black", cex = 0.5, pch = 20)
plot(pred.gbm2)
points(pts[dsp == 1,], col = "black", cex = 0.5, pch = 20)
plot(pred.rf1)
points(pts[dsp == 1,], col = "black", cex = 0.5, pch = 20)
points(pts[dsp == 0,], col = "blue", cex = 0.5, pch = 20)



#### GLM
glm.formula <- makeFormula("dsp", data[, 2:8], type = "quadratic")

glm.formula <- dsp ~ 1 + BO21_tempmean_ss + I(BO21_tempmean_ss^2) + BO21_salinitymean_ss + 
        BO21_lightbotmax_bdmean + I(BO21_lightbotmax_bdmean^2) + 
        BO_ph + I(BO_ph^2) + BO21_chlomean_ss + I(BO21_chlomean_ss^2) + 
        BO21_silicatemean_ss + I(BO21_silicatemean_ss^2) + BO21_dissoxmean_ss + 
        I(BO21_dissoxmean_ss^2)

glm1 <- glm(glm.formula, binomial(link = 'logit'),
            data = data, weights = wgt,
            mustart = rep(0.5, nrow(data)),
            control = glm.control(epsilon = 1e-08, maxit = 50, trace = FALSE))


glm.enh <- MASS::stepAIC(glm1, glm.formula, data = data,
                   direction = "both", trace = FALSE, 
                   weights = wgt, 
                   steps = 10000, mustart = rep(0.5, nrow(data)))


pred.glm1 <- predict(env, glm.enh)

thr.glm <- threshold(eval <- dismo::evaluate(data[dsp$dsp == 1,2:8], data[dsp$dsp == 0,2:8], glm1))

pred.glm1 <- thr(pred.glm1, thr.glm$spec_sens)

partial_response(glm.enh, data[,2:8])

par(mfrow = c(2,2))
predplots(current = env, model = gbm1, future = r85, samp = 1000, thresh = thr1$spec_sens)




gbm.current <- predict(env, gbm2)
gbm.fut     <- predict(r24, gbm2)

gbm.current <- thr(gbm.current, thr2$spec_sens)
gbm.fut     <- thr(gbm.fut   , thr2$spec_sens)

par(mfrow = c(2,2))
plot(gbm.current)
plot(bio.lyva.cur)
plot(gbm.fut)
plot(bio.lyva.fut)
