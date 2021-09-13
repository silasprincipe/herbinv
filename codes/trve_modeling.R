#### Modelling of Herbivorous invertebrates of coral reefs ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2021

### Modeling - Sea-urchins ###

#### Echinometra lucunter models ####

# Climatic variables with regular SST variable ----

# Load packages ----
library(dismo)
library(fs)
library(PresenceAbsence)

# Set seed for replicability
set.seed(2932)

# Set output folder
out.f <- "modeling/regular/trve"
dir_create(out.f)
dir_create(paste0(out.f, "/", c("proj", "models", "eval")))

# Creates a log file
write(paste("Modeling procedure started at:", date()),
      file = paste0(out.f, "/trve_log_file.txt"))

# GBM ----
# Train the initial model and find the optimal one
gbm.m <- gbm.step(sp.data$data, gbm.x = 2:ncol(sp.data$data), gbm.y = 1,
                  fold.vector = sp.data$blocks, n.folds = 5,
                  site.weights = sp.data$weigths,
                  keep.fold.fit = T, learning.rate = 0.002,
                  tree.complexity = 5)

# See the response curves
par(mfrow = c(2,3))
resp.curves(gbm.m, sp.data$data[,2:ncol(sp.data$data)], method = "stat3")

# See the variable importance
gbm.m$contributions

var.importance(sp.data$data, model = gbm.m)

# Get other stats of the blocks CV
cv.res <- data.frame(matrix(nrow = 5, ncol = 5))
colnames(cv.res) <- c("AUC", "Spec", "Sens", "TSS", "Boyce")

for (k in 1:5) {
        
        pred <- gbm.m$fold.fit[sp.data$blocks == k]
        
        pred <- exp(pred)/(1 + exp(pred))
        
        obs <- sp.data$data[,1][sp.data$blocks == k]
        
        td <- data.frame(id = 1,
                         obs = obs,
                         pred = pred)
        
        thresh <- optimal.thresholds(td,
                                     opt.methods = "MaxSens+Spec")$pred
        
        cmat <- cmx(td, threshold = thresh)
        
        cv.res[k,1] <- PresenceAbsence::auc(td)[,1]
        
        cv.res[k,2] <- PresenceAbsence::specificity(cmat)[,1]
        
        cv.res[k,3] <- PresenceAbsence::sensitivity(cmat)[,1]
        
        cv.res[k,4] <- (cv.res[k,2] + cv.res[k,3]) - 1
        
        cv.res[k,5] <- ecospat::ecospat.boyce(pred,
                                                   pred[obs == 1],
                                                   PEplot = F)$Spearman.cor
        
        
}

# See the results
cv.res
apply(cv.res, 2, mean)

# Cross-validation for latitudinal blocks
cv.res.lat <- data.frame(matrix(nrow = 5, ncol = 5))
colnames(cv.res.lat) <- c("AUC", "Spec", "Sens", "TSS", "Boyce")

gbm.m.lat <- list()

for (k in 1:5) {
        
        gbm.m.lat[[k]] <- gbm.fixed(
                sp.data$data[sp.data$lat_blocks != k,],
                gbm.x = 2:ncol(sp.data$data), gbm.y = 1,
                tree.complexity = gbm.m$interaction.depth,
                learning.rate = gbm.m$shrinkage,
                n.trees = gbm.m$n.trees,
                site.weights = sp.data$weigths[sp.data$lat_blocks != k]) 
        
        test.data <- sp.data$data[sp.data$lat_blocks == k,]
        
        pred <- predict(gbm.m.lat[[k]], test.data[,2:ncol(test.data)],
                        type = "response")
        
        
        td <- data.frame(id = 1,
                         obs = test.data[,1],
                         pred = pred)
        
        thresh <- optimal.thresholds(td,
                                     opt.methods = "MaxSens+Spec")$pred
        
        cmat <- cmx(td, threshold = thresh)
        
        cv.res.lat[k,1] <- PresenceAbsence::auc(td)[,1]
        
        cv.res.lat[k,2] <- PresenceAbsence::specificity(cmat)[,1]
        
        cv.res.lat[k,3] <- PresenceAbsence::sensitivity(cmat)[,1]
        
        cv.res.lat[k,4] <- (cv.res.lat[k,2] + cv.res.lat[k,3]) - 1
        
        cv.res.lat[k,5] <- ecospat::ecospat.boyce(pred,
                                              pred[test.data[,1] == 1],
                                              PEplot = F)$Spearman.cor
        
        rm(test.data, pred, td, thresh, cmat)
        
}

# See the results
cv.res.lat
apply(cv.res.lat, 2, mean)

# Get final model data
final.result <- list(gbm.m, cv.res, cv.res.lat)
names(final.result) <- c("model", "block_cv", "latitudinal_cv")

# Get final model results
td <- data.frame(id = 1,
                 obs = sp.data$data[,1],
                 pred = predict(gbm.m, sp.data$data, type = "response"))

thresh <- optimal.thresholds(td,
                             opt.methods = "MaxSens+Spec")$pred

final.result$tss_thresh <- thresh

final.result$boyce <- ecospat::ecospat.boyce(td$pred,
                                             td$pred[sp.data$data[,1] == 1],
                                             PEplot = F)

final.result$varimp <- var.importance(sp.data$data, model = gbm.m, get.mean = F)


# Get predictions
cur.full <- pred.sdm(env, gbm.m)
cur.thre <- pred.sdm(env, gbm.m, final.result$tss_thresh) 

writeRaster(cur.full, paste0(out.f, "/proj/GBM_cur_trve.tif"))
writeRaster(cur.thre, paste0(out.f, "/proj/GBM_cur_trve_TSS.tif"))

fut.models <- list(r12, r24, r37, r85)

for (i in 1:4) {
        
        ssp <- paste0("ssp", c(126, 245, 370, 585))[i]
        
        fut.full <- pred.sdm(fut.models[[i]], gbm.m)
        fut.thre <- pred.sdm(fut.models[[i]], gbm.m, final.result$tss_thresh) 
        
        writeRaster(fut.full, paste0(out.f, "/proj/GBM_", ssp, "_trve.tif"))
        writeRaster(fut.thre, paste0(out.f, "/proj/GBM_", ssp, "_trve_TSS.tif"))
}

# Plot to see
par(mfrow = c(2, 3))
files <- list.files(paste0(out.f, "/proj/"), pattern = "TSS", full.names = T)
lapply(files, function(x){plot(raster(x))})

# Save model data
save.mod(final.result, out.f, "GBM", "trve")

rm.obj <- ls()
rm.obj <- rm.obj[!rm.obj %in% c("env", "r12", "r24", "r37", "r85",
                                "var.importance", "resp.curves", "save.mod",
                                "pred.sdm", "sp.data", "out.f", "env.alt")]
rm(list = rm.obj)

# RF ----
# Train the initial model and find the optimal one
library(caret)
library(randomForest)

# Create a list to be used in the tunning
indexlist <- list(seq(1:nrow(sp.data$data))[sp.data$blocks != 1],
                  seq(1:nrow(sp.data$data))[sp.data$blocks != 2],
                  seq(1:nrow(sp.data$data))[sp.data$blocks != 3],
                  seq(1:nrow(sp.data$data))[sp.data$blocks != 4],
                  seq(1:nrow(sp.data$data))[sp.data$blocks != 5])

fitControl <- trainControl(
        method = "cv",
        number = 5,
        index = indexlist)

rf1 <- train(dsp ~ ., data = sp.data$data, 
             method = "rf", 
             trControl = fitControl,
             ## This last option is actually one
             ## for gbm() that passes through
             verbose = T)

mtry <- rf1$bestTune$mtry

# Get full model
rf.m <- randomForest(x = sp.data$data[, 2:ncol(sp.data$data)],
                     y = sp.data$data[, 1],
                     mtry = mtry)

# See the response curves
par(mfrow = c(2,3))
resp.curves(rf.m, sp.data$data[,2:ncol(sp.data$data)], method = "stat3")

# See the variable importance
rf.m$importance

var.importance(sp.data$data, model = rf.m)

# Get other stats of the blocks CV
cv.res <- rf.crossval(sp.data$data, sp.data$blocks, mtry)

cv.res

cv.res.lat <- rf.crossval(sp.data$data, sp.data$lat_blocks, mtry)

cv.res.lat

# Get final model data
final.result <- list(rf.m, cv.res, cv.res.lat)
names(final.result) <- c("model", "block_cv", "latitudinal_cv")

# Get final model results
td <- data.frame(id = 1,
                 obs = sp.data$data[,1],
                 pred = predict(rf.m, sp.data$data[,2:ncol(sp.data$data)],
                                type = "response"))

thresh <- optimal.thresholds(td,
                             opt.methods = "MaxSens+Spec")$pred

final.result$tss_thresh <- thresh

final.result$boyce <- ecospat::ecospat.boyce(fit = td$pred,
                                             obs = td$pred[sp.data$data[,1] == 1],
                                             PEplot = F)

final.result$varimp <- var.importance(sp.data$data, model = rf.m, get.mean = F)


# Get predictions
cur.full <- pred.sdm(env, rf.m)
cur.thre <- pred.sdm(env, rf.m, final.result$tss_thresh) 

writeRaster(cur.full, paste0(out.f, "/proj/RF_cur_trve.tif"))
writeRaster(cur.thre, paste0(out.f, "/proj/RF_cur_trve_TSS.tif"))

fut.models <- list(r12, r24, r37, r85)

for (i in 1:4) {
        
        ssp <- paste0("ssp", c(126, 245, 370, 585))[i]
        
        fut.full <- pred.sdm(fut.models[[i]], rf.m)
        fut.thre <- pred.sdm(fut.models[[i]], rf.m, final.result$tss_thresh) 
        
        writeRaster(fut.full, paste0(out.f, "/proj/RF_", ssp, "_trve.tif"))
        writeRaster(fut.thre, paste0(out.f, "/proj/RF_", ssp, "_trve_TSS.tif"))
}

# Plot to see
par(mfrow = c(2, 3))
files <- list.files(paste0(out.f, "/proj/"), pattern = "TSS", full.names = T)
files <- files[grep("RF", files)]
lapply(files, function(x){plot(raster(x))})

# Save model data
save.mod(final.result, out.f, "RF", "trve")

rm.obj <- ls()
rm.obj <- rm.obj[!rm.obj %in% c("env", "r12", "r24", "r37", "r85",
                                "var.importance", "resp.curves", "save.mod",
                                "pred.sdm", "sp.data", "out.f", "env.alt")]
rm(list = rm.obj)



# GAM ----
source("functions/modeling_functions.R")

gam.formula <- paste("dsp ~", paste0('s(', colnames(sp.data$data[,2:ncol(sp.data$data)]),', k = 5)', collapse=' + '))

gam.m <- mgcv::gam(as.formula(gam.formula),
                   data = sp.data$data,
                   family = binomial(link = 'logit'), 
                   weights = sp.data$weigths,
                   select=TRUE, method='REML')

summary(gam.m)

par(mfrow = c(2,3))
resp.curves(gam.m, sp.data$data[,2:ncol(sp.data$data)], method = "stat3")

# See the variable importance
var.importance(sp.data$data, model = gam.m)


# Get other stats of the blocks CV
cv.res <- gam.crossval(gam.m$formula, sp.data$blocks, sp.data$data, sp.data$weigths)

cv.res

cv.res.lat <- gam.crossval(gam.m$formula, sp.data$lat_blocks, sp.data$data, sp.data$weigths)

cv.res.lat


# Get final model data
final.result <- list(gam.m, cv.res, cv.res.lat)
names(final.result) <- c("model", "block_cv", "latitudinal_cv")

# Get final model results
td <- data.frame(id = 1,
                 obs = sp.data$data[,1],
                 pred = predict(gam.m, sp.data$data[,2:ncol(sp.data$data)],
                                type = "response"))

thresh <- optimal.thresholds(td,
                             opt.methods = "MaxSens+Spec")$pred

final.result$tss_thresh <- thresh

final.result$boyce <- ecospat::ecospat.boyce(fit = td$pred,
                                             obs = td$pred[sp.data$data[,1] == 1],
                                             PEplot = F)

final.result$varimp <- var.importance(sp.data$data, model = gam.m, get.mean = F)


# Get predictions
cur.full <- pred.sdm(env, gam.m)
cur.thre <- pred.sdm(env, gam.m, final.result$tss_thresh) 

writeRaster(cur.full, paste0(out.f, "/proj/GAM_cur_trve.tif"))
writeRaster(cur.thre, paste0(out.f, "/proj/GAM_cur_trve_TSS.tif"))

fut.models <- list(r12, r24, r37, r85)

for (i in 1:4) {
        
        ssp <- paste0("ssp", c(126, 245, 370, 585))[i]
        
        fut.full <- pred.sdm(fut.models[[i]], gam.m)
        fut.thre <- pred.sdm(fut.models[[i]], gam.m, final.result$tss_thresh) 
        
        writeRaster(fut.full, paste0(out.f, "/proj/GAM_", ssp, "_trve.tif"))
        writeRaster(fut.thre, paste0(out.f, "/proj/GAM_", ssp, "_trve_TSS.tif"))
}

# Plot to see
par(mfrow = c(2, 3))
files <- list.files(paste0(out.f, "/proj/"), pattern = "TSS", full.names = T)
files <- files[grep("GAM", files)]
lapply(files, function(x){plot(raster(x))})

# Save model data
save.mod(final.result, out.f, "GAM", "trve")

rm.obj <- ls()
rm.obj <- rm.obj[!rm.obj %in% c("env", "r12", "r24", "r37", "r85",
                                "var.importance", "resp.curves", "save.mod",
                                "pred.sdm", "sp.data", "out.f", "env.alt")]
rm(list = rm.obj)


# GLM ----
source("functions/modeling_functions.R")

glm.formula <- dsp ~ BO21_lightbotmax_bdmean + I(BO21_lightbotmax_bdmean^2) +
        BO_ph + I(BO_ph^2) + BO21_chlomean_ss + I(BO21_chlomean_ss^2) +
        BO21_silicatemean_ss + I(BO21_silicatemean_ss^2) +
        BO21_salinitymean_ss + I(BO21_salinitymean_ss^2) +
        BO21_tempmean_ss + I(BO21_tempmean_ss^2)

glm.m <- glm(as.formula(glm.formula),
                   data = sp.data$data,
                   family = binomial(link = 'logit'), 
                   weights = sp.data$weigths)

glm.m.step <- step(glm.m)

summary(glm.m.step)

par(mfrow = c(2,3))
resp.curves(glm.m.step, sp.data$data[,2:ncol(sp.data$data)], method = "stat3")

glm.formula <- dsp ~ BO21_lightbotmax_bdmean + BO_ph + I(BO_ph^2) + 
        I(BO21_chlomean_ss^2) + BO21_silicatemean_ss + I(BO21_tempmean_ss^2)

glm.m <- glm(as.formula(glm.formula),
             data = sp.data$data,
             family = binomial(link = 'logit'), 
             weights = sp.data$weigths)

summary(glm.m)

par(mfrow = c(2,3))
resp.curves(glm.m, sp.data$data[,2:ncol(sp.data$data)], method = "stat3")

# See the variable importance
var.importance(sp.data$data, model = glm.m)


# Get other stats of the blocks CV
cv.res <- glm.crossval(glm.m$formula, sp.data$blocks, sp.data$data, sp.data$weigths)

cv.res

cv.res.lat <- glm.crossval(glm.m$formula, sp.data$lat_blocks, sp.data$data, sp.data$weigths)

cv.res.lat


# Get final model data
final.result <- list(glm.m, cv.res, cv.res.lat)
names(final.result) <- c("model", "block_cv", "latitudinal_cv")

# Get final model results
td <- data.frame(id = 1,
                 obs = sp.data$data[,1],
                 pred = predict(glm.m, sp.data$data[,2:ncol(sp.data$data)],
                                type = "response"))

thresh <- optimal.thresholds(td,
                             opt.methods = "MaxSens+Spec")$pred

final.result$tss_thresh <- thresh

final.result$boyce <- ecospat::ecospat.boyce(fit = td$pred,
                                             obs = td$pred[sp.data$data[,1] == 1],
                                             PEplot = F)

final.result$varimp <- var.importance(sp.data$data, model = glm.m, get.mean = F)


# Get predictions
cur.full <- pred.sdm(env, glm.m)
cur.thre <- pred.sdm(env, glm.m, final.result$tss_thresh) 

writeRaster(cur.full, paste0(out.f, "/proj/GLM_cur_trve.tif"))
writeRaster(cur.thre, paste0(out.f, "/proj/GLM_cur_trve_TSS.tif"))

fut.models <- list(r12, r24, r37, r85)

for (i in 1:4) {
        
        ssp <- paste0("ssp", c(126, 245, 370, 585))[i]
        
        fut.full <- pred.sdm(fut.models[[i]], glm.m)
        fut.thre <- pred.sdm(fut.models[[i]], glm.m, final.result$tss_thresh) 
        
        writeRaster(fut.full, paste0(out.f, "/proj/GLM_", ssp, "_trve.tif"))
        writeRaster(fut.thre, paste0(out.f, "/proj/GLM_", ssp, "_trve_TSS.tif"))
}

# Plot to see
par(mfrow = c(2, 3))
files <- list.files(paste0(out.f, "/proj/"), pattern = "TSS", full.names = T)
files <- files[grep("GLM", files)]
lapply(files, function(x){plot(raster(x))})

# Save model data
save.mod(final.result, out.f, "GLM", "trve")

rm.obj <- ls()
rm.obj <- rm.obj[!rm.obj %in% c("env", "r12", "r24", "r37", "r85",
                                "var.importance", "resp.curves", "save.mod",
                                "pred.sdm", "sp.data", "out.f", "env.alt")]
rm(list = rm.obj)

### END of trve modeling
write(paste("Modeling procedure concluded at:", date()),
      file = paste0(out.f, "/trve_log_file.txt"), append = T)

write(capture.output(sessionInfo()),
      file = paste0(out.f, "/trve_log_file.txt"), append = T)
