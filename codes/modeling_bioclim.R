rm(list = ls())
#### Modelling of Herbivorous invertebrates of coral reefs ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2021

### Run modeling of all species ###

# Climatic variables with regular SST variable

# Load packages and define settings ----
library(dismo)
library(randomForest)
library(PresenceAbsence)
library(fs)
library(mgcv)
library(ecospat)

# Set seed for replicability
set.seed(2932)

# Set output folder
dir_create("modeling/bioclim")
mainf <- "modeling/bioclim/"


# Load needed functions ----
source("functions/modeling_functions.r")
source("functions/resp_curves_adapt.r")


# Load enviromental data ----
source("functions/Varload.r")
env <- var.load(folder = "crop_layers/", "env_layers.txt", "2_100")
r12 <- var.load(folder = "proj_layers/ssp126/", "env_layers.txt", "2_100")
r24 <- var.load(folder = "proj_layers/ssp245/", "env_layers.txt", "2_100")
r37 <- var.load(folder = "proj_layers/ssp370/", "env_layers.txt", "2_100")
r85 <- var.load(folder = "proj_layers/ssp585/", "env_layers.txt", "2_100")

bcc <- raster("data/env/bioclim_layers/mtemp_warmq_current_hr.tif")
bc1 <- extend(raster("data/env/bioclim_layers/mtemp_warmq_ssp126_hr.tif"), env)
bc2 <- extend(raster("data/env/bioclim_layers/mtemp_warmq_ssp245_hr.tif"), env)
bc3 <- extend(raster("data/env/bioclim_layers/mtemp_warmq_ssp370_hr.tif"), env)
bc8 <- extend(raster("data/env/bioclim_layers/mtemp_warmq_ssp585_hr.tif"), env)

bcc <- mask(bcc, env$BO21_tempmean_ss)
bc1 <- mask(bc1, env$BO21_tempmean_ss)
bc2 <- mask(bc2, env$BO21_tempmean_ss)
bc3 <- mask(bc3, env$BO21_tempmean_ss)
bc8 <- mask(bc8, env$BO21_tempmean_ss)

names(bcc) <- names(bc1) <- names(bc2) <- names(bc3) <- names(bc8) <- c("tempwq")

env <- stack(bcc, env)
r12 <- stack(bc1, r12)
r24 <- stack(bc2, r24)
r37 <- stack(bc3, r37)
r85 <- stack(bc8, r85)

rm(bcc, bc1, bc2, bc3, bc8)

env <- dropLayer(env, "BO21_tempmean_ss")
r12 <- dropLayer(r12, "BO21_tempmean_ss")
r24 <- dropLayer(r24, "BO21_tempmean_ss")
r37 <- dropLayer(r37, "BO21_tempmean_ss")
r85 <- dropLayer(r85, "BO21_tempmean_ss")

proj.rast <- list(env, r12, r24, r37, r85)
names(proj.rast) <- c("current", "r12", "r24", "r37", "r85")

# Define species list
species <- c("lyva", "eclu", "trve")

# Set output folder
# dir_create(c(paste0(mainf, species, "/models"),
#              paste0(mainf, species, "/evals")))
dir_create(paste0(mainf, species, "/evals"))


#### Start modeling ----

for (j in species) {
        
        cat("Species", j, "started. \n")
        
        outf <- paste0(mainf, j, "/")
        
        # List all files
        files <- list.files(paste0("data/", j), full.names = T)
        files <- files[grep("pa", files)]
        files <- files[-grep("folds", files)]
        
        # Run for all files
        for (i in 1:length(files)) {
                
                # Prepare data
                pts <- read.csv(files[i])
                
                blocks <- read.csv(gsub(".csv", "_folds.csv", files[i]))
                blocks <- blocks[,1]
                
                latblocks <- read.csv(gsub(".csv", "_latfolds.csv", files[i]))
                latblocks <- latblocks[,1]
                
                sp.data <- data.frame(cbind(dsp = pts$dsp,
                                            extract(env, pts[,1:2])))
                
                wgt <- wt(pts$dsp)
                
                fname <- gsub(".csv", "",
                              gsub(paste0("data/", j, "/"), "", files[i]))
                
                cat(fname, "started. \n")
                
                cv.la.table <- cv.bl.table <- 
                        data.frame(matrix(nrow = 0, ncol = 5))
                
                        
                        
                #### GBM ----
                cat("Runing GBM. \n")
                gbm.m <- gbm.auto(sp.data, wgt, blocks)
                
                # Block cross-validation
                cv.res <- data.frame(matrix(nrow = 5, ncol = 4))
                
                for (k in 1:5) {
                        
                        pred <- gbm.m$fold.fit[blocks == k]
                        
                        pred <- exp(pred)/(1 + exp(pred))
                        
                        obs <- sp.data[,1][blocks == k]
                        
                        cv.res[k,] <- model.metrics(obs, pred)
                        
                }
                
                cv.bl.table <- rbind(cv.bl.table, cbind(cv.res, model = "GBM"))
                # write.csv(cv.res, paste0(outf, "/evals/", fname,
                #                          "_gbm_block.csv"))
                
                # Latitudinal cross-validation
                cv.res[,] <- NA
                
                for (k in 1:5) {
                        
                        gbm.m.lat <- gbm.fixed(
                                sp.data[latblocks != k,],
                                gbm.x = 2:length(sp.data), gbm.y = 1,
                                tree.complexity = gbm.m$interaction.depth,
                                learning.rate = gbm.m$shrinkage,
                                n.trees = gbm.m$n.trees,
                                site.weights = wgt[latblocks != k]) 
                        
                        pred <- predict(gbm.m.lat, sp.data[latblocks == k,],
                                        type = "response")
                        
                        obs <- sp.data[,1][latblocks == k]
                        
                        cv.res[k,] <- model.metrics(obs, pred)
                        
                }
                
                cv.la.table <- rbind(cv.la.table, cbind(cv.res, model = "GBM"))
                # write.csv(cv.res, paste0(outf, "/evals/", fname,
                #                          "_gbm_latblock.csv"))
                rm(gbm.m.lat)
                
                # Get threshold
                gbm.thr <- get.thresh(sp.data$dsp,
                                      exp(gbm.m$fit)/(1 + exp(gbm.m$fit)))
                
                # Get predictions
                gbm.p <- pred.sdm(proj.rast, gbm.m, gbm.thr, outf, fname, "GBM")
                
                # Get variable importance
                gbm.vimp <- var.importance(sp.data, model = gbm.m, n = 30)
                # write.csv(gbm.vimp, paste0(outf, "/evals/", fname,
                #                          "_gbm_vimp.csv"))
                
                # # Save model
                # saveRDS(gbm.m, file = paste0(outf, "/models/", fname,
                #                              "_gbm_model.rds"))
                
                
                #### RF ----
                cat("Runing RF. \n")
                rf.m <- randomForest(dsp ~ ., data = sp.data,
                                     ntree = 1000, nodesize = 20)
                
                # Block CV
                cv.res[,] <- NA
                
                for (k in 1:5) {
                        
                        train <- sp.data[blocks != k,]
                        test <- sp.data[blocks == k,]
                        
                        rf.cv <- update(rf.m, data = train)
                        
                        pred <- predict(rf.cv, test, type = "response")
                        
                        cv.res[k,] <- model.metrics(test$dsp, pred)
                }
                
                cv.bl.table <- rbind(cv.bl.table, cbind(cv.res, model = "RF"))
                # write.csv(cv.res, paste0(outf, "/evals/", fname,
                #                          "_rf_block.csv"))
                
                # Lat block CV
                cv.res[,] <- NA
                
                for (k in 1:5) {
                        
                        train <- sp.data[latblocks != k,]
                        test <- sp.data[latblocks == k,]
                        
                        rf.cv <- update(rf.m, data = train)
                        
                        pred <- predict(rf.cv, test, type = "response")
                        
                        cv.res[k,] <- model.metrics(test$dsp, pred)
                }
                
                cv.la.table <- rbind(cv.la.table, cbind(cv.res, model = "RF"))
                # write.csv(cv.res, paste0(outf, "/evals/", fname,
                #                          "_rf_latblock.csv"))
                
                
                # Get threshold
                rf.thr <- get.thresh(sp.data$dsp,
                                     predict(rf.m, sp.data, type = "response"))
                
                # Get predictions
                rf.p <- pred.sdm(proj.rast, rf.m, rf.thr, outf, fname, "RF")
                
                # Get variable importance
                rf.vimp <- var.importance(sp.data, model = rf.m, n = 30)
                # write.csv(rf.vimp, paste0(outf, "/evals/", fname,
                #                            "_rf_vimp.csv"))
                
                
                
                #### GLM ----
                cat("Runing GLM. \n")
                glm.m <- glm(as.formula(paste("dsp ~ ", 
                                        paste(names(sp.data)[2:length(sp.data)],
                                        paste0("+ I(",
                                              names(sp.data)[2:length(sp.data)],
                                               "^2)"), collapse = "+"))),
                                    data = sp.data, weights = wgt,
                                    family = "binomial")
                
                #glm.m <- step(glm.m)
                
                # Block CV
                cv.res[,] <- NA
                
                for (k in 1:5) {
                        
                        train <- sp.data[blocks != k,]
                        test <- sp.data[blocks == k,]
                        
                        glm.cv <- update(glm.m, data = train,
                                         weights = wgt[blocks != k])
                        
                        pred <- predict(glm.cv, test, type = "response")
                        
                        cv.res[k,] <- model.metrics(test$dsp, pred)
                }
                
                cv.bl.table <- rbind(cv.bl.table, cbind(cv.res, model = "GLM"))
                # write.csv(cv.res, paste0(outf, "/evals/", fname,
                #                          "_glm_block.csv"))
                
                # Lat block CV
                cv.res[,] <- NA
                
                for (k in 1:5) {
                        
                        train <- sp.data[latblocks != k,]
                        test <- sp.data[latblocks == k,]
                        
                        glm.cv <- update(glm.m, data = train,
                                         weights = wgt[latblocks != k])
                        
                        pred <- predict(glm.cv, test, type = "response")
                        
                        cv.res[k,] <- model.metrics(test$dsp, pred)
                }
                
                cv.la.table <- rbind(cv.la.table, cbind(cv.res, model = "GLM"))
                # write.csv(cv.res, paste0(outf, "/evals/", fname,
                #                          "_glm_latblock.csv"))
                
                # Get threshold
                glm.thr <- get.thresh(sp.data$dsp,
                                      predict(glm.m, sp.data, type = "response"))
                
                # Get predictions
                glm.p <- pred.sdm(proj.rast, glm.m, glm.thr, outf, fname, "GLM")
                
                # Get variable importance
                glm.vimp <- var.importance(sp.data, model = glm.m, n = 30)
                # write.csv(glm.vimp, paste0(outf, "/evals/", fname,
                #                            "_glm_vimp.csv"))
                
                
                
                #### GAM ----
                cat("Runing GAM. \n")
                gam.m <- mgcv::gam(as.formula(paste("dsp ~",
                                 paste("s(", names(sp.data)[2:length(sp.data)],
                                          ",k=4)", collapse = "+", sep = ""))),
                                   data = sp.data,
                                   family = "binomial", 
                                   weights = wgt)
                
                # Block CV
                cv.res[,] <- NA
                
                for (k in 1:5) {
                        
                        train <- sp.data[blocks != k,]
                        test <- sp.data[blocks == k,]
                        
                        gam.cv <- update(gam.m, data = train,
                                         weights = wgt[blocks != k])
                        
                        pred <- predict(gam.cv, test, type = "response")
                        
                        cv.res[k,] <- model.metrics(test$dsp, pred)
                }
                
                cv.bl.table <- rbind(cv.bl.table, cbind(cv.res, model = "GAM"))
                # write.csv(cv.res, paste0(outf, "/evals/", fname,
                #                          "_gam_block.csv"))
                
                # Lat CV
                cv.res[,] <- NA
                
                for (k in 1:5) {
                        
                        train <- sp.data[latblocks != k,]
                        test <- sp.data[latblocks == k,]
                        
                        gam.cv <- update(gam.m, data = train,
                                         weights = wgt[latblocks != k])
                        
                        pred <- predict(gam.cv, test, type = "response")
                        
                        cv.res[k,] <- model.metrics(test$dsp, pred)
                }
                
                cv.la.table <- rbind(cv.la.table, cbind(cv.res, model = "GAM"))
                # write.csv(cv.res, paste0(outf, "/evals/", fname,
                #                          "_gam_latblock.csv"))
                
                
                # Get threshold
                gam.thr <- get.thresh(sp.data$dsp,
                                      predict(gam.m, sp.data, type = "response"))
                
                # Get predictions
                gam.p <- pred.sdm(proj.rast, gam.m, gam.thr, outf, fname, "GAM")
                
                # Get variable importance
                gam.vimp <- var.importance(sp.data, model = gam.m, n = 30)
                # write.csv(gam.vimp, paste0(outf, "/evals/", fname,
                #                            "_gam_vimp.csv"))
                
                
                #### Get models response curves ----
                cat("Getting response curves. \n")
                label <- c("Mean SST (°C)",
                           "Mean pH",
                           "Mean Chl-a (mg/m³)", 
                           "Mean salinity", 
                           "Mean depth (m)")
                
                
                
                jpeg(paste0(outf, "/evals/", fname, "_gbm_infrespcurv.jpg"),
                     units = "px", width = 1400, height = 780, res = 200,
                     quality = 100)
                par(mfrow = c(2, 3), mar = c(4, 4, 1, 1))
                inflated.response(gbm.m, sp.data[,2:length(sp.data)],
                                  label = label)
                dev.off()
                
                jpeg(paste0(outf, "/evals/", fname, "_rf_infrespcurv.jpg"),
                     units = "px", width = 1400, height = 780, res = 200,
                     quality = 100)
                par(mfrow = c(2, 3), mar = c(4, 4, 1, 1))
                inflated.response(rf.m, sp.data[,2:length(sp.data)],
                                  label = label)
                dev.off()
                
                jpeg(paste0(outf, "/evals/", fname, "_glm_infrespcurv.jpg"),
                     units = "px", width = 1400, height = 780, res = 200,
                     quality = 100)
                par(mfrow = c(2, 3), mar = c(4, 4, 1, 1))
                inflated.response(glm.m, sp.data[,2:length(sp.data)],
                                  label = label)
                dev.off()
                
                jpeg(paste0(outf, "/evals/", fname, "_gam_infrespcurv.jpg"),
                     units = "px", width = 1400, height = 780, res = 200,
                     quality = 100)
                par(mfrow = c(2, 3), mar = c(4, 4, 1, 1))
                inflated.response(gam.m, sp.data[,2:length(sp.data)],
                                  label = label)
                dev.off()
                
                
                #### Creates an ensemble ----
                cat("Getting ensembles. \n")
                ens.model <- gbm.p + rf.p + glm.p + gam.p
                ens.model <- ens.model >= 3
                writeRaster(ens.model,
                            paste0(paste0(outf, fname, "_ensemble_"),
                                   names(proj.rast), ".tif"),
                            bylayer = T)
                
                
                #### Get ens. variable importance ----
                ens.vimp <- (gbm.vimp + rf.vimp + glm.vimp + gam.vimp)/4
                write.csv(ens.vimp, paste0(outf, "/evals/", fname,
                                           "_ensemble_vimp.csv"))
                
                #### Get MESS map for this species
                cat("Getting final info... \n")
                for (k in 1:length(proj.rast)) {
                        epts <- rasterToPoints(proj.rast[[k]])
                        epts <- na.omit(epts)
                        extr <- extract(proj.rast[[k]], pts[,1:2])
                        mess.en <- ecospat.mess(epts, cbind(pts[,1:2], extr),
                                                w="default")
                        writeRaster(rasterFromXYZ(mess.en[,c(1,2,3)]),
                                    paste0(outf, "/evals/", fname, "_",
                                           names(proj.rast)[k],
                                           "_mess.tif"))
                        writeRaster(rasterFromXYZ(mess.en[,c(1,2,4)]),
                                    paste0(outf, "/evals/", fname, "_",
                                           names(proj.rast)[k],
                                           "_messw.tif"))
                        writeRaster(rasterFromXYZ(mess.en[,c(1,2,5)]),
                                    paste0(outf, "/evals/", fname, "_",
                                           names(proj.rast)[k],
                                           "_messneg.tif"))
                }
                
                #### Get Boyce values for all models
                boyce <- data.frame(model = c("GBM", "RF", "GLM", "GAM"),
                                    boyce = NA)
                
                for (k in 1:4) {
                        
                        pred <- predict(list(gbm.m, rf.m, glm.m, gam.m)[[k]],
                                        sp.data,
                                        type = "response")
                        
                        boyce[k,2] <- ecospat.boyce(pred,
                                                    pred[sp.data$dsp == 1],
                                                    PEplot = F)$Spearman.cor
                }
                
                write.csv(boyce, paste0(outf, "/evals/", fname,
                                        "_boyce.csv"))
                
                #### Save CV results
                names(cv.la.table) <- names(cv.bl.table) <- 
                        c("AUC", "Specificity", "Sensitivity", "TSS", "Model")
                
                write.csv(cv.bl.table, paste0(outf, "/evals/", fname,
                                              "_blockcv.csv"))
                
                write.csv(cv.la.table, paste0(outf, "/evals/", fname,
                                              "_latblockcv.csv"))
                
                #### Save variable importance results
                vimp.final <- rbind(
                        cbind(gbm.vimp, model = "GBM"),
                        cbind(rf.vimp, model = "RF"),
                        cbind(glm.vimp, model = "GLM"),
                        cbind(gam.vimp, model = "GAM")
                )
                
                write.csv(vimp.final, paste0(outf, "/evals/", fname,
                                             "_vimp.csv"))
                
                
                cat(fname, "finished. \n")
                
        }
        
        cat("Species", j, "finished. \n")
}

system("shutdown -s")