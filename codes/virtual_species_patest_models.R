#### Modelling of Herbivorous invertebrates of coral reefs ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2021

### Virtual species modeling ###
### PA method test ###

# Load packages and define settings ----
library(dismo)
library(randomForest)
library(PresenceAbsence)
library(fs)
library(MASS)
library(virtualspecies)

# Set seed for replicability
set.seed(2932)

# Load needed functions ----
source("functions/vsp_modeling.r")

# Load enviromental data ----
source("functions/Varload.r")
env <- var.load(folder = "crop_layers/", "env_layers.txt", "2_100")
r12 <- var.load(folder = "proj_layers/ssp126/", "env_layers.txt", "2_100")
r24 <- var.load(folder = "proj_layers/ssp245/", "env_layers.txt", "2_100")
r37 <- var.load(folder = "proj_layers/ssp370/", "env_layers.txt", "2_100")
r85 <- var.load(folder = "proj_layers/ssp585/", "env_layers.txt", "2_100")

# proj.rast <- list(env, r12, r24, r37, r85)
# names(proj.rast) <- c("current", "r12", "r24", "r37", "r85")
proj.rast <- list(env)
names(proj.rast) <- "current"

# Folder settings ----
setwd("data/vsp")

# Set output folder
dir_create("modeling/models")
out.f <- "modeling/"

# Get environmental data to use with the respcurve function ----
env.dat <- data.frame(var = names(env),
                      min = minValue(env),
                      max = maxValue(env),
                      mean = as.numeric(cellStats(env, stat = "mean")))

# Get original response curves
orig.respcurv <- data.frame(matrix(nrow = 50, ncol = length(env.dat$var)))

for (z in 1:ncol(orig.respcurv)) {
        orig.respcurv[,z] <- seq(env.dat$min[z], env.dat$max[z], length = 50)
}

names(orig.respcurv) <- env.dat$var

orig.respcurv$BO_ph <- dnorm(orig.respcurv$BO_ph, mean = 8.1, sd = 0.2)
orig.respcurv$BO21_chlomean_ss <- 
        logisticFun(orig.respcurv$BO21_chlomean_ss, beta = 0.8, alpha = 0.1)
# orig.respcurv$BO21_silicatemean_ss <- 
#         logisticFun(orig.respcurv$BO21_silicatemean_ss, beta = 25, alpha = 1) 
orig.respcurv$BO21_salinitymean_ss <- 
        dnorm(orig.respcurv$BO21_salinitymean_ss, mean = 35, sd = 10)

orig.rc.vsp1 <- orig.rc.vsp2 <- orig.respcurv

orig.rc.vsp1$BO21_tempmean_ss <- 
        custnorm(orig.rc.vsp1$BO21_tempmean_ss, mean = 25, diff = 10, prob = 0.99)
orig.rc.vsp2$BO21_tempmean_ss <- 
        custnorm(orig.rc.vsp2$BO21_tempmean_ss, mean = 26.5, diff = 7.5, prob = 0.99)


orig.rc.vsp1 <- apply(orig.rc.vsp1, 2, function(x){
        (x - min(x))/(max(x) - min(x))
})
orig.rc.vsp2 <- apply(orig.rc.vsp2, 2, function(x){
        (x - min(x))/(max(x) - min(x))
})

# Load original rasters
vsp1 <- raster("vsp1_pa.tif")
vsp2 <- raster("vsp2_pa.tif")

# Future projections
vsp1.fut <- stack(list.files(pattern = "vsp1_r"))
vsp2.fut <- stack(list.files(pattern = "vsp2_r"))


#### Start modeling ----

# List all files
files <- list.files("padata", full.names = T)
files <- files[-grep("folds", files)]

# Creates list to hold internal CV
cv.total <- data.frame(matrix(nrow = length(files), ncol = 8))
names(cv.total) <- c("AUC", "Spec", "Sens", "TSS",
                     paste0(c("AUC", "Spec", "Sens", "TSS"), "_SD"))

cv.list <- list(cv.total, cv.total, cv.total, cv.total)
names(cv.list) <- c("GBM", "RF", "GLM", "GAM")
rm(cv.total)

# Creates a list to hold metrics of ensemble
ens.metrics <- data.frame(matrix(nrow = length(files), ncol = 8))
names(ens.metrics) <- c("model", "jac", "I", "D",
                        "tss", "auc", "spec", "sens")

# Creates a list to hold respcurves results
rc.results <- list()

# Creates a list to hold resp curves spearman cor
ens.spearm <- data.frame(matrix(nrow = length(files), ncol = nlayers(env)))
names(ens.spearm) <- names(env)

# Creates a list to hold jaccard comparison of future
ens.jacc <- data.frame(matrix(nrow = length(files), ncol = length(proj.rast)))
names(ens.jacc) <- names(proj.rast)

# Run for all files
for (i in 1:length(files)) {
        
        pts <- read.csv(files[i])
        
        blocks <- read.csv(gsub(".csv", "_folds.csv", files[i]))
        blocks <- blocks[,1]
        
        sp.data <- data.frame(cbind(dsp = pts$dsp, extract(env, pts[,1:2])))
        
        wgt <- wt(pts$dsp)
        
        fname <- gsub(".csv", "", gsub("padata/", "", files[i]))
        
        cat(fname, "started. \n")
        
        #### GBM ----
        gbm.m <- gbm.auto(sp.data, wgt, blocks)
        
        cv.res <- data.frame(matrix(nrow = 5, ncol = 4))
        
        for (k in 1:5) {
                
                pred <- gbm.m$fold.fit[blocks == k]
                
                pred <- exp(pred)/(1 + exp(pred))
                
                obs <- sp.data[,1][blocks == k]
                
                cv.res[k,] <- model.metrics(obs, pred)
                
        }
        
        cv.list$GBM[i,1:4] <- apply(cv.res, 2, mean)
        cv.list$GBM[i,5:8] <- apply(cv.res, 2, sd)
        
        gbm.thr <- get.thresh(sp.data$dsp, exp(gbm.m$fit)/(1 + exp(gbm.m$fit)))
        
        #### RF ----
        rf.m <- randomForest(dsp ~ ., data = sp.data)
        
        cv.res[,] <- NA
        
        for (k in 1:5) {
                
                train <- sp.data[blocks != k,]
                test <- sp.data[blocks == k,]
                
                rf.cv <- update(rf.m, data = train)
                
                pred <- predict(rf.cv, test, type = "response")
                
                cv.res[k,] <- model.metrics(test$dsp, pred)
        }
        
        cv.list$RF[i,1:4] <- apply(cv.res, 2, mean)
        cv.list$RF[i,5:8] <- apply(cv.res, 2, sd)
        
        rf.thr <- get.thresh(sp.data$dsp,
                             predict(rf.m, sp.data, type = "response"))
        
        #### GLM ----
        
        glm.m <- glm(paste("dsp ~",
                           paste("poly(", names(sp.data)[2:length(sp.data)], ",2)",
                                 collapse = "+", sep = "")),
                     data = sp.data, weights = wgt, family = "binomial")
        
        # glm.m <- glm(paste("dsp ~",
        #                    paste(names(sp.data)[2:length(sp.data)], " + I(",
        #                          names(sp.data)[2:length(sp.data)], "^2)",
        #                          collapse = "+", sep = "")),
        #              data = sp.data, weights = wgt, family = "binomial")
        
        glm.step <- step(glm.m)
        
        cv.res[,] <- NA
        
        for (k in 1:5) {
                
                train <- sp.data[blocks != k,]
                test <- sp.data[blocks == k,]
                
                glm.cv <- update(glm.step, data = train,
                                 weights = wgt[blocks != k])
                
                pred <- predict(glm.cv, test, type = "response")
                
                cv.res[k,] <- model.metrics(test$dsp, pred)
        }
        
        cv.list$GLM[i,1:4] <- apply(cv.res, 2, mean)
        cv.list$GLM[i,5:8] <- apply(cv.res, 2, sd)
        
        glm.thr <- get.thresh(sp.data$dsp,
                              predict(glm.step, sp.data, type = "response"))
        
        #### GAM ----
        gam.m <- mgcv::gam(as.formula(paste("dsp ~",
                                            paste("s(", names(sp.data)[2:length(sp.data)],
                                                  ",k=4)", collapse = "+", sep = ""))),
                           data = sp.data,
                           family = binomial(link = 'logit'), 
                           weights = wgt)
        
        cv.res[,] <- NA
        
        for (k in 1:5) {
                
                train <- sp.data[blocks != k,]
                test <- sp.data[blocks == k,]
                
                gam.cv <- update(gam.m, data = train,
                                 weights = wgt[blocks != k])
                
                pred <- predict(gam.cv, test, type = "response")
                
                cv.res[k,] <- model.metrics(test$dsp, pred)
        }
        
        cv.list$GAM[i,1:4] <- apply(cv.res, 2, mean)
        cv.list$GAM[i,5:8] <- apply(cv.res, 2, sd)
        
        gam.thr <- get.thresh(sp.data$dsp,
                              predict(gam.m, sp.data, type = "response"))
        
        ### Get response curves ----
        if (i == 1) {
                rc.results[[1]] <- resp.curve(gbm.m, fname, env.dat)
                rc.results[[2]] <- resp.curve(rf.m, fname, env.dat)
                rc.results[[3]] <- resp.curve(glm.m, fname, env.dat)
                rc.results[[4]] <- resp.curve(gam.m, fname, env.dat)
        } else{
                rc.results[[1]] <- rbind(rc.results[[1]],
                                         resp.curve(gbm.m, fname, env.dat))
                rc.results[[2]] <- rbind(rc.results[[2]],
                                         resp.curve(rf.m, fname, env.dat))
                rc.results[[3]] <- rbind(rc.results[[3]],
                                         resp.curve(glm.m, fname, env.dat))
                rc.results[[4]] <- rbind(rc.results[[4]],
                                         resp.curve(gam.m, fname, env.dat))
        }
        
        ### Make projections ----
        pred.gbm <- pred.sdm(proj.rast, gbm.m, gbm.thr)
        pred.rf <- pred.sdm(proj.rast, rf.m, rf.thr)
        pred.glm <- pred.sdm(proj.rast, glm.m, glm.thr)
        pred.gam <- pred.sdm(proj.rast, gam.m, gam.thr)
        
        ### Creates an ensemble ----
        ens.model <- pred.gbm + pred.rf + pred.glm + pred.gam
        ens.model <- ens.model >= 3
        
        if (grepl("vsp1", fname)) {
                ens.metrics[i,2:8] <- compare.models(ens.model[[1]], vsp1)
                ens.metrics[i,1] <- fname
        } else {
                ens.metrics[i,2:8] <- compare.models(ens.model[[1]], vsp2)
                ens.metrics[i,1] <- fname
        }
        
        ### Get ensemble response curves ----
        ens.resp.curve <- (resp.curve(gbm.m, fname, env.dat)[,1:nlayers(env)] +
                resp.curve(rf.m, fname, env.dat)[,1:nlayers(env)] +
                resp.curve(glm.m, fname, env.dat)[,1:nlayers(env)] +
                resp.curve(gam.m, fname, env.dat)[,1:nlayers(env)])/4
        
        for (z in 1:ncol(ens.resp.curve)) {
                if (grepl("vsp1", fname)) {
                        ens.spearm[i,z] <- cor(ens.resp.curve[,z],
                                               orig.rc.vsp1[,z],
                                               method = "spearman") 
                } else {
                        ens.spearm[i,z] <- cor(ens.resp.curve[,z],
                                               orig.rc.vsp2[,z],
                                               method = "spearman")
                }
                
        }
        
        # ### Get jaccard values for the current and future projections
        # if (grepl("vsp1", fname)) {
        #         ens.jacc[i,1] <- jaccard.eval(ens.model[[1]], vsp1)
        #         for (z in 1:4) {
        #                 ens.jacc[i, (z + 1)] <-
        #                         jaccard.eval(ens.model[[(z + 1)]], vsp1.fut[[z]])
        #         }
        # } else {
        #         ens.jacc[i,1] <- jaccard.eval(ens.model[[1]], vsp2)
        #         for (z in 1:4) {
        #                 ens.jacc[i, (z + 1)] <-
        #                         jaccard.eval(ens.model[[(z + 1)]], vsp2.fut[[z]])
        #         }
        # }
        # 
        saveRDS(list(gbm.m, rf.m, glm.step, gam.m, 
                     gbm.thr, rf.thr, glm.thr, gam.thr),
                file = paste0("modeling/models/", fname, ".rds"))
}

# Save cross-validation results
cv.list$GBM$model <- gsub(".csv", "", gsub("padata/", "", files))
cv.list$RF$model <- gsub(".csv", "", gsub("padata/", "", files))
cv.list$GLM$model <- gsub(".csv", "", gsub("padata/", "", files))
cv.list$GAM$model <- gsub(".csv", "", gsub("padata/", "", files))

write.csv(cv.list$GBM, paste0(out.f, "gbm_cv.csv"), row.names = F)
write.csv(cv.list$RF, paste0(out.f, "rf_cv.csv"), row.names = F)
write.csv(cv.list$GLM, paste0(out.f, "glm_cv.csv"), row.names = F)
write.csv(cv.list$GAM, paste0(out.f, "gam_cv.csv"), row.names = F)

# Save response curves results
write.csv(rc.results[[1]], paste0(out.f, "gbm_respcurv.csv"), row.names = F)
write.csv(rc.results[[2]], paste0(out.f, "rf_respcurv.csv"), row.names = F)
write.csv(rc.results[[3]], paste0(out.f, "glm_respcurv.csv"), row.names = F)
write.csv(rc.results[[4]], paste0(out.f, "gam_respcurv.csv"), row.names = F)

# Save ensemble metrics results
ens.jacc$model <- ens.spearm$model <- ens.metrics$model
write.csv(ens.metrics, paste0(out.f, "ens_metrics.csv"), row.names = F)
# write.csv(ens.jacc, paste0(out.f, "ens_jaccard.csv"), row.names = F)
write.csv(ens.spearm, paste0(out.f, "ens_respcurves_cor.csv"), row.names = F)

setwd("../..")