#### Modelling of Herbivorous invertebrates of coral reefs ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2021

### Modeling - Sea-urchins ###
# Climatic variables with regular SST variable

# Load packages and define settings ----
library(dismo)
library(randomForest)
library(PresenceAbsence)
library(crayon)
library(fs)
library(MASS)
library(caret)

# Set seed for replicability
set.seed(2932)

# Set output folder
dir_create("modeling/regular")
out.f <- "modeling/regular/"

# Creates a log file
write(paste("Modeling procedure started at:", date()),
      file = paste0(out.f, "log_file.txt"))




# Prepare data for modeling ----

# Define species list
sp.list <- c("lyva", "eclu", "trve")

# Load environmental data
source("functions/Varload.r")
env <- var.load(folder = "crop_layers/", layers = "env_layers.txt")

r12 <- var.load(folder = "proj_layers/ssp126/", layers = "env_layers.txt")

r24 <- var.load(folder = "proj_layers/ssp245/", layers = "env_layers.txt")

r37 <- var.load(folder = "proj_layers/ssp370/", layers = "env_layers.txt")

r85 <- var.load(folder = "proj_layers/ssp585/", layers = "env_layers.txt")

### Open species data in a list
# First we create a function to open all species data
prep.data <- function(species){
        
        # Load species data
        pts <- read.csv(paste0("data/", species, "/", species, "_pts.csv"))
        
        dsp <- read.csv(paste0("data/", species, "/", species, "_dsp.csv"))
        dsp[is.na(dsp$dsp),] <- 0 # this data is in the same format to use with
                                  # BIOMOD, where PA is treated as NA
        
        # Load split table (block cross validation)
        blocks <- read.csv(paste0("data/", species, "/",
                                  species, "_foldindex.csv"))
        blocks <- blocks[,1]
        
        # The same for latitudinal blocks
        blocks.lat <- read.csv(paste0("data/", species, "/",
                                      species, "_foldindex_lat.csv"))
        blocks.lat <- blocks.lat[,1]
        
        # Extract the data for the points
        data <- raster::extract(env, pts)
        
        # Bind with prsence/absence data
        data <- cbind(dsp, data)
        
        # Creates a weigthed response
        wt <- function(resp, prev = 0.5){
                
                pres <- sum(resp)
                
                absc <- length(resp) - pres
                
                wgth <- rep(1, length(resp))
                
                wgth[which(resp > 0)] <- (prev * absc)/(pres * (1 - prev))
                
                wgth <- round(wgth)
                
                wgth
        }
        
        wt.dsp <- wt(dsp$dsp)
        
        # Generate a list with all data
        full <- list(
                "points" = pts,
                "pa" = dsp,
                "blocks" = blocks,
                "lat_blocks" = blocks.lat,
                "data" = data,
                "weigths" = wt.dsp
        )
        
        full
}

# Then run for the species list
sp.data <- lapply(sp.list, prep.data)

names(sp.data) <- sp.list



# Run models in loop ----
for (i in 1:length(sp.list)) {
        
        species <- sp.list[i]
        
        fin.f <- paste0(out.f, species, "/")
        
        # Creates needed folders
        dir_create(paste0(fin.f, c("models", "projections", "eval")))
        
        
        cat("Modeling of", species, "started at", date(), "\n")
        
        write(paste("Modeling of", species, "started at", date()),
              file = paste0(out.f, "log_file.txt"),
              append = TRUE)
        
        ##### Boosted Regression Trees ----
        
        cat(green("===== Fitting GBM model ===== \n"))
        
        ## Find optimal model
        gbm.m <- gbm.step(sp.data[[species]]$data, gbm.x = 2:8, gbm.y = 1,
                          fold.vector = sp.data[[species]]$blocks, n.folds = 5,
                          site.weights = sp.data[[species]]$weigths,
                          keep.fold.fit = T)
        
        
        
        ## Cross validate for latitudinal blocks
        cat("Latitudinal blocks CV \n")
        
        cv.res <- data.frame(matrix(nrow = 5, ncol = 4))
        colnames(cv.res) <- c("AUC", "Spec", "Sens", "TSS")
        
        gbm.m.lat <- list()
        
        for (k in 1:5) {
                
                gbm.m.lat[[k]] <- gbm.fixed(
                                sp.data[[species]]$
                                        data[sp.data[[species]]$
                                                     lat_blocks != k,],
                                gbm.x = 2:8, gbm.y = 1,
                                tree.complexity = gbm.m$interaction.depth,
                                learning.rate = gbm.m$shrinkage,
                                n.trees = gbm.m$n.trees,
                                site.weights = sp.data[[species]]$
                                        weigths[sp.data[[species]]$
                                                        lat_blocks != k]) 
                
                test.data <- sp.data[[species]]$data[sp.data[[species]]$
                                                             lat_blocks == k,]
                
                pred <- predict(gbm.m.lat[[k]], test.data)
                
                
                td <- data.frame(id = 1,
                                 obs = test.data$dsp,
                                 pred = pred)
                
                thresh <- optimal.thresholds(td,
                                             opt.methods = "MaxSens+Spec")$pred
                
                cmat <- cmx(td, threshold = thresh)
                
                cv.res[k,1] <- PresenceAbsence::auc(td)[,1]
                
                cv.res[k,2] <- PresenceAbsence::specificity(cmat)[,1]
                
                cv.res[k,3] <- PresenceAbsence::sensitivity(cmat)[,1]
                
                cv.res[k,4] <- cv.res[k,2] + cv.res[k,3] - 1
                
                rm(test.data, pred, td, thresh, cmat)

        }
        
        print(round(cv.res, 2))
        
        cat("Mean for latitudinal blocks \n AUC  Spec Sens TSS \n", 
            green(round(apply(cv.res, 2, mean),2)), "\n")
        
        
        pred.models <- function(model, m.code){
                
                cur.pred <- predict(env, model)
                writeRaster(cur.pred,
                            paste0(fin.f, "projections/", m.code, "_current.tif"))
                
                proj.list <- list(r12, r24, r37, r85)
                scenarios <- c("ssp126", "ssp245", "ssp370", "ssp585")
                
                for (i in 1:4) {
                        fut.pred <- predict(proj.list[[i]], model)
                        writeRaster(fut.pred,
                                    paste0(fin.f, "projections/",
                                           m.code, "_", scenarios[i], ".tif"))
                }
                
        }
        
        
        cat(blue("## Generating projections ## \n"))
        
        pred.models(gbm.m, "GBM")
        
        cat(blue("## Saving models results ## \n"))
        
        saveRDS(gbm.m, file = paste0(fin.f, "models/",
                                     species, "_gbm_model.RDS"))
        
        write.csv(gbm.m$cv.roc.matrix,
                  paste0(fin.f, "eval/", species, "_gbm_cv.csv"))
        
        write.csv(cv.res,
                  paste0(fin.f, "eval/", species, "_gbm_latcv.csv"))
        
        
        
        
        
        #### Random Forest ----
        cat(green("===== Fitting RF model ===== \n"))
        
        sumfun <- function(data, lev = NULL, model = NULL){
                td <- data.frame(
                        id = 1,
                        obs = data$obs,
                        pred = data$pred
                )
                
                #out <- PresenceAbsence::auc(td)$AUC
                
                rocObject <- try(PresenceAbsence::auc(td)$AUC)
                out <- if (inherits(rocObject, "try-error")){NA
                        }else{rocObject}
                
                names(out) <- "AUC"
                
                out
        }
        
        # Create a list to be used in the tunning
        indexlist <- list(seq(1:nrow(sp.data[[species]]$data))
                                                [sp.data[[species]]$blocks != 1],
                          seq(1:nrow(sp.data[[species]]$data))
                                                [sp.data[[species]]$blocks != 2],
                          seq(1:nrow(sp.data[[species]]$data))
                                                [sp.data[[species]]$blocks != 3],
                          seq(1:nrow(sp.data[[species]]$data))
                                                [sp.data[[species]]$blocks != 4],
                          seq(1:nrow(sp.data[[species]]$data))
                                                [sp.data[[species]]$blocks != 5])
        
        fitControl <- trainControl(## 10-fold CV
                method = "cv",
                number = 5,
                index = indexlist,
                classProbs = TRUE,
                summaryFunction = sumfun)
        
        rf1 <- train(dsp ~ ., data = sp.data$lyva$data, 
                     method = "rf", 
                     trControl = fitControl,
                     ## This last option is actually one
                     ## for gbm() that passes through
                     verbose = T,
                     metric = "AUC")
        
        mtry <- rf1$bestTune$mtry
        
        
        # rf.train <- expand.grid(ntree = seq(500, 800, 100),
        #                         mtry = c(2, 4, 7))
        # 
        # 
        # rfTuning <- function(index, data, ntree, mtry){
        #         
        #         rmse <- rep(NA, 5)
        #         
        #         for (i in 1:5) {
        #                 
        #                 train <- data[index == i, ]
        #                 test <- data[-index != i, ]
        #                 
        #                 rf <- randomForest(x = train[, 2:ncol(train)],
        #                                    y = train[, 1],
        #                                    ntree = ntree,
        #                                    mtry = mtry)
        #                 
        #                 pred <- predict(rf, newdata = test[, 2:8])
        #                 
        #                 
        #         }
        #         
        #         final.rmse <- mean(rmse)
        #         
        #         final.rmse
        # }
        # 
        # 
        # results <- rep(NA, nrow(rf.train))
        # 
        # 
        # for (j in 1:nrow(rf.train)) {
        #         results[j] <- rfTuning(sp.data[[species]]$blocks,
        #                                sp.data[[species]]$data,
        #                                ntree = rf.train[j,1],
        #                                mtry = rf.train[j,2])
        # }
        # 
        # mtry <- rf.train[which(results == max(results)),2]
        
        
        rf.crossval <- function(data, index, mtry){
                
                cv.res <- data.frame(matrix(nrow = 5, ncol = 4))
                colnames(cv.res) <- c("AUC", "Spec", "Sens", "TSS")
                
                for (i in 1:5) {
                        
                        train <- data[index != i, ]
                        test <- data[index == i, ]
                        
                        rf <- randomForest(x = train[, 2:ncol(train)],
                                           y = train[, 1],
                                           mtry = mtry)
                        
                        pred <- predict(rf, newdata = test[, 2:8])
                        
                        td <- data.frame(id = 1,
                                         obs = test[,1],
                                         pred = pred)
                        
                        thresh <- optimal.thresholds(td,
                                                     opt.methods = "MaxSens+Spec")$pred
                        
                        cmat <- cmx(td, threshold = thresh)
                        
                        cv.res[i,1] <- PresenceAbsence::auc(td)[,1]
                
                        cv.res[i,2] <- PresenceAbsence::specificity(cmat)[,1]
                
                        cv.res[i,3] <- PresenceAbsence::sensitivity(cmat)[,1]
                
                        cv.res[i,4] <- cv.res[i,2] + cv.res[i,3] - 1
                        
                }
                
                cv.res
                
                
        }
        
        rf.cv.res <- rf.crossval(sp.data[[species]]$data,
                                sp.data[[species]]$blocks,
                                mtry = mtry)
        
        print(round(rf.cv.res, 2))
        
        cat("Mean for blocks \n AUC  Spec Sens TSS \n", 
            green(round(apply(rf.cv.res, 2, mean),2)), "\n")
        
        
        rf.lat.res <- rf.crossval(sp.data[[species]]$data,
                                  sp.data[[species]]$lat_blocks,
                                  mtry = mtry)
        
        print(round(rf.lat.res, 2))
        
        cat("Mean for latitudinal blocks \n AUC  Spec Sens TSS \n", 
            green(round(apply(rf.lat.res, 2, mean),2)), "\n")
        
        # Generates final full model
        rf.final <- randomForest(x = sp.data[[species]]$data[,
                                        2:ncol(sp.data[[species]]$data)],
                                 y = sp.data[[species]]$data[, 1],
                                 mtry = mtry,
                                 strata = factor(c(0,1)))
        
        
        cat(blue("## Generating projections ## \n"))
        
        pred.models(rf.final, "RF")
        
        cat(blue("## Saving models results ## \n"))
        
        saveRDS(rf.final, file = paste0(fin.f, "models/",
                                     species, "_rf_model.RDS"))
        
        write.csv(rf.cv.res,
                  paste0(fin.f, "eval/", species, "_rf_cv.csv"))
        
        write.csv(rf.lat.res,
                  paste0(fin.f, "eval/", species, "_rf_latcv.csv"))
        
        
        
        
        
        #### GLM ----
        cat(green("===== Fitting GLM model ===== \n"))
        
        # creates the formula
        glm.formula <- dsp ~ 1 + BO21_tempmean_ss + I(BO21_tempmean_ss^2) + BO21_salinitymean_ss + 
                BO21_lightbotmax_bdmean + I(BO21_lightbotmax_bdmean^2) + 
                BO_ph + I(BO_ph^2) + BO21_chlomean_ss + I(BO21_chlomean_ss^2) + 
                BO21_silicatemean_ss + I(BO21_silicatemean_ss^2) + BO21_dissoxmean_ss + 
                I(BO21_dissoxmean_ss^2)
        
        glm.m <- glm(glm.formula, binomial(link = 'logit'),
                    data = sp.data[[species]]$data,
                    weights = sp.data[[species]]$weigths,
                    mustart = rep(0.5, nrow(sp.data[[species]]$data)),
                    control = glm.control(epsilon = 1e-08, maxit = 50, trace = FALSE))
        
        glm.best <- stepAIC(glm.m, glm.formula, 
                            data = sp.data[[species]]$data,
                            direction = "both", trace = FALSE, 
                            weights = sp.data[[species]]$weigths, 
                            steps = 10000,
                            mustart = rep(0.5, nrow(sp.data[[species]]$data)))
        
        
        glm.cv <- function(best.formula, index, data, wgt){
                
                cv.res <- data.frame(matrix(nrow = 5, ncol = 4))
                colnames(cv.res) <- c("AUC", "Spec", "Sens", "TSS")
                
                for (i in 1:5) {
                        
                        train <- data[index != i, ]
                        test <- data[index == i, ]

                        form <- as.character(best.formula)
                        
                        glm.temp <- glm(formula = as.formula(paste(form[2], form[1], form[3])),
                                        family = binomial(link = 'logit'),
                                        data = train,
                                        weights = wgt[index != i],
                                        mustart = rep(0.5, nrow(train)),
                                        control = glm.control(epsilon = 1e-08,
                                                              maxit = 50,
                                                              trace = FALSE))
                        
                        pred <- predict(glm.temp, newdata = test[, 2:8])
                        
                        td <- data.frame(id = 1,
                                         obs = test[,1],
                                         pred = pred)
                        
                        thresh <- optimal.thresholds(td,
                                        opt.methods = "MaxSens+Spec")$pred
                        
                        cmat <- cmx(td, threshold = thresh)
                        
                        cv.res[i,1] <- PresenceAbsence::auc(td)[,1]
                        
                        cv.res[i,2] <- PresenceAbsence::specificity(cmat)[,1]
                        
                        cv.res[i,3] <- PresenceAbsence::sensitivity(cmat)[,1]
                        
                        cv.res[i,4] <- cv.res[i,2] + cv.res[i,3] - 1
                        
                }
                
                cv.res
                
        }
        
        glm.block.cv <- glm.cv(best.formula = glm.best$formula,
                                index = sp.data[[species]]$blocks,
                                data = sp.data[[species]]$data,
                                wgt = sp.data[[species]]$weigths)
        
        print(round(glm.block.cv, 2))
        
        cat("Mean for blocks \n AUC  Spec Sens TSS \n", 
            green(round(apply(glm.block.cv, 2, mean),2)), "\n")
        
        glm.lat.cv <- glm.cv(best.formula = glm.best$formula,
                               index = sp.data[[species]]$lat_blocks,
                               data = sp.data[[species]]$data,
                               wgt = sp.data[[species]]$weigths)
        
        print(round(glm.lat.cv, 2))
        
        cat("Mean for latitudinal blocks \n AUC  Spec Sens TSS \n", 
            green(round(apply(glm.lat.cv, 2, mean),2)), "\n")
        
        
        
        cat(blue("## Generating projections ## \n"))
        
        pred.models(glm.best, "GLM")
        
        cat(blue("## Saving models results ## \n"))
        
        saveRDS(glm.best, file = paste0(fin.f, "models/",
                                        species, "_glm_model.RDS"))
        
        write.csv(glm.block.cv,
                  paste0(fin.f, "eval/", species, "_glm_cv.csv"))
        
        write.csv(glm.lat.cv,
                  paste0(fin.f, "eval/", species, "_glm_latcv.csv"))
        
        
        
        #### GAM ----
        #cat(green("===== Fitting GAM model ===== \n"))
        
        write(paste("Modeling of", species, "finished at:", date()),
              file = paste0(out.f, "log_file.txt"), append = T)
        
}


# Closes log file
write(paste("Modeling procedure finished at:", date()),
      file = paste0(out.f, "log_file.txt"), append = T)


itens <- ls()
itens <- itens[itens != "sp.data"]
rm(list = c(itens, "itens"))

##### Threshold by maxTSS, MTP and p10
gbm.m <- readRDS("modeling/regular/eclu/models/eclu_gbm_model.RDS")

glm.m <- readRDS("modeling/regular/eclu/models/eclu_glm_model.RDS")

rf.m <- readRDS("modeling/regular/eclu/models/eclu_rf_model.RDS")

gbm.r <- raster("modeling/regular/eclu/projections/GBM_current.tif")
glm.r <- raster("modeling/regular/eclu/projections/GLM_current.tif")
rf.r <- raster("modeling/regular/eclu/projections/RF_current.tif")

gbm.mtp <- sdm_threshold(gbm.r, sp.data$eclu$points[sp.data$eclu$pa == 1,])
gbm.p10 <- sdm_threshold(gbm.r, sp.data$eclu$points[sp.data$eclu$pa == 1,],
                         type = "p10")

glm.mtp <- sdm_threshold(glm.r, sp.data$eclu$points[sp.data$eclu$pa == 1,])
glm.p10 <- sdm_threshold(glm.r, sp.data$eclu$points[sp.data$eclu$pa == 1,],
                         type = "p10")

rf.mtp <- sdm_threshold(rf.r, sp.data$eclu$points[sp.data$eclu$pa == 1,])
rf.p10 <- sdm_threshold(rf.r, sp.data$eclu$points[sp.data$eclu$pa == 1,],
                         type = "p10")

thr.tss <- function(model, data, r){
  
  #pred <- predict(model, data[,2:ncol(data)])
  
  thr <- threshold(evaluate(data[data[,1] == 1,],
                            data[data[,1] == 0,],
                            model))
  
  if (thr$spec_sens < 0) {
    r[r >= thr$spec_sens] <- 1
    r[r < thr$spec_sens] <- 0
  } else {
    r[r < thr$spec_sens] <- 0
    r[r >= thr$spec_sens] <- 1
  }
  
  r
  
}

gbm.tss <- thr.tss(gbm.m, sp.data$eclu$data, gbm.r)
glm.tss <- thr.tss(glm.m, sp.data$eclu$data, glm.r)
rf.tss  <- thr.tss(rf.m, sp.data$eclu$data, rf.r )


par(mfrow = c(3, 3))
lapply(list(gbm.mtp, gbm.p10, gbm.tss, glm.mtp, glm.p10, glm.tss, rf.mtp, rf.p10, rf.tss), plot)

st <- stack(gbm.mtp, gbm.p10, gbm.tss, glm.mtp, glm.p10, glm.tss, rf.mtp, rf.p10, rf.tss)

#sdmvis::sdm_leaflet(st, pts = cbind(sp.data$eclu$points, sp.data$eclu$pa))
