#### Modelling of Herbivorous invertebrates of coral reefs ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2021

### Modeling - Sea-urchins ###
### Functions ###

### Functions to perform the training of models, cross-validation and
### latitudinal cross-validation

### GBM ----
gbm.modeling <- function(data, blocks, lat.blocks, weights, ...){
        
        # Train the initial model and find the optimal one
        gbm.m <- gbm.step(data, gbm.x = 2:ncol(data), gbm.y = 1,
                          fold.vector = blocks, n.folds = 5,
                          site.weights = weights,
                          keep.fold.fit = T, ...)
        
        # Get other stats of the blocks CV
        cv.res.init <- data.frame(matrix(nrow = 5, ncol = 5))
        colnames(cv.res.init) <- c("AUC", "Spec", "Sens", "TSS", "Boyce")
        
        for (k in 1:5) {
                
                pred <- gbm.m$fold.fit[blocks == k]
                
                pred <- exp(pred)/(1 + exp(pred))
                
                obs <- data[,1][blocks == k]
                
                td <- data.frame(id = 1,
                                 obs = obs,
                                 pred = pred)
                
                thresh <- optimal.thresholds(td,
                                             opt.methods = "MaxSens+Spec")$pred
                
                cmat <- cmx(td, threshold = thresh)
                
                cv.res.init[k,1] <- PresenceAbsence::auc(td)[,1]
                
                cv.res.init[k,2] <- PresenceAbsence::specificity(cmat)[,1]
                
                cv.res.init[k,3] <- PresenceAbsence::sensitivity(cmat)[,1]
                
                cv.res.init[k,4] <- (cv.res.init[k,2] + cv.res.init[k,3]) - 1
                
                cv.res.init[k,5] <- ecospat::ecospat.boyce(pred,
                                                      pred[obs == 1],
                                                      PEplot = F)$Spearman.cor
                
                
        }
        
        # Cross-validation for latitudinal blocks
        cv.res <- data.frame(matrix(nrow = 5, ncol = 5))
        colnames(cv.res) <- c("AUC", "Spec", "Sens", "TSS", "Boyce")
        
        gbm.m.lat <- list()
        
        for (k in 1:5) {
                
                gbm.m.lat[[k]] <- gbm.fixed(
                        data[lat.blocks != k,],
                        gbm.x = 2:ncol(data), gbm.y = 1,
                        tree.complexity = gbm.m$interaction.depth,
                        learning.rate = gbm.m$shrinkage,
                        n.trees = gbm.m$n.trees,
                        site.weights = weights[lat.blocks != k]) 
                
                test.data <- data[lat.blocks == k,]
                
                pred <- predict(gbm.m.lat[[k]], test.data[,2:ncol(data)],
                                type = "response")
                
                
                td <- data.frame(id = 1,
                                 obs = test.data[,1],
                                 pred = pred)
                
                thresh <- optimal.thresholds(td,
                                             opt.methods = "MaxSens+Spec")$pred
                
                cmat <- cmx(td, threshold = thresh)
                
                cv.res[k,1] <- PresenceAbsence::auc(td)[,1]
                
                cv.res[k,2] <- PresenceAbsence::specificity(cmat)[,1]
                
                cv.res[k,3] <- PresenceAbsence::sensitivity(cmat)[,1]
                
                cv.res[k,4] <- (cv.res[k,2] + cv.res[k,3]) - 1
                
                cv.res[k,5] <- ecospat::ecospat.boyce(pred,
                                                      pred[test.data[,1] == 1],
                                                      PEplot = F)$Spearman.cor
                
                rm(test.data, pred, td, thresh, cmat)
                
        }
        
        final.result <- list(gbm.m, cv.res.init, cv.res)
        names(final.result) <- c("model", "block_cv", "latitudinal_cv")
        
        # Get final model results
        td <- data.frame(id = 1,
                         obs = data[,1],
                         pred = gbm.m$fitted)
        
        thresh <- optimal.thresholds(td,
                                     opt.methods = "MaxSens+Spec")$pred
        
        final.result$tss_thresh <- thresh
        
        final.result$boyce <- ecospat::ecospat.boyce(gbm.m$fitted,
                                                     gbm.m$fitted[data[,1]== 1],
                                                     PEplot = F)
        
        final.result
        
}



### RF ----
rf.modeling <- function(data, blocks, lat.blocks, ...){
        
        # Create a list to be used in the tunning
        indexlist <- list(seq(1:nrow(data))[blocks != 1],
                          seq(1:nrow(data))[blocks != 2],
                          seq(1:nrow(data))[blocks != 3],
                          seq(1:nrow(data))[blocks != 4],
                          seq(1:nrow(data))[blocks != 5])
        
        fitControl <- trainControl(## 10-fold CV
                method = "cv",
                number = 5,
                index = indexlist)
        
        rf1 <- train(dsp ~ ., data = data, 
                     method = "rf", 
                     trControl = fitControl,
                     ## This last option is actually one
                     ## for gbm() that passes through
                     verbose = T)
        
        mtry <- rf1$bestTune$mtry
        
        
        rf.crossval <- function(data, index, mtry, ...){
                
                cv.res <- data.frame(matrix(nrow = 5, ncol = 5))
                colnames(cv.res) <- c("AUC", "Spec", "Sens", "TSS", "Boyce")
                
                for (k in 1:5) {
                        
                        train <- data[index != k, ]
                        test <- data[index == k, ]
                        
                        rf <- randomForest(x = train[, 2:ncol(train)],
                                           y = train[, 1],
                                           mtry = mtry,
                                           ...)
                        
                        pred <- predict(rf, newdata = test[, 2:ncol(test)],
                                        type = "response")
                        
                        td <- data.frame(id = 1,
                                         obs = test[,1],
                                         pred = pred)
                        
                        thresh <- optimal.thresholds(td,
                                            opt.methods = "MaxSens+Spec")$pred
                        
                        cmat <- cmx(td, threshold = thresh)
                        
                        cv.res[k,1] <- PresenceAbsence::auc(td)[,1]
                        
                        cv.res[k,2] <- PresenceAbsence::specificity(cmat)[,1]
                        
                        cv.res[k,3] <- PresenceAbsence::sensitivity(cmat)[,1]
                        
                        cv.res[k,4] <- (cv.res[k,2] + cv.res[k,3]) - 1
                        
                        cv.res[k,5] <- ecospat::ecospat.boyce(pred,
                                                  pred[test[,1] == 1],
                                                  PEplot = F)$Spearman.cor
                        
                        rm(td, thresh, cmat)
                        
                }
                
                cv.res
                
        }
        
        rf.block.cv <- rf.crossval(data, blocks, mtry, ...) 
        
        rf.lat.cv <- rf.crossval(data, lat.blocks, mtry, ...) 
        
        
        # Get final model
        rf.m <- randomForest(x = data[, 2:ncol(data)],
                             y = data[, 1],
                             mtry = mtry)
        
        
        final.result <- list(rf.m, rf.block.cv, rf.lat.cv)
        names(final.result) <- c("model", "block_cv", "latitudinal_cv")
        
        # Get final model results
        td <- data.frame(id = 1,
                         obs = data[,1],
                         pred = predict(rf.m, data[,2:ncol(data)],
                                        type = "response"))
        
        thresh <- optimal.thresholds(td,
                                     opt.methods = "MaxSens+Spec")$pred
        
        final.result$tss_thresh <- thresh
        
        final.result$boyce <- ecospat::ecospat.boyce(td$pred,
                                                     td$pred[data[,1]== 1],
                                                     PEplot = F)
        
        final.result
        
}



### GLM ----
glm.modeling <- function(start.formula, data, blocks, lat.blocks, weigths, ...){
        
        wt <- weigths
        
        glm.m <- glm(formula(start.formula), binomial(link = 'logit'),
                     data = data,
                     weights = weigths,
                     mustart = rep(0.5, nrow(data)),
                     control = glm.control(epsilon = 1e-08,
                                           maxit = 50,
                                           trace = FALSE))
        
        glm.best <- stepAIC(glm.m, formula(start.formula), 
                            data = data,
                            direction = "both", trace = FALSE, 
                            weights = weigths, 
                            steps = 10000,
                            mustart = rep(0.5, nrow(data)))
        
        
        
        glm.crossval <- function(best.formula, index, dat, wgt, ...){
                
                cv.res <- data.frame(matrix(nrow = 5, ncol = 5))
                colnames(cv.res) <- c("AUC", "Spec", "Sens", "TSS", "Boyce")
                
                for (k in 1:5) {
                        
                        train <- dat[index != k, ]
                        test <- dat[index == k, ]
                        
                        form <- as.character(best.formula)
                        
                        glm.temp <- glm(formula = as.formula(paste(form[2],
                                                                   form[1],
                                                                   form[3])),
                                        family = binomial(link = 'logit'),
                                        data = train,
                                        weights = wgt[index != k],
                                        ...)
                        
                        pred <- predict(glm.temp, newdata = test[, 2:ncol(test)],
                                        type = "response")
                        
                        td <- data.frame(id = 1,
                                         obs = test[,1],
                                         pred = pred)
                        
                        thresh <- optimal.thresholds(td,
                                        opt.methods = "MaxSens+Spec")$pred
                        
                        cmat <- cmx(td, threshold = thresh)
                        
                        cv.res[k,1] <- PresenceAbsence::auc(td)[,1]
                        
                        cv.res[k,2] <- PresenceAbsence::specificity(cmat)[,1]
                        
                        cv.res[k,3] <- PresenceAbsence::sensitivity(cmat)[,1]
                        
                        cv.res[k,4] <- (cv.res[k,2] + cv.res[k,3]) - 1
                        
                        cv.res[k,5] <- ecospat::ecospat.boyce(pred,
                                                pred[test[,1] == 1],
                                                PEplot = F)$Spearman.cor
                        
                        rm(td, thresh, cmat)
                        
                }
                
                cv.res
                
        }
        
        
        glm.block.cv <- glm.crossval(glm.best$formula,
                                     blocks, data, wgt = wt)
        
        
        glm.lat.cv <- glm.crossval(glm.best$formula,
                                     lat.blocks, data, wgt = wt)
        
        
        
        final.result <- list(glm.best, glm.block.cv, glm.lat.cv)
        names(final.result) <- c("model", "block_cv", "latitudinal_cv")
        
        # Get final model results
        td <- data.frame(id = 1,
                         obs = data[,1],
                         pred = predict(glm.best, data[,2:ncol(data)],
                                        type = "response"))
        
        thresh <- optimal.thresholds(td,
                                     opt.methods = "MaxSens+Spec")$pred
        
        final.result$tss_thresh <- thresh
        
        final.result$boyce <- ecospat::ecospat.boyce(td$pred,
                                                     td$pred[data[,1]== 1],
                                                     PEplot = F)
        
        final.result
        
}



### GAM ----
gam.modeling <- function(start.formula, data, blocks, lat.blocks, weigths, ...){
        
        gam.formula <- formula(start.formula)
        
        gam.m <- mgcv::gam(gam.formula,
                          data = data,
                          family = binomial(link = 'logit'), 
                          weights = weigths,
                          select=TRUE, method='REML', ...)
        
        
        gam.crossval <- function(best.formula, index, data, wgt, ...){
                
                cv.res <- data.frame(matrix(nrow = 5, ncol = 5))
                colnames(cv.res) <- c("AUC", "Spec", "Sens", "TSS", "Boyce")
                
                for (k in 1:5) {
                        
                        train <- data[index != k, ]
                        test <- data[index == k, ]
                        
                        #form <- as.character(best.formula)
                        
                        gam.temp <- mgcv::gam(best.formula,
                                              data = train,
                                              family = binomial(link = 'logit'))
                        
                        pred <- predict(gam.temp, newdata = test[, 2:ncol(test)],
                                        type = "response")
                        
                        td <- data.frame(id = 1,
                                         obs = test[,1],
                                         pred = pred)
                        
                        thresh <- optimal.thresholds(td,
                                                     opt.methods = "MaxSens+Spec")$pred
                        
                        cmat <- cmx(td, threshold = thresh)
                        
                        cv.res[k,1] <- PresenceAbsence::auc(td)[,1]
                        
                        cv.res[k,2] <- PresenceAbsence::specificity(cmat)[,1]
                        
                        cv.res[k,3] <- PresenceAbsence::sensitivity(cmat)[,1]
                        
                        cv.res[k,4] <- (cv.res[k,2] + cv.res[k,3]) - 1
                        
                        cv.res[k,5] <- ecospat::ecospat.boyce(pred,
                                                              pred[test[,1] == 1],
                                                              PEplot = F)$Spearman.cor
                        
                        rm(td, thresh, cmat)
                        
                }
                
                cv.res
                
        }
        
        
        gam.block.cv <- gam.crossval(gam.m$formula,
                                     blocks, data, weights)
        
        
        gam.lat.cv <- gam.crossval(gam.m$formula,
                                   lat.blocks, data, weights)
        
        
        
        final.result <- list(gam.m, gam.block.cv, gam.lat.cv)
        names(final.result) <- c("model", "block_cv", "latitudinal_cv")
        
        # Get final model results
        td <- data.frame(id = 1,
                         obs = data[,1],
                         pred = predict(gam.m, data[,2:ncol(data)],
                                        type = "response"))
        
        thresh <- optimal.thresholds(td,
                                     opt.methods = "MaxSens+Spec")$pred
        
        final.result$tss_thresh <- thresh
        
        final.result$boyce <- ecospat::ecospat.boyce(td$pred,
                                                     td$pred[data[,1]== 1],
                                                     PEplot = F)
        
        final.result
        
}



# A melhorar
#trocar método de obter TSS, thresholds, etc >>> usar o do DISMO?


## Prediction function ----
pred.sdm <- function(rast, model, thresh = NULL){
        
        pred <- predict(rast, model, type = "response")
        
        if (!is.null(thresh)) {
                pred <- calc(pred, function(x){
                        x[x < thresh] <- 0
                        x[x >= thresh] <- 1
                        x
                })
        }
        
        pred
}


## Save model stuff ----
save.mod <- function(results, folder, model.name, sp.name){
        
        saveRDS(res, file = paste0(folder, "/models/",
                                   model.name, "_", sp.name,  ".rds"))
        
        write.csv(results$block_cv,
                  paste0(folder, "/eval/",
                         model.name, "_", sp.name, "_blockcv.csv"))
        
        write.csv(results$latitudinal_cv,
                  paste0(folder, "/eval/",
                         model.name, "_", sp.name, "_latcv.csv"))
        
        write.csv(results$varimp,
                  paste0(folder, "/eval/",
                         model.name, "_", sp.name, "_varimp.csv"))
        
        df <- data.frame(boyce = results$boyce$Spearman.cor,
                         tss_thresh = results$tss_thresh)
        
        write.csv(df,
                  paste0(folder, "/eval/",
                         model.name, "_", sp.name, "_boyce_tss.csv"))
        
        
}
