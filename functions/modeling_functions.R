#### Modelling of Herbivorous invertebrates of coral reefs ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2021

### Modeling - Sea-urchins ###
### Functions ###

### Functions to perform the cross-validation, predictions and save model data

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