#### Modelling of Herbivorous invertebrates of coral reefs ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2021

### Modeling functions ###

# Get weigthed response function
wt <- function(resp, prev = 0.5){
        
        pres <- sum(resp)
        
        absc <- length(resp) - pres
        
        wgth <- rep(1, length(resp))
        
        wgth[which(resp > 0)] <- (prev * absc)/(pres * (1 - prev))
        
        wgth <- round(wgth)
        
        wgth
}

# Function to train GBM auto
gbm.auto <- function(sp.data, weigths, blocks){
        
        gbm.m <- gbm.step(sp.data,
                          gbm.x = 2:length(sp.data), gbm.y = 1,
                          fold.vector = blocks, n.folds = 5,
                          site.weights = weigths,
                          keep.fold.fit = T, plot.main = F)
        
        if (gbm.m$n.trees < 1000) {
                cat("Reducing learning rate.")
                
                lr <- gbm.m$shrinkage
                
                tn <- gbm.m$n.trees
                
                while (tn < 1000) {
                        lr <- lr/2
                        
                        if (lr < 0.001) {
                                cat("Learning rate reached less than 0.001.
                                    Fitting a final model with lr = 0.001")
                                
                                lr <- 0.001
                                
                                gbm.m <- gbm.step(sp.data,
                                                  gbm.x = 2:length(sp.data), gbm.y = 1,
                                                  fold.vector = blocks, n.folds = 5,
                                                  site.weights = weigths,
                                                  learning.rate = lr,
                                                  keep.fold.fit = T,
                                                  plot.main = F)
                                
                                tn <- 1000
                                
                        } else {
                                
                                gbm.m <- gbm.step(sp.data,
                                                  gbm.x = 2:length(sp.data), gbm.y = 1,
                                                  fold.vector = blocks, n.folds = 5,
                                                  site.weights = weigths,
                                                  learning.rate = lr,
                                                  keep.fold.fit = T,
                                                  plot.main = F)
                                
                                tn <- gbm.m$n.trees
                        }
                }
        }
        
        return(gbm.m)
}

# Function to get predictions
pred.sdm <- function(rast, model, thresh = NULL,
                     outfold = NULL, filenam = NULL, code = NULL){
        
        for (j in 1:length(rast)) {
                
                pred <- predict(rast[[j]], model, type = "response")
                
                if (j == 1) {
                        outstack <- pred
                } else {
                        outstack <- stack(outstack, pred)
                }
        }
        
        names(outstack) <- names(rast)
        
        if (!is.null(outfold)) {
                writeRaster(outstack,
                            paste0(paste0(outfold, filenam, "_", 
                                          tolower(code), "_"),
                                   names(outstack), ".tif"),
                            bylayer = T)   
        }
        
        if (!is.null(thresh)) {
                outstack <- calc(outstack, function(x){
                        x[x < thresh] <- 0
                        x[x >= thresh] <- 1
                        x
                })
                
                names(outstack) <- names(rast)
                
                if (!is.null(outfold)) {
                        writeRaster(outstack,
                                    paste0(paste0(outfold, filenam, "_", 
                                                  tolower(code), "_"),
                                           names(outstack), "_bin.tif"),
                                    bylayer = T)   
                }
        }
        
        outstack
}

# Metrics
model.metrics <- function(original, predicted, raster = F){
        
        if (isTRUE(raster)) {
                td <- data.frame(id = 1,
                                 obs = rasterToPoints(original)[,3],
                                 pred = rasterToPoints(predicted)[,3])
                
                cmat <- cmx(td, threshold = 0.5)
        } else {
                td <- data.frame(id = 1,
                                 obs = original,
                                 pred = predicted)
                
                thresh <- optimal.thresholds(td,
                                             opt.methods = "MaxSens+Spec")$pred
                
                cmat <- cmx(td, threshold = thresh)
        }
        
        au <- PresenceAbsence::auc(td)[,1]
        
        spec <- PresenceAbsence::specificity(cmat)[,1]
        
        sens <- PresenceAbsence::sensitivity(cmat)[,1]
        
        tss <- (spec + sens) - 1
        
        return(c(au, spec, sens, tss))
}

# Get threshold
get.thresh <- function(observed, predicted){
        td <- data.frame(id = 1,
                         obs = observed,
                         pred = predicted)
        
        thr <- optimal.thresholds(td,
                                      opt.methods = "MaxSens+Spec")$pred
        thr
}

# Variable importance
var.importance <- function(data, cols = 2:length(data), model, n = 5,
                           get.mean = F){
        
        results <- matrix(nrow = n, ncol = length(cols))
        
        colnames(results) <- colnames(data[,cols])
        
        for (i in 1:n) {
                
                for (k in cols) {
                        
                        x <- data
                        
                        x[,k] <- x[sample.int(nrow(data)), k]
                        
                        pred <- predict(model, data[, cols], type = "response")
                        
                        pred.new <- predict(model, x[, cols], type = "response")
                        
                        if (min(cols) == 1) {
                                results[i, k] <- (1 - cor(pred, pred.new))
                        } else{
                                results[i, (k-1)] <- (1 - cor(pred, pred.new))
                        }
                        
                        
                }
        }
        
        if (isTRUE(get.mean)) {
                round(apply(results, 2, mean),3)
        }else{
                results
        }
}
