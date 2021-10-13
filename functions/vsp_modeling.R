#### Modelling of Herbivorous invertebrates of coral reefs ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2021

### Virtual species modeling ###
### Modeling functions ###

# Function to get response curves
resp.curve <- function(model, mname, env.dat){
        
        results <- data.frame(matrix(nrow = 50, ncol = nrow(env.dat)))
        
        all.data <- results
        
        names(all.data) <- env.dat$var
        
        for (j in 1:nrow(env.dat)) {
                all.data[,j] <- env.dat$mean[j]
        }
        
        for (z in 1:nrow(env.dat)) {
                
                df <- all.data
                
                df[,z] <- seq(env.dat$min[z], env.dat$max[z], length = 50)
                
                results[,z] <- predict(model, df, type = "response")
        }
        
        names(results) <- env.dat$var
        
        results$filename = mname
        
        results
}

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
                          keep.fold.fit = T)
        
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
                                                  keep.fold.fit = T)
                                
                                tn <- 1000
                                
                        } else {
                                
                                gbm.m <- gbm.step(sp.data,
                                                  gbm.x = 2:length(sp.data), gbm.y = 1,
                                                  fold.vector = blocks, n.folds = 5,
                                                  site.weights = weigths,
                                                  learning.rate = lr,
                                                  keep.fold.fit = T)
                                
                                tn <- gbm.m$n.trees
                        }
                }
        }
        
        return(gbm.m)
}

# Function to get predictions
pred.sdm <- function(rast, model, thresh = NULL){
        
        for (j in 1:length(rast)) {
                
                pred <- predict(rast[[j]], model, type = "response")
                
                if (!is.null(thresh)) {
                        pred <- calc(pred, function(x){
                                x[x < thresh] <- 0
                                x[x >= thresh] <- 1
                                x
                        })
                }
                
                if (j == 1) {
                        outstack <- pred
                } else {
                        outstack <- stack(outstack, pred)
                }
        }
        
        names(outstack) <- names(rast)
        
        outstack
}

# Jaccard index
jaccard.eval <- function(rast1, rast2){
        combination <- rast1 + rast2
        intersection <- combination == 2
        union <- combination >= 1
        return(sum(rasterToPoints(intersection)[,3])/
                       sum(rasterToPoints(union)[,3]))
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
        
        # if (isTRUE(raster)) {
        #         return(c(au, spec, sens, tss))     
        # } else{
        #         boyce <- ecospat::ecospat.boyce(predicted,
        #                                               predicted[original == 1],
        #                                               PEplot = F)$Spearman.cor
        #         
        #         return(c(au, spec, sens, tss, boyce))
        # }
        
        return(c(au, spec, sens, tss))
}

# Compare models
compare.models <- function(predicted, orig){

        jacc <- jaccard.eval(orig, predicted)
        
        im <- nicheOverlap(orig, predicted, stat = "I")
        dm <- nicheOverlap(orig, predicted, stat = "D")
        
        met <- model.metrics(orig, predicted, raster = T)
        
        return(c(jacc, im, dm, met))
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