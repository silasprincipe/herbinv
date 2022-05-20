sloo.cv <- function(model, forms, ips, pts, runs = 20, radius = 100){
        # @model = inlabru model
        # @forms = formula to the model
        # @ips = inlabru integration points
        # @pts = location points
        # @runs = number of runs of the LOOCV
        # @radius = radius of the buffer to remove points spatially
        
        library(progress)
        
        # Get points to leave out
        p <- sample(1:length(pts), runs)
        
        cv.f <- as.formula(paste0(
                "~ exp(",
                paste(if (is.null(names(
                        model$summary.random))) {
                        model$names.fixed
                } else{
                        c(model$names.fixed,
                          names(model$summary.random))
                }
                , collapse = "+"), ")"
        ))
        
        res <- data.frame(lq = rep(NA, runs), me = NA, uq = NA)
        
        cat("Starting loop. \n")
        
        for (i in 1:runs) {
                lop <- pts[p[i],]
                buf <- buffer(lop, radius)
                
                test <- intersect(pts, buf)
                train <- over(pts, buf)
                train <- pts[is.na(train),]
                
                ms <- lgcp(forms, train,
                           ips = ips,
                           domain = list(coordinates = mesh),
                           options = list(control.inla = list(int.strategy = "eb"),
                                          inla.mode = "classic"))
                
                Score_joint <- 
                        predict(ms, int_test_prod, 
                                ~ sum(dpois(count,
                                            lambda = weight * exp(Intercept + sal + sst + spatial),
                                            log = T)))
                
                res$lq[i] <- Score_joint$q0.025
                res$uq[i] <- Score_joint$q0.975
                res$me[i] <- Score_joint$mean
                
                #pb$tick()
        } 
        
        Score_joint <- 
                predict(model, int_test_prod, 
                        ~ sum(dpois(count,
                                    lambda = weight * exp(Intercept + sal + sst + spatial),
                                    log = T)))
        
        res <- rbind(cbind(res, cv = 1:nrow(res)),
                     cbind(lq = Score_joint$q0.025,
                           me = Score_joint$mean,
                           uq = Score_joint$q0.975,
                           cv = 0))
        
        # Residuals
        # resid <- pred.f-p
        # # RMSE
        # rmse <- sqrt(mean(resid^2))
        res
}

#pb <- progress_bar$new(total = 2)
cvres <- sloo.cv(m[[2]], forms[[1]], ips, pts, runs = 2, radius = 700)


# Loop through sightings. Find closest mesh locations and distance value.
int_ind <- rep(NA, length(pts))
int_test_prod <- ips
int_test_prod$count <- 0
# loop through the sightings
for(i in 1:length(pts)){
        int_ind[i] <- as.numeric(apply(gDistance(int_test_prod, pts[i,],
                                                 byid=T),
                                       1,FUN = function(x){which(x==min(x))}))
        # Update the count by 1
        int_test_prod$count[int_ind[i]] <- int_test_prod$count[int_ind[i]] + 1
}