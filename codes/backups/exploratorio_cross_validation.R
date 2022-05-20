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
        
        res <- c()
        
        pb <- progress_bar$new(total = runs)
        
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
                
                pred <- predict(ms, pts, cv.f)
                b <- 1-exp(-pred$mean)
                
                res <- c(res, pred)
                
                pred.f <- predict(model, test)
                
                res.f <- c(res, pred.f)
                
                pb$tick()
        } 
        
        # Residuals
        resid <- pred.f-p
        # RMSE
        rmse <- sqrt(mean(resid^2))
        
}