#### Occurrence Thinner function, based on Verbruggen et. al work
# Adapted by Silas C. Principe
        
occthinner <- function(pts,          # Input coordinates (df with lon/lat)
                       kdensity,     # Kernel density raster
                       nr = 10,      # Number of replicates
                       t1 = 0,       # Threshold 1 (default to 0)
                       t2 = 1,       # Threshold 2 (default to 1)
                       cols = NULL,  # The name or index of lon/lat cols
                       save = FALSE, # Should files be saved at disk?
                       outf = NULL,  # In case files should be saved, where?
                       verbose = TRUE# Should info be printed?
                       ){
        
        # Set lon/lat col names
        if (is.null(cols)) {
                cols <- c("longitude", "latitude")
        }
        
        # Normalize kernel density raster
        kdensity <- setMinMax(kdensity)
        
        kdensity <- (kdensity - minValue(kdensity))/
                (maxValue(kdensity)-minValue(kdensity))
        
        # Creates a list to hold results
        thinned.list <- list()
        
        if (verbose) {
                cat("Running thinning procedure", nr, "times. \n")
        }
        
        # Run for nr replicates
        for (k in 1:nr) {
                
                # Creates an index vector
                sel.pts <- rep(NA, nrow(pts))
                
                # Run for each point
                for (i in 1:nrow(pts)) {
                        
                        # Sample a value between seq
                        ### Need confirmation with Verbruggen about the step(by)
                        #rd <- sample(seq(t1, t2, by = 0.1), 1)
                        rd <- runif(1, min = t1, max = t2)
                        
                        # Extract density value in the point
                        kde <- extract(kdensity,
                                       pts[i, cols])
                        
                        # If density value less than the sampled, point selected
                        # else, point will be excluded
                        sel.pts[i] <- ifelse(kde <= rd, 1, 0)
                        
                }
                
                # Get valid points
                thinned.list[[k]] <- pts[sel.pts == 1, ]
        }
        
        if (save) {
                if (verbose) {
                        cat("Saving files. \n")
                }
                
                if(is.null(outf)){
                        outf <- ""
                } else {
                        if (!dir.exists(outf)) {
                                dir.create(outf)
                        }
                        outf <- paste0(outf, "/")
                }
                
                for (z in 1:nr) {
                        write.csv(thinned.list[[z]],
                                  paste0(outf, "occ_thinned_", z, ".csv"))
                }
        }
        
        if (verbose) {
                cat("Number of rows in the original file:", nrow(pts), "\n")
                cat("Number of rows retained in each replicate: \n")
                cat(paste(1:nr, unlist(lapply(thinned.list, nrow)),
                          collapse = "\n", sep = "|"), "\n", sep = "")
        }
        
        return(thinned.list)
}


#### Example
# library(raster)
# 
# b.pts <- read.csv("occ_thinner_test.csv")
# 
# kde.raster <- raster("occ_thinner_kdetest.tif")
# 
# thinned <- occthinner(b.pts, kde.raster, t1 = 0.4, t2 = 1, cols = 1:2,
#                       save = T, outf = "pasta2/")
# 
# lapply(thinned, nrow)
# 
# par(mfrow = c(1,3))
# plot(est.raster)
# points(b.pts, pch = 20, cex = 0.5)
# plot(est.raster)
# points(thinned[[1]], pch = 20, cex = 0.5)
# plot(est.raster)
# points(thinned[[2]], pch = 20, cex = 0.5)
####