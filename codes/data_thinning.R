#### Modelling of Herbivorous invertebrates of coral reefs ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2021

### Thinning of occurrence data ###

# Thinning is a procedure to reduce the effects of sampling bias
# Here we use the method from Verbruggen H. (2012) 
# [OccurrenceThinner version 1.04. github.com/hverbruggen/OccurrenceThinner]
# but using an adapted code in R

# We also apply a spatial thinning to reduce possible effects of spatial auto-
# correlation. We apply a filter so no points occur in the adjacent cells.

# Load functions/libraries ----
source("functions/occurrence_thinner.R")
library(KernSmooth)
library(raster)
library(spThin)
set.seed(2932)

# Creates a function to generate kernel rasters
genkernel <- function(sp.data){
        
        est <- bkde2D(sp.data, 
                      bandwidth=c(3,3), 
                      gridsize=c(840,1020),
                      range.x=list(c(-99,-29),c(-42.5,42.5)))
        
        est$fhat[est$fhat < 0.00001] <- 0
        est.raster <- raster(list(x=est$x1, y=est$x2, z=est$fhat))
        projection(est.raster) <- CRS("+init=epsg:4326")
        xmin(est.raster) <- -99
        xmax(est.raster) <- -29
        ymin(est.raster) <- -42.5
        ymax(est.raster) <- 42.5
        
        est.raster
}

# Run in loop for all three species ----
species <- c("lyva", "eclu", "trve")

for (i in species) {
        
        # Define folder
        folder <- paste0("data/", i, "/")
        
        # Load species data
        pts <- read.csv(paste0(folder, i, "_cell.csv"))
        
        # Kernel thinning ----
        # Generate kernel map
        spkn <- genkernel(pts[,1:2])
        
        # Thin occurrences
        pts.th <- occthinner(pts, spkn, t1 = 0.5, t2 = 1, cols = 1:2)
        
        # Select the one with the highest number of occurrences retained
        pts.th.sel <- pts.th[[which(unlist(lapply(pts.th, nrow)) ==
                                max(unlist(lapply(pts.th, nrow))))[1]]]
        
        lonlat <- names(pts.th.sel)[1:2]
        
        # Spatial thinning ----
        pts.spth <- thin(cbind(pts.th.sel, sp = "sp1"),
                         lat.col = lonlat[2], long.col = lonlat[1],
                         spec.col = "sp",
                         thin.par = 10, reps = 10, locs.thinned.list.return = T,
                         write.files = F, write.log.file = F)
        
        # Select the one with the highest number of occurrences retained
        pts.spth.sel <- pts.spth[[which(unlist(lapply(pts.spth, nrow)) ==
                                    max(unlist(lapply(pts.spth, nrow))))[1]]] 
        
        # Save files ----
        write.csv(pts.th.sel,  paste0(folder, i, "_knthin.csv"), row.names = F)
        write.csv(pts.spth.sel, paste0(folder, i, "_spthin.csv"), row.names = F)
}

# END