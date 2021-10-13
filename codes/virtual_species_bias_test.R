#### Modelling of Herbivorous invertebrates of coral reefs ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2021

### Testing sampling bias and autocorrelation removal methods ###
# In this file we will just generate the thinned datasets

# Load functions/libraries ----
source("functions/occurrence_thinner.R")
library(KernSmooth)
library(raster)
library(spThin)
set.seed(2932)

setwd("data/vsp/")

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

### Virtual Species 1 ----

# Load sampling bias files
vsp1.ln <- read.csv("vsp1_ptpolbiassamp_ln.csv")
vsp1.hn <- read.csv("vsp1_ptpolbiassamp_hn.csv")

vsp1.ln <- vsp1.ln[, 1:2]
vsp1.hn <- vsp1.hn[, 1:2]

# Generate kernel rasters
vsp1.ln.kn <- genkernel(vsp1.ln[,1:2])
vsp1.hn.kn <- genkernel(vsp1.hn[,1:2])

# Thin occurrences
vsp1.ln.thin <- occthinner(vsp1.ln, vsp1.ln.kn, t1 = 0.5, t2 = 1, cols = 1:2)
vsp1.hn.thin <- occthinner(vsp1.hn, vsp1.hn.kn, t1 = 0.5, t2 = 1, cols = 1:2)

vsp1.ln.tsel <- vsp1.ln.thin[[which(unlist(lapply(vsp1.ln.thin, nrow)) ==
                                   max(unlist(lapply(vsp1.ln.thin, nrow))))[1]]]
vsp1.hn.tsel <- vsp1.hn.thin[[which(unlist(lapply(vsp1.hn.thin, nrow)) ==
                                   max(unlist(lapply(vsp1.hn.thin, nrow))))[1]]]


# Spatial thinning
vsp1.ln.tselsp <- thin(cbind(vsp1.ln.tsel, sp = "vsp1"),
                       lat.col = "y", long.col = "x", spec.col = "sp",
                       thin.par = 10, reps = 10, locs.thinned.list.return = T,
                       write.files = F, write.log.file = F)

vsp1.hn.tselsp <- thin(cbind(vsp1.hn.tsel, sp = "vsp1"),
                       lat.col = "y", long.col = "x", spec.col = "sp",
                       thin.par = 10, reps = 10, locs.thinned.list.return = T,
                       write.files = F, write.log.file = F)


vsp1.ln.tselsp <- vsp1.ln.tselsp[[which(unlist(lapply(vsp1.ln.tselsp, nrow)) ==
                                 max(unlist(lapply(vsp1.ln.tselsp, nrow))))[1]]] 
vsp1.hn.tselsp <- vsp1.hn.tselsp[[which(unlist(lapply(vsp1.hn.tselsp, nrow)) ==
                                 max(unlist(lapply(vsp1.hn.tselsp, nrow))))[1]]] 

# Save files

write.csv(vsp1.ln.tsel, "vsp1_kernthin_ln.csv", row.names = F)
write.csv(vsp1.hn.tsel, "vsp1_kernthin_hn.csv", row.names = F)

write.csv(vsp1.ln.tselsp, "vsp1_spatialthin_ln.csv", row.names = F)
write.csv(vsp1.hn.tselsp, "vsp1_spatialthin_hn.csv", row.names = F)

### Virtual Species 2 ----

# Load sampling bias files
vsp2 <- read.csv("vsp2_ptpolbiassamp.csv")
vsp2 <- vsp2[, 1:2]

# Generate kernel rasters
vsp2.kn <- genkernel(vsp2[,1:2])

# Thin occurrences
vsp2.thin <- occthinner(vsp2, vsp2.kn, t1 = 0.5, t2 = 1, cols = 1:2)

vsp2.tsel <- vsp2.thin[[which(unlist(lapply(vsp2.thin, nrow)) ==
                                max(unlist(lapply(vsp2.thin, nrow))))[1]]]

# Spatial thinning
vsp2.tselsp <- thin(cbind(vsp2.tsel, sp = "vsp2"),
                       lat.col = "y", long.col = "x", spec.col = "sp",
                       thin.par = 10, reps = 10, locs.thinned.list.return = T,
                       write.files = F, write.log.file = F)

vsp2.tselsp <- vsp2.tselsp[[which(unlist(lapply(vsp2.tselsp, nrow)) ==
                                 max(unlist(lapply(vsp2.tselsp, nrow))))[1]]] 

# Save files
write.csv(vsp2.tsel, "vsp2_kernthin.csv", row.names = F)

write.csv(vsp2.tselsp, "vsp2_spatialthin.csv", row.names = F)

setwd("../..")
