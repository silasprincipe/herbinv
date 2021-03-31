eclu <- raster("data/oisst/frequency_eclu_thresh.tif")
lyva <- raster("data/oisst/frequency_lyva_thresh.tif")
trve <- raster("data/oisst/frequency_trve_thresh.tif")

val <- values(eclu)
val[is.nan(val)] <- NA
val <- max(val, na.rm = T)

eclu <- eclu/val
lyva <- lyva/val
trve <- trve/val

calcx <- function(x){
        x[x<0.8] <- 0
        x[x >=0.8] <- 1
        x
}

eclu <- calc(eclu, calcx)
lyva <- calc(lyva, calcx)
trve <- calc(trve, calcx)

epts <- read.csv("data/eclu/eclu_cell_noaa.csv")[,1:2]
lpts <- read.csv("data/lyva/lyva_cell_noaa.csv")[,1:2]
tpts <- read.csv("data/trve/trve_cell_noaa.csv")[,1:2]


plot(eclu)
points(epts, col = 'black', pch = 20, cex = 0.8)
plot(lyva)
points(lpts, col = 'black', pch = 20, cex = 0.8)
plot(trve)
points(tpts, col = 'black', pch = 20, cex = 0.8)
