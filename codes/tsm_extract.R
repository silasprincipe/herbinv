library(raster)
library(tidyverse)

#Load ensemble

eclu <- raster("eclu/proj_current_eclu_iht5/individual_projections/eclu_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")

eclu[eclu < 667] <- 0
eclu[eclu >= 667] <- 1

plot(eclu)

# Extract data

data <- data.frame(cbind(xyFromCell(eclu, 1:ncell(eclu)), "values" = values(eclu)))

data$values[is.na(data$values)] <- 0

data <- data[data$values == 1,]

data.s <- sample_n(data, size = 1000)

points(data.s[,1:2], col = 'blue', pch = 20, cex = 0.5)

# Extract temperature

sst <- raster("data/env/crop_layers/BO21_tempmax_ss.tif")

sst.data <- raster::extract(sst, data.s[,1:2])

hist(sst.data)

summary(sst.data)

eclu.opt <- 29.4

abline(v = eclu.opt, col = 'red')


tsm <- sst.data - eclu.opt

hist(tsm)

tsm.thresh <- cbind(data.s, tsm)

tsm.thresh$values[tsm.thresh$tsm > 0] <- "D"
tsm.thresh$values[tsm.thresh$tsm <= 0] <- "O"




plot(eclu)
points(tsm.thresh[tsm.thresh$values == "D",1:2], col = 'red', pch = 20, cex = 0.5)
points(tsm.thresh[tsm.thresh$values == "O",1:2], col = 'blue', pch = 20, cex = 0.5)




# If warming tolerance

# Extract temperature

sst <- raster("data/env/crop_layers/BO21_tempmax_ss.tif")

sst.data <- raster::extract(sst, data.s[,1:2])

hist(sst.data)

summary(sst.data)

eclu.opt <- 36

abline(v = eclu.opt, col = 'red')


tsm <- sst.data - eclu.opt

hist(tsm)

tsm.thresh <- cbind(data.s, tsm)

tsm.thresh$values[tsm.thresh$tsm > 0] <- "D"
tsm.thresh$values[tsm.thresh$tsm <= 0] <- "O"




plot(eclu)
points(tsm.thresh[tsm.thresh$values == "D",1:2], col = 'red', pch = 20, cex = 0.5)
points(tsm.thresh[tsm.thresh$values == "O",1:2], col = 'blue', pch = 20, cex = 0.5)

library(rgdal)
library(rasterVis)
library(RColorBrewer)



ggplot(tsm.thresh)+geom_point(aes(x = x, y = y, color = tsm))+
        scale_colour_gradientn(colours = terrain.colors(10))
