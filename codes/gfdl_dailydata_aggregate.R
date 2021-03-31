### GFDL data aggregate
library(raster)

files <- list.files("data/oisst/daily_data", full.names = T)

daily.sst <- stack(files)

plot(daily.sst[[1:4]])

writeRaster(daily.sst, "data/oisst/daily_ssp85_gfdlnoaa.tif")
