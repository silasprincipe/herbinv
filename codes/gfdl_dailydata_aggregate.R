#### Modelling of coral reef herbivorous invertebrates ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2021

# Aggregate daily data on GFDL-NOAA climate projections of SST #
# Data was first interpolated using another code and in parallel
# To supply this data to GEE we need to aggregate it in a tif file,
# being each day a layer/band

# Load library ----
library(raster)

# List files ----
files <- list.files("data/oisst/daily_data", full.names = T)

# Stack files and save ----
daily.sst <- stack(files)

plot(daily.sst[[1:4]]) # We can plot some to verify

writeRaster(daily.sst, "data/oisst/daily_ssp85_gfdlnoaa.tif")

##END