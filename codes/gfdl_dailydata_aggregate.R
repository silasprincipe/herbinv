#### Modelling of coral reef herbivorous invertebrates ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2021

# Aggregate daily data on GFDL-NOAA climate projections of SST #
# Data was first interpolated using another code and in parallel
# To supply this data to GEE we need to aggregate it in a tif file,
# being each day a layer/band

# Load library ----
library(raster)

# List files ----
files <- list.files("data/oisst/daily_data", full.names = T, pattern = ".tif")

# Load bath layer ----
bath <- raster("data/env/bath_layers/bath_2_100_oisst.tif")

for (i in 1:11) {
        
        base <- seq(from = 365, to = 4015, by = 365)
        years <- seq(2090, 2100)
        
        # Stack files, mask and save ----
        if (i == 1) {
                daily.sst <- stack(files[1:base[i]])     
        } else{
                daily.sst <- stack(files[(base[i - 1]+1):base[i]])
        }
        
        daily.sst <- mask(daily.sst, bath)
        
        #plot(daily.sst[[1:4]]) # We can plot some to verify
        
        writeRaster(daily.sst,
                    paste0("data/oisst/daily_ssp85_gfdlnoaa_",
                           years[i], ".tif"),
                    overwrite = T) 
        
        print(names(daily.sst))
        
}

##END