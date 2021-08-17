#### Modelling of coral reef herbivorous invertebrates ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2021

### Threshold maps by max and min temperature experienced by sea urchins ###

# Load libraries and data ----
library(raster)
library(progress)
library(fs)

load("data/sst_limits/allspecies_oisst_thvalues.RData")

### All species in current period ----
for (j in 1:3) {
        
        # List species
        species <- c("lyva", "eclu", "trve")[j]
        
        # Load NOAA SST files
        current <- list.files("data/sst_limits/noaa", full.names = T)
        current <- lapply(current, brick)
        current <- stack(current)
        
        # Rotate to make -180 to 180
        current <- rotate(current)
        
        # Load bath layer (to mask NOAA file)
        m <- raster("data/env/bath_layers/bath_2_100_oisst.tif")
        
        # Create a progress bar to monitor proccess (approx 1 h per species)
        pb <- progress_bar$new(total = nlayers(current),
                               format = "[:bar] :percent :eta")
        
        # Load threshold values obtained by "sst_limits_threshold.R"
        maxlim <- thresholds[[species]]$maximum_limit
        minlim <- thresholds[[species]]$minimum_limit
        
        # Create a raster stack to receive the file
        thr <- stack()
        
        # Run for all layers
        for (i in 1:nlayers(current)) {
                
                # Load, crop and mask
                r <- current[[i]]
                
                nam <- names(r)
                
                r <- crop(r, m)
                
                r <- mask(r, m)
                
                # Threshold according to limits
                r <- calc(r, function(x){
                        x <- ifelse(x > maxlim,
                                    0, ifelse(x < minlim, 0, 1))
                        x
                })
                
                # In the first run stack, in the others sum
                if (i == 1) {
                        thr <- stack(thr, r)
                } else{
                        thr <- thr + r
                }
                
                pb$tick() # tick the progress bar
        }
        
        # Divide by the number of layers (convert it to decimal)
        thr <- thr/nlayers(current)
        
        # save file
        writeRaster(thr, paste0("data/sst_limits/",
                                species, "_current_thresh.tif"))
        
}

rm(list = ls()) # clean the environment

###  All species in future period ----

# Load threshold data
load("data/sst_limits/allspecies_oisst_thvalues.RData")

# Run in loop for all SSP scenarios
for (k in 1:4) {
        
        ssp <- c("ssp126", "ssp245", "ssp370", "ssp585")[k]
        
        # Run in loop for all species
        for (j in 1:3) {
                
                species <- c("lyva", "eclu", "trve")[j]
                
                # Load files
                clayers <- list.files(paste0("data/sst_limits/cmip6/", ssp,
                                             "/mean"),
                                      full.names = T)
                
                # Set a progress bar
                pb <- progress_bar$new(total = length(clayers),
                                       format = "[:bar] :percent :eta")
                
                # Get thresholds
                maxlim <- thresholds[[species]]$maximum_limit
                minlim <- thresholds[[species]]$minimum_limit
                
                # Creates a stack to hold data
                thr <- stack()
                
                for (i in 1:length(clayers)) {
                        
                        # Load layer
                        r <- raster(clayers[i])
                        
                        nam <- names(r)
                        
                        # Apply threshold
                        r <- calc(r, function(x){
                                x <- ifelse(x > maxlim,
                                            0, ifelse(x < minlim, 0, 1))
                                x
                        })
                        
                        # stack/sum layer
                        if (i == 1) {
                                thr <- stack(thr, r)
                        } else{
                                thr <- thr + r
                        }
                        
                        pb$tick() # tick the progress bar
                }
                
                # Divide by the number o layers
                thr <- thr/length(clayers)
                
                # Save file
                writeRaster(thr, paste0("data/sst_limits/",
                                        species, "_", ssp, "_thresh.tif"))
                
        }
}

### END of code