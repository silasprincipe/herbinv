#### Modelling of coral reef herbivorous invertebrates ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2021

### Sea-urchins data conversion to 1 per cell - NOAA version ###
# This data will be used in GEE for other part of the work
# It uses as base the OISST data from NOAA in 0.25° resolution

# Load packages ----
library(tidyverse)
library(raster)

# Convert species points to one per cell ----

# List species
sp.codes <- c("eclu", "lyva", "trve")

# Load environmental layers
env <- raster(list.files("data/temperature/noaa_oisst/", full.names = T))

env <- rotate(env)

plot(env)

bath <- raster("data/env/bath_layers/bath_2_100_oisst.tif")

env <- crop(env, bath)
env <- mask(env, bath)

# Create a function to convert to 1 point per cell
to.cell <- function(species){
        
        sp <- read.csv(paste0("data/", species, "/", species, "_final.csv"))
        
        #Remove points out of the area/bathymetry
        sp2 <- sp
        coordinates(sp2) <- ~ decimalLongitude + decimalLatitude
        crs(sp2) <- crs(env)
        out.p <- raster::extract(env, sp2)
        sp <- sp[!is.na(out.p),]
        rm(sp2, out.p)
        
        #Create a clean raster
        r <- env
        values(r) <- NA
        
        #Put 1 in presence cells
        r[cellFromXY(r, sp[,2:3])] <- 1
        
        #Convert raster to dataframe
        data <- as.data.frame(coordinates(r))
        data$sp <- as.vector(values(r))
        
        #Change column names
        colnames(data) <- c("decimalLongitude", "decimalLatitude",
                            as.character(species))
        
        #Remove NAs
        sp <- data %>% filter(!is.na(data[,3]))
        
        sp
}

# Apply to datasets
cell.data <- lapply(sp.codes, to.cell)

# Save datasets ----
for (i in 1:3) {
        
        write.csv(cell.data[[i]], paste0("data/",
                                          sp.codes[i],
                                          "/", sp.codes[i],
                                          "_cell_noaa.csv"),
                  row.names = F)
}

#END of code