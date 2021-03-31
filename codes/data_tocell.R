#### Modelling of Herbivorous invertebrates of coral reefs ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2020

### Sea-urchins data conversion to 1 per cell ###

# Load packages ----
library(tidyverse)
library(raster)

# Load functions ----
source("functions/varload.r")

# Convert species points to one per cell ----

# List species
sp.codes <- c("eclu", "lyva", "trve")

# Load environmental layers
env <- var.load(folder = "crop_layers/", 
                layers = "env_layers.txt")

env <- mean(env) # We do this as we just need the cells for which we have info

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
        r <- raster(nrow=1020, ncol=1452,
                    xmn = -99, xmx= 22,
                    ymn= -42.5, ymx = 42.5)
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
                                          "_cell.csv"),
                  row.names = F)
}

#END of code