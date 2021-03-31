#### Modelling of Herbivorous invertebrates of coral reefs ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2020

### Block generation - Sea-urchins ###

# Load packages ----
library(raster)
library(blockCV)
library(sf)


# Load functions and data----
#Load environmental layers
source("functions/Varload.r")
env <-
        var.load(folder = "crop_layers/", layers = "env_layers.txt",
                bath = "2_200")

#Set seed for replicability
set.seed(2932)

# Establish species codes and block size ----

#Establish the best block size
sac <- spatialAutoRange(
        rasterLayer = env,
        sampleNumber = 5000,
        doParallel = TRUE,
        showPlots = TRUE
)

chos.range <- sac$range

#Define species
sp.codes <- c("eclu", "lyva", "trve")

# Create function for block generation ----

block.gen <- function(species){
        
        #Load files
        pts <-
                read.csv(paste("data/", species, "/", species, "_pts.csv", sep =
                                       ""))
        
        dsp <-
                read.csv(paste("data/", species, "/", species, "_dsp.csv", sep =
                                       ""))
        
        #Prepare data.table
        sp.data <- pts
        sp.data$occ <- dsp$dsp
        sp.data$occ[is.na(sp.data$occ)] <-
                0 #Block_CV need binary numbers 0/1, and biomod NA for PA.
        
        # Make a SpatialPointsDataFrame object from data.frame
        pa_data <-
                st_as_sf(
                        sp.data,
                        coords = c("decimalLongitude", "decimalLatitude"),
                        crs = crs(env)
                )
        
        # spatial blocking by specified range with random assignment
        
        sb <- spatialBlock(
                speciesData = pa_data,
                species = "occ",
                rasterLayer = env,
                theRange = chos.range,
                k = 5,
                selection = "random",
                iteration = 100,
                biomod2Format = TRUE,
                showBlocks = FALSE
        )
        
        #Create a file with the split table for BIOMOD2
        split.table <- as.data.frame(sb$biomodTable)
        
        #Save split table
        write.csv(
                split.table,
                paste("data/", species, "/", species,
                      "_splittable.csv", sep = ""),
                row.names = F
        )
        
        #save index for BIOMOD Tuning
        write.table(sb$foldID,
                    paste("data/", species, "/", species,
                          "_foldindex.csv", sep = ""))
}

# Apply for each species ----

lapply(sp.codes, block.gen)

#END of code