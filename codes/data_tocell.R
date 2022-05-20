#### Modelling of coral reef herbivorous invertebrates ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2021

### Sea-urchins data conversion to 1 per cell ###
# This file includes both the regular version (0.083 resolution)
# and the NOAA version (0.25 resolution). See below.

# Load packages ----
library(raster)

# Regular version ----

# List species
sp.codes <- c("eclu", "lyva", "trve")

# Load environmental layers (already prepared using bath layer as mask)
base <- raster("data/env/crop_layers/BO21_tempmean_ss.tif")

# Create a function to convert to 1 point per cell
to.cell <- function(species, plot = T){
        
        sp <- read.csv(paste0("data/", species, "/", species, "_final.csv"))
        
        ## For Tripneustes ventricosus we add a point in Fernando de Noronha
        ## According to field collection done by the authors
        if (species == "trve") {
                sp <- rbind(sp, data.frame(
                        scientificName = "trve",
                        decimalLongitude = -32.43,
                        decimalLatitude = -3.8
                ))
        }
        
        if (plot) {
                par(mfrow = c(1,2))
                plot(base, col = "grey", legend = F,
                     main = paste(species, "original"))
                points(sp[,2:3], pch = 16, cex = 0.7, col = "blue")
        }
        
        n <- nrow(sp)
        
        #Remove points out of the area/bathymetry
        out.p <- raster::extract(base, sp[,2:3])
        sp <- sp[!is.na(out.p),]
        
        #Create a clean raster
        r <- base
        r[] <- NA
        
        #Put 1 in presence cells
        r[cellFromXY(r, sp[,2:3])] <- 1
        
        #Convert raster to dataframe
        data <- rasterToPoints(r)
        
        #Change column names
        colnames(data) <- c("decimalLongitude", "decimalLatitude",
                            species)
        
        if (plot) {
                plot(base, col = "grey", legend = F,
                     main = "thinned");points(data[,1:2], pch = 16, cex = 0.7,
                                               col = "blue")
        }
        
        cat(species, "had", n, "points.", nrow(data), "remained. \n")
        
        data
}

# Apply to datasets
cell.data <- lapply(sp.codes, to.cell)

# Save datasets
for (i in 1:3) {
        
        write.csv(cell.data[[i]], paste0("data/",
                                          sp.codes[i],
                                          "/", sp.codes[i],
                                          "_cell.csv"),
                  row.names = F)
}

rm(base, cell.data)


# NOAA version (resolution of 0.25) ----
# This data will be used in other part of the work
# It uses as base the OISST data from NOAA in 0.25° resolution

# Load environmental layers
env <- raster(list.files("data/sst_limits/noaa/", full.names = T)[1])
env <- rotate(env)

plot(env)
bath <- raster("data/env/bath_layers/bath_2_100_oisst.tif")

env <- crop(env, bath)
base <- mask(env, bath)
plot(base)

# Apply to datasets
cell.data <- lapply(sp.codes, to.cell)

# Save datasets
for (i in 1:3) {
        
        write.csv(cell.data[[i]], paste0("data/",
                                         sp.codes[i],
                                         "/", sp.codes[i],
                                         "_cell_noaa.csv"),
                  row.names = F)
}

#END of code