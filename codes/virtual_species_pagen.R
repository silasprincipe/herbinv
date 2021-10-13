#### Modelling of Herbivorous invertebrates of coral reefs ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2021

### Testing pseudo-absence generation methods ###
# In this file we will generate the pseudo-absence datasets with three methods
# 1. Random with small buffer
# 2. Random with larger buffer
# 3. Mahalanobis distance (pre-model)

# Load packages and functions ----
library(raster)
library(dismo)
library(adehabitatHS)
library(rgeos)
source("functions/Varload.r")

# Set seed for replicability
set.seed(2932)

# Load environmental layers
env <- var.load(folder = "crop_layers/", "env_layers.txt", "2_100")

#Ensure all layers have the same NA values
env <- mask(env, env[[1]])

# Create a mask for the sampling
mask.env <- env[[1]]

# Load presence files ----
setwd("data/vsp/")

# Creates a list with all data
files <- list.files(pattern = ".csv")
dlist <- lapply(files, read.csv)

# Get names for files
files <- gsub(".csv", "", files)

names(dlist) <- files

# Create a file to hold the pseudo-absence data
if (!dir.exists("padata")) {dir.create("padata")}
setwd("padata")

#### 1. Random with small buffer ----
getrandom.sb <- function(name, data.list, mask.env){
        
        sp.data <- data.list[[name]]
        
        sp.data <- sp.data[,1:2]
        
        names(sp.data) <- c("x", "y")
        
        coordinates(sp.data) <- ~ x + y
        crs(sp.data) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
        
        x <- raster::buffer(sp.data, width=20000)
        
        # Set all raster cells inside the buffer to NA
        x <- mask(mask.env, x, inverse = T)
        
        sb.pts <- gridSample(rasterToPoints(x)[,1:2], x, chess = "white")
        
        
        for (i in 1:5) {
                
                sb.pts.2t <- data.frame(sb.pts[sample(1:nrow(sb.pts),
                                                      (length(sp.data)*2)),])
                sb.pts.10t <- data.frame(sb.pts[sample(1:nrow(sb.pts),
                                                       (length(sp.data)*10)),])
                
                sp.pts <- data.frame(sp.data@coords)
                sp.pts$dsp <- 1
                
                sb.pts.2t$dsp <- 0
                sb.pts.10t$dsp <- 0
                
                write.csv(rbind(sp.pts, sb.pts.2t),
                          paste0(name, "_sbpa_2t_", i, ".csv"), row.names = F)
                
                write.csv(rbind(sp.pts, sb.pts.10t),
                          paste0(name, "_sbpa_10t_", i, ".csv"), row.names = F)
        }
        
        name
}

lapply(files, getrandom.sb, dlist, mask.env)


#### 2. Random with larger buffer ----
getrandom.lb <- function(name, data.list, mask.env){
        
        sp.data <- data.list[[name]]
        
        sp.data <- sp.data[,1:2]
        
        names(sp.data) <- c("x", "y")
        
        coordinates(sp.data) <- ~ x + y
        crs(sp.data) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
        
        x <- raster::buffer(sp.data, width=200000)
        
        # Set all raster cells inside the buffer to NA
        x <- mask(mask.env, x, inverse = T)
        
        lb.pts <- gridSample(rasterToPoints(x)[,1:2], x, chess = "white")
        
        for (i in 1:5) {
                
                lb.pts.2t <- data.frame(lb.pts[sample(1:nrow(lb.pts),
                                                      (length(sp.data)*2)),])
                lb.pts.10t <- data.frame(lb.pts[sample(1:nrow(lb.pts),
                                                       (length(sp.data)*10)),])
                
                sp.pts <- data.frame(sp.data@coords)
                sp.pts$dsp <- 1
                
                lb.pts.2t$dsp <- 0
                lb.pts.10t$dsp <- 0
                
                write.csv(rbind(sp.pts, lb.pts.2t),
                          paste0(name, "_lbpa_2t_", i, ".csv"), row.names = F)
                
                write.csv(rbind(sp.pts, lb.pts.10t),
                          paste0(name, "_lbpa_10t_", i, ".csv"), row.names = F)
        }
        
        name
}

lapply(files, getrandom.lb, dlist, mask.env)


#### 3. Mahalonobis distance ----

# Convert environmental layers to SpatialGrid
env.data <- as(env, "SpatialPixelsDataFrame")

#Prepare data to Mahalanobis
slot(env.data, "data") <-
        dudi.pca(slot(env.data, "data"), scannf = FALSE)$tab

getmahala <- function(name, data.list, mahala.data, maskr){
        
        sp.data <- data.list[[name]]
        
        sp.data <- sp.data[,1:2]
        
        names(sp.data) <- c("x", "y")
        
        coordinates(sp.data) <- ~ x + y
        crs(sp.data) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
        
        #Correct CRS
        crs(mahala.data) <- crs(sp.data)
        
        #Habitat suitability mapping with Mahalonobis distance
        hsm <- mahasuhab(mahala.data, sp.data, type = "probability")
        
        # 3: Remove cells from the suitable area of Mahalonobis distance
        c.mh <- data.frame(cbind(hsm@coords, hsm@data))
        c.mh <- c.mh[c.mh$MD > 0.01, ] # We ignore points with very low values 
        # of suitability (less than 1%)
        maskr[cellFromXY(maskr, c.mh[, 1:2])] <- NA
        
        x <- raster::buffer(sp.data, width = 40000)
        
        # Set all raster cells inside the buffer to NA
        x <- mask(maskr, x, inverse = T)
        
        ma.pts <- gridSample(rasterToPoints(x)[,1:2], x, chess = "white")
        
        for (i in 1:5) {
                
                ma.pts.2t <- data.frame(ma.pts[sample(1:nrow(ma.pts),
                                                      (length(sp.data)*2)),])
                ma.pts.10t <- data.frame(ma.pts[sample(1:nrow(ma.pts),
                                                       (length(sp.data)*10)),])
                
                sp.pts <- data.frame(sp.data@coords)
                sp.pts$dsp <- 1
                
                ma.pts.2t$dsp <- 0
                ma.pts.10t$dsp <- 0
                
                write.csv(rbind(sp.pts, ma.pts.2t),
                          paste0(name, "_mapa_2t_", i, ".csv"), row.names = F)
                
                write.csv(rbind(sp.pts, ma.pts.10t),
                          paste0(name, "_mapa_10t_", i, ".csv"), row.names = F)
        }
        
        name
}

lapply(files, getmahala, dlist, env.data, maskr = mask.env)

### Generate blocks for cross validation
library(blockCV)
library(sf)

files <- list.files()

# Calculate best block size
sac <- spatialAutoRange(
        rasterLayer = env,
        sampleNumber = 5000,
        doParallel = F,
        showPlots = F
)

chos.range <- sac$range

genblocks <- function(file){
        
        sp <- read.csv(file)
        
        pa_data <- st_as_sf(sp,
                            coords = c("x", "y"),
                            crs = crs(env))
        
        sb <- spatialBlock(
                speciesData = pa_data,
                species = "dsp",
                rasterLayer = env,
                theRange = chos.range,
                k = 5,
                selection = "random",
                iteration = 5,
                biomod2Format = TRUE,
                showBlocks = FALSE
        )
        
        #save index for BIOMOD Tuning (if needed - not used on our case)
        write.table(sb$foldID,
                    paste0(gsub(".csv", "", file), "_folds.csv"),
                    row.names = F)
}

lapply(files, genblocks)

setwd("../../..")
