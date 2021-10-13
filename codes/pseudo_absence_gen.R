#### Modelling of Herbivorous invertebrates of coral reefs ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2021

### Pseudo-absences and block generation ###

# Here we first generate the pseudo-absences using a pre-model (Mahalanobis)
# Then we generate blocks for the cross-validation. We cross-validate using
# two methods of blocking: random blocks and latitudinal blocks.
# For more info see the methods section.

# Pseudo-absences generation =====

# Load libraries/functions ----
library(raster)
library(dismo)
library(adehabitatHS)
library(rgeos)
source("functions/Varload.r")

# Set seed for replicability
set.seed(2932)

# Load environmental layers ----
env <- var.load(folder = "crop_layers/", "env_layers.txt", "2_100")

#Ensure all layers have the same NA values
env <- mask(env, env[[1]])

# Prepare data to mahalanobis
# Convert environmental layers to SpatialGrid
env.data <- as(env, "SpatialPixelsDataFrame")

slot(env.data, "data") <-
        dudi.pca(slot(env.data, "data"), scannf = FALSE)$tab

# Set pa generation parameters ----
pan <- 10 # The multiplier for the number of pseudo-absences generated
parn <- 10 # The number of pseudo-absences datasets to be generated

# Run for all species in loop ----
species <- c("lyva", "eclu", "trve")

for (i in species) {
        
        # Load species data
        sp.data <- read.csv(paste0("data/", i, "/", i, "_knthin.csv"))
        
        sp.data <- sp.data[,1:2]
        
        # Ensure consistency
        names(sp.data) <- c("decimalLongitude", "decimalLatitude")
        
        coordinates(sp.data) <- ~ decimalLongitude + decimalLatitude
        crs(sp.data) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
        
        #Correct CRS
        crs(env.data) <- crs(sp.data)
        
        #Habitat suitability mapping with Mahalonobis distance
        hsm <- mahasuhab(env.data, sp.data, type = "probability")
        # crs(hsm) <- crs(sp.data)
        # mtp <- sp::over(sp.data, hsm)
        # mtp <- min(na.omit(mtp))
        # p10 <- ceiling(length(mtp$MD) * 0.9)
        # thresh <- rev(sort(mtp$MD))[p10]
        
        # 3: Remove cells from the suitable area of Mahalonobis distance
        c.mh <- data.frame(cbind(hsm@coords, hsm@data))
        c.mh <- c.mh[c.mh$MD > 0.01, ] # We ignore points with very low values 
        # of suitability (less than 1%)
        
        # Create a mask for the sampling
        mask.env <- env[[1]]
        mask.env[cellFromXY(mask.env, c.mh[, 1:2])] <- NA
        
        # Creates a buffer around presences
        x <- raster::buffer(sp.data, width = 40000)
        
        # Set all raster cells inside the buffer to NA
        x <- mask(mask.env, x, inverse = T)
        
        # Ensure that adjacent cells are not sampled
        ma.pts <- gridSample(rasterToPoints(x)[,1:2], x, chess = "white")
        
        # Generate parn number of datasets
        for (z in 1:parn) {
                
                # Sample pa points
                ma.pts.samp <- data.frame(ma.pts[sample(1:nrow(ma.pts),
                                                       (length(sp.data)*pan)),])
                
                # Get original points
                sp.pts <- data.frame(sp.data@coords)
                sp.pts$dsp <- 1
                
                # Set pa as 0
                ma.pts.samp$dsp <- 0
                
                names(ma.pts.samp) <- names(sp.pts)
                
                # Write file
                write.csv(rbind(sp.pts, ma.pts.samp),
                          paste0("data/", i, "/", i, "_pa_", z, ".csv"),
                          row.names = F)
                
        }
        
        cat(i, "pseudo-absence generated. \n")
}


# Block-cross validation folds generation =====
# Load libraries ----
library(blockCV)
library(sf)
library(sp)

# Calculate best block size ----
sac <- spatialAutoRange(
        rasterLayer = env,
        sampleNumber = 5000,
        doParallel = F,
        showPlots = F
)

chos.range <- sac$range

# Run in loop for all species ----
for (i in species) {
        files <- list.files(paste0("data/", i))
        files <- files[grep("pa", files)]
        files <- files[-grep("folds", files)]
        
        for (z in 1:length(files)) {
                # Read file
                sp <- read.csv(paste0("data/", i, "/", files[z]))
                
                # Convert to sf
                pa_data <- st_as_sf(sp,
                                    coords = c("decimalLongitude",
                                               "decimalLatitude"),
                                    crs = crs(env))
                
                ##### Get random blocks folds -----
                sb <- spatialBlock(
                        speciesData = pa_data,
                        species = "dsp",
                        rasterLayer = env,
                        theRange = chos.range,
                        k = 5,
                        selection = "random",
                        iteration = 100,
                        biomod2Format = FALSE,
                        showBlocks = FALSE
                )
                
                #save fold index
                write.table(sb$foldID,
                            paste0("data/", i, "/",
                                   gsub(".csv", "", files[z]), "_folds.csv"),
                            row.names = F)
                
                #### Get latitudinal folds ----
                # r <- raster(nrow = 15, ncol = 1)
                # r[] <- c(rep(1:5, 3))
                # extent(r) <- extent(-99,-29,-42.5,42.5)
                # pols <- rasterToPolygons(r)
                # 
                # coordinates(sp) <- ~decimalLongitude + decimalLatitude
                # crs(pols) <- crs(sp)
                # 
                # latf <- over(sp, pols)
                # 
                # cat("===== Latitudinal blocks: \n")
                # cat("TR0", "TR1", "TE0", "TE1", "\n")
                # for (j in 1:5) {
                #         cat(length(sp$dsp[sp$dsp == 0 & latf[,1] != j]),
                #             length(sp$dsp[sp$dsp == 1 & latf[,1] != j]),
                #             length(sp$dsp[sp$dsp == 0 & latf[,1] == j]),
                #             length(sp$dsp[sp$dsp == 1 & latf[,1] == j]), "\n")
                # }
                # 
                # #save fold index
                # write.table(latf[,1],
                #             paste0("data/", i, "/",
                #                    gsub(".csv", "", files[z]), "_latfolds.csv"),
                #             row.names = F)
                
                sbl <- spatialBlock(
                        speciesData = pa_data,
                        species = "dsp",
                        rasterLayer = env,
                        rows = 15,
                        k = 5,
                        selection = "random",
                        iteration = 100,
                        biomod2Format = FALSE,
                        showBlocks = FALSE
                )
                
                #save fold index
                write.table(sbl$foldID,
                            paste0("data/", i, "/",
                                   gsub(".csv", "", files[z]), "_latfolds.csv"),
                            row.names = F)
                
        }
        
        cat(i, "block generation done. \n")
}

# END