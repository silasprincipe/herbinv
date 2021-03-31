#### Modelling of Herbivorous invertebrates of coral reefs ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2020

### Pseudo-absence generation - Sea-urchins ###

# Load packages ----
library(tidyverse)
library(raster)
library(adehabitatHS)
library(dismo)
library(sf)
library(rnaturalearth)


# Load functions and data----
source("functions/varload.r")
world <- ne_countries(scale = "medium", returnclass = "sf")

# Define species codes and load env data ----

#Set seed for replicability
set.seed(2932)

#Define species acronyms
sp.codes <- c("eclu", "lyva", "trve")

#Load environmental data
env <-
        var.load(folder = "crop_layers/",
                layers = "env_layers.txt",
                bath = "2_100")

#Correct bath layer to have the same NA values than others
bath <- mask(env$bath_2_100, env$BO_calcite)
env$bath_2_100 <- bath

#Convert environmental layers to SpatialGrid
sp.data <- as(env, "SpatialPixelsDataFrame")

#Prepare data to Mahalonobis
slot(sp.data, "data") <-
        dudi.pca(slot(sp.data, "data"), scannf = FALSE)$tab


# Create function to generate pseudo-absence ----

pseudo.abs <- function(species){
        
        #Define path to PA preparing
        path <- paste("data/", species, "/", sep = "")
        
        #Read species data
        sp <- read.csv(paste(path, species, "_cell.csv", sep = ""))
        
        #Convert occurrence points to SpatialPoints
        occs <-
                SpatialPointsDataFrame(sp[sp[, 3] == 1, 1:2], 
                                       data.frame(sp[sp[, 3] == 1, ]))
        
        #Correct CRS
        crs(occs) <- crs(sp.data)
        
        #Habitat suitability mapping with Mahalonobis distance
        hsm <- mahasuhab(sp.data, occs, type = "probability")
        
        # Save Mahalonobis habitat suitability maps (can be skipped if wanted)
        hsm_sf <- st_as_sf(hsm)
        st_write(hsm_sf,
                 paste0("gis/hsm_", species, ".shp"),
                 driver = "ESRI Shapefile",
                 append = FALSE)
        
        # Get only lon lat info
        sp.p <- sp[sp[, 3] == 1, 1:2]
        
        #Create a buffer around presence points (20km)
        buf1 <- circles(sp.p, lonlat = T, d = 20000)
        
        
        ###Get the P10 threshold - code based on the one from Cecina B. Morrow
        # https://babichmorrowc.github.io/post/2019-04-12-sdm-threshold/
        
        #Separate species data
        sp.t <- sp
        sp.t <- sp.t[sp.t[,3] == 1,1:2]
        
        #Create an empty raster
        tr <- raster(nrow = 1020, ncol = 840, xmn= -99, xmx = -29, ymn= -42.5, ymx=42.5)
        values(tr) <- NA
        
        #Include values from Mahalonobis model on the raster
        tr[cellFromXY(tr, hsm@coords)] <- hsm@data$MD
        
        #Extract Mahalonobis values from the occ. points
        occPredVals <- raster::extract(tr, sp.t)
        
        #Get the 10th percentile value
        if(length(occPredVals) < 10){
                p10 <- floor(length(occPredVals) * 0.9)
        } else {
                p10 <- ceiling(length(occPredVals) * 0.9)
        }
        
        #Get the threshold
        thresh <- rev(sort(occPredVals))[p10]
        
        ###END of P10 threshold
        
        
        #Remove cells from the suitable area of Mahalonobis distance
        c.mh <- data.frame(cbind(hsm@coords, hsm@data))
        c.mh <- c.mh[c.mh$MD >= thresh, ] 
        env[cellFromXY(env, c.mh[, 1:2])] <- NA
        
        #Sample random points in the environmental area
        samp <-
                as.data.frame(sampleRandom(
                        env,
                        size = 10000,
                        na.rm = TRUE,
                        xy = TRUE
                )[, 1:2])
        
        #Convert to spatial points
        samp.sp <- SpatialPoints(coords = samp[, 1:2],
                                 proj4string = CRS('+proj=longlat +datum=WGS84'))
        
        #Remove points that are in the buffered area (close to presences)
        out <- over(samp.sp, geometry(buf1))
        
        samp.sel <- samp[is.na(out), ]
        
        #Sample final points
        samp.final <-
                sample_n(samp.sel, length(sp.p$decimalLongitude) * 2)
        
        #Create the dataset in the final format
        samp.final$dsp <- NA
        samp.final$PA1 <- TRUE
        
        colnames(samp.final) <- c("decimalLongitude",
                                  "decimalLatitude",
                                  "dsp",
                                  "PA1")
        
        
        #Create datasets from biblio/obis/gbif (true dataset)
        pres <- samp.final[0, ]
        
        pres <- bind_rows(pres, sp.p)
        pres$dsp <- 1
        pres$PA1 <- TRUE
        
        #Bind all datasets
        pa.table <- bind_rows(pres, samp.final)
        
        #Write final file with full table
        write.csv(pa.table,
                  paste(path, species, "_pa_table.csv", sep = ""),
                  row.names = F)
        
        #Write separate final files (for use with biomod2)
        pts <- pa.table[, 1:2]
        dsp <- as.data.frame(pa.table[, 3])
        colnames(dsp) <- "dsp"
        set <- data.frame(pa.table[, 4])
        colnames(set) <- "PA1"
        
        write.csv(pts,
                  paste(path, species, "_pts.csv", sep = ""),
                  row.names = F)
        write.csv(dsp,
                  paste(path, species, "_dsp.csv", sep = ""),
                  row.names = F)
        write.csv(set,
                  paste(path, species, "_set.csv", sep = ""),
                  row.names = F)
        
        ret <- list()
        
        ret[[1]] <- table(pa.table$dsp[pa.table$PA1 == T], useNA = "ifany")
        
        ret[[2]] <- ggplot(data = world) +
                geom_sf() +
                xlim(-99,-29) +
                ylim(-42.5, 42.5)+
                geom_point(data = pa.table, aes(x = decimalLongitude,
                                                y = decimalLatitude,
                                                color = dsp))+
                theme_classic()
        
        names(ret) <- c("table", "plot")
        
        ret

}

# Apply pseudo-absence generation to each dataset ----

t.results <- lapply(sp.codes, pseudo.abs)

#Verify
t.results[[3]]$table

t.results[[3]]$plot #just change the number to verify each

#END of code