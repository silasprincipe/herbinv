# Get TSM (based on maximum temp) for sea-urchins

# Load libraries ----
library(raster)
library(tidyverse)

# Establish species
species <- c("eclu")
round.code <- "ih_c1"

# Establish Tmax for the species
tmax <- c(36)

# Establish threshold function
thresh <- function(x){
        x[x < 667] <- 0
        x[x >= 667] <- 1
        x
}

# Import SDM result (current)
cur <- raster(paste0(species, "/proj_current_", species, "_", round.code,
                     "/individual_projections/",
                     species, "_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img"))

# Apply threshold
cur <- calc(cur, thresh)

cur <- na.omit(cbind("values" = values(cur), data.frame(xyFromCell(cur, 1:ncell(cur)))))
cur <- cur[cur$values == 1,]

c.samp <- sample_n(cur, size = 1000)

# Load SST data (bio-oracle)
sst <- raster("data/env/crop_layers/BO21_tempmax_ss.tif")

# Extract temp data
s.sst <- raster::extract(sst, c.samp[,2:3])
s.sst <- s.sst - tmax
s.sst <- cbind(s.sst, c.samp[,2:3])

# Plot to see
library(rnaturalearth)
library(rnaturalearthdata)

world <- ne_countries(scale = "medium", returnclass = "sf")

ggplot(data = world) +
        geom_sf()+
        coord_sf(xlim = c(-99, -29), ylim = c(-42.5, 42.5))+
        geom_point(data = s.sst, aes(x = x, y = y, color = s.sst))+
        scale_color_continuous(type = "viridis")

# Future sst temperatures
sst.126 <- raster("data/env/proj_layers/ssp126/BO21_tempmax_ss.tif")
sst.245 <- raster("data/env/proj_layers/ssp245/BO21_tempmax_ss.tif")
sst.585 <- raster("data/env/proj_layers/ssp585/BO21_tempmax_ss.tif")

# Extract temp data
s.sst.126 <- raster::extract(sst.126, c.samp[,2:3])
s.sst.126 <- s.sst.126 - tmax
s.sst.126 <- cbind(s.sst.126, c.samp[,2:3])

# Extract temp data
s.sst.245 <- raster::extract(sst.245, c.samp[,2:3])
s.sst.245 <- s.sst.245 - tmax
s.sst.245 <- cbind(s.sst.245, c.samp[,2:3])

# Extract temp data
s.sst.585 <- raster::extract(sst.585, c.samp[,2:3])
s.sst.585 <- s.sst.585 - tmax
s.sst.585 <- cbind(s.sst.585, c.samp[,2:3])

ggplot(data = world) +
        geom_sf()+
        coord_sf(xlim = c(-99, -29), ylim = c(-42.5, 42.5))+
        geom_point(data = s.sst.126, aes(x = x, y = y, color = s.sst.126))+
        scale_color_continuous(type = "viridis")

ggplot(data = world) +
        geom_sf()+
        coord_sf(xlim = c(-99, -29), ylim = c(-42.5, 42.5))+
        geom_point(data = s.sst.245, aes(x = x, y = y, color = s.sst.245))+
        scale_color_continuous(type = "viridis")

ggplot(data = world) +
        geom_sf()+
        coord_sf(xlim = c(-99, -29), ylim = c(-42.5, 42.5))+
        geom_point(data = s.sst.585, aes(x = x, y = y, color = s.sst.585))+
        scale_color_continuous(type = "viridis")
