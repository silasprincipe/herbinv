#### Modelling of coral reef herbivorous invertebrates ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2022 ##

## Barrier Model under the Log-Gaussian Cox Process framework ##

## Prepare data to use with models ##
# - Reproject environmental data and save
# - Prepare INLA mesh and save
# - Get integration points and save
# - Prepare full covariate dataset and save


# Load packages and define settings ----
# Modeling
library(INLA)
library(inlabru)
# Spatial data
library(raster)
library(rgeos)
library(gstat)
# Plotting and utilities
library(ggplot2)
# Settings
set.seed(2932) # Replicability of sampling
spatt <- theme_classic() # Theme for better ploting


# Load environmental and species data ----
source("functions/Varload.r")
lays <- c("salinitymean", "tempmean", "ph", "chlomean", "phosphatemax")
biocl <- c("warmm", "coldm")

env <- var.load(layers = lays, bioclim = biocl)

r12 <- var.load(folder = "proj_layers/ssp126/", layers = lays,
                bioclim = biocl, bioclimfut = "ssp126")
r24 <- var.load(folder = "proj_layers/ssp245/", layers = lays,
                bioclim = biocl, bioclimfut = "ssp245")
r37 <- var.load(folder = "proj_layers/ssp370/", layers = lays,
                bioclim = biocl, bioclimfut = "ssp370")

# Load study area shapefile
starea <- shapefile("gis/starea.shp")

# Mask environmental layers
env <- mask(env, starea)
r12 <- mask(r12, starea)
r24 <- mask(r24, starea)
r37 <- mask(r37, starea)

# Change names for easier handling
names(env) <- c("ph", "chl", "pho", "sal", "sst", "coldm", "warmm")
names(r37) <- names(r24) <- names(r12) <- names(env)

# Put chlorophyll in log scale
env$chl <- log(env$chl)
r12$chl <- log(r12$chl)
r24$chl <- log(r24$chl)
r37$chl <- log(r37$chl)

# Scale variables
# Future variables are scaled according to the current period
# Get mean and SD values for current period:
m <- cellStats(env, 'mean')
sd <- cellStats(env, 'sd')

# Use auto scaling for current
env <- scale(env)

# Manually scale future based on current period
r12 <- (r12 - m)/sd
r24 <- (r24 - m)/sd
r37 <- (r37 - m)/sd

# Reproject everything to the Equal Area Lambers projection
proj <- "+proj=laea +lat_0=0 +lon_0=-70 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs"

env <- stack(projectRaster(env, crs = proj))
r12 <- stack(projectRaster(r12, env))
r24 <- stack(projectRaster(r24, env))
r37 <- stack(projectRaster(r37, env))

# Save for use
writeRaster(env, bylayer = T, filename = paste0("data/env/ready_layers/",
                                                names(env), "_cur.tif"),
            overwrite = T)

writeRaster(r12, bylayer = T, filename = paste0("data/env/ready_layers/",
                                                names(env), "_r12.tif"),
            overwrite = T)

writeRaster(r24, bylayer = T, filename = paste0("data/env/ready_layers/",
                                                names(env), "_r24.tif"),
            overwrite = T)

writeRaster(r37, bylayer = T, filename = paste0("data/env/ready_layers/",
                                                names(env), "_r37.tif"),
            overwrite = T)

# Load distance to coast layer and scale
dist <- raster("data/env/crop_layers/distcoast.tif")
dist <- projectRaster(dist, env)
dist <- scale(dist)

# Save scaled version
writeRaster(dist, "data/env/ready_layers/distcoast.tif", overwrite = T)

# Include in the env stack
names(dist) <- "dist"
env <- stack(env, dist)

# Reproject study area
starea <- spTransform(starea, CRS(proj))

# Simplify study area
starea <- gSimplify(gBuffer(starea, width=0.3*111), tol=20)



# Prepare INLA mesh ----

# Degree to km multiplier
dgk <- 111

# Creates mesh and plot
mesh <- inla.mesh.2d(
  boundary = starea,
  max.edge = c(0.5, 3) * dgk,
  cutoff = 0.1 * dgk,
  crs = crs(env),
  offset = c(0.1, 6) * dgk
)

ggplot()+gg(mesh)+coord_fixed()+spatt;mesh$n

# Get integration points
ips <- ipoints(starea, mesh)

ggplot()+gg(mesh)+gg(ips)+coord_fixed()+spatt

# Get data for integration points using IDW ----
getd <- function(rast, ip){
        rast <- extend(rast, (extent(min(mesh$loc[,1]),
                                     max(mesh$loc[,1]),
                                     min(mesh$loc[,2]),
                                     max(mesh$loc[,2])))+
                               c(-200, 200, -200, 200))
        
        epts <- extract(rast, ip)
        epts <- data.frame(epts, coordinates(ip))
        
        tofill <- epts[is.na(epts[,1]),]
        
        epts <- data.frame(rasterToPoints(rast))
        
        coordinates(epts) <- ~x+y
        coordinates(tofill) <- ~x+y
        
        for (i in 1:nlayers(rast)) {
                
                mod <- gstat(formula = as.formula(paste(names(epts)[i],"~ 1")),
                             data = epts, nmax = 12)
                
                pred <- predict(mod, tofill)
                
                rast[[i]][cellFromXY(rast[[i]], coordinates(tofill))] <- pred$var1.pred
        }
        
        rast
}

env.e <- getd(env, ip = ips)

# Convert to SpatialPixels*
env.e <- as(env.e, "SpatialPixelsDataFrame")


# Save everything ----
lgcp.data <- list(
        # Environmental layers expanded data
        env.e = env.e,
        # Mesh
        mesh = mesh,
        # Study area shape
        starea = starea
)

saveRDS(lgcp.data, file = "data/lgcp_data.rds")
