#### Modelling of coral reef herbivorous invertebrates ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2022 ##

## Barrier Model under the Log-Gaussian Cox Process framework ##

## Species: Lytechinus variegatus

# We run the model for each species separately. That way we can optimize the
# process accordingly.




# Load packages and define settings ----
# Modeling
library(INLA)
library(inlabru)
# Spatial data
library(raster)
library(rgeos)
library(gstat)
# Cross-validation
library(blockCV)
library(pROC)
# Plotting and utilities
library(ggplot2)
library(fs)
# Settings
set.seed(2932) # Replicability of sampling
spatt <- theme_classic() # Theme for better ploting



# Load environmental and species data ----
source("functions/Varload.r")
lays <- c("salinitymean", "tempmean", "ph", "chlomean")
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
names(env) <- c("ph", "chl", "sal", "sst", "coldm", "warmm")
names(r37) <- names(r24) <- names(r12) <- names(env)

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

# Load species data
species <- "lyva"

pts <- read.csv(paste0("data/", species, "/", species, "_cell.csv"))[,1:2]
coordinates(pts) <- ~ decimalLongitude + decimalLatitude
crs(pts) <- crs(env)
colnames(pts@coords) <- c("x", "y")

# Reproject everything to the Equal Area Lambers projection
proj <- "+proj=laea +lat_0=0 +lon_0=-70 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs"

env <- projectRaster(env, crs = proj)
r12 <- projectRaster(r12, env)
r24 <- projectRaster(r24, env)
r37 <- projectRaster(r37, env)

starea <- spTransform(starea, CRS(proj))

pts <- spTransform(pts, CRS(proj))



# Prepare INLA mesh ----

# Degree to km multiplier
dgk <- 111

# Creates mesh and plot
mesh <- inla.mesh.2d(boundary = starea,
                     max.edge = c(1.2, 20)*dgk,
                     cutoff = 0.25*dgk, crs = crs(env))

ggplot()+gg(mesh)+coord_fixed()+spatt;mesh$n



# Prepare the barrier model ----
# Based on Bakka, H., J. Vanhatalo, J. Illian, D. Simpson, and H. Rue. 2016.
# “Accounting for Physical Barriers in Species Distribution Modeling with
# Non-Stationary Spatial Random Effects.” ArXiv preprint arXiv:1608.03787.
# Norwegian University of Science; Technology, Trondheim, Norway.
# https://haakonbakkagit.github.io/btopic128.html

# Get number of triangles of the mesh
tl <- length(mesh$graph$tv[,1])

# Create a matrix with the central coordinates of each triangle
tri.cord <- matrix(0, tl, 2)

for(i in 1:tl){
        # Take the vertex of triangles
        temp <- mesh$loc[mesh$graph$tv[i, ], ]
        
        # Compute center of each triangle
        tri.cord[i,] <- colMeans(temp)[c(1,2)]
}

# Convert to SpatialPoints
tri.cord <- SpatialPoints(tri.cord)

# Get intersection between mesh points and study area
crs(tri.cord) <- crs(starea)
intersec <- over(starea, tri.cord, returnList = T)

# Remove the study area triangles from all to obtain the barrier ones
intersec <- unlist(intersec)
barrier.triangles <- setdiff(1:tl, intersec)

# Create a barrier polygon, i.e. all the polygons composing the "islands"
poly.barrier <- inla.barrier.polygon(mesh, barrier.triangles)



# Matérn models ----
# Get range and sigma
size <- diff(range(mesh$loc[,2]))
range0 <- round(size / 2)
sigma <- 3

# Stationary model for comparison (if wanted)
spde <- inla.spde2.pcmatern(mesh,
                            prior.range = c(range0, 0.1),
                            prior.sigma = c(sigma, 0.01))

# Barrier model
b.model <- inla.barrier.pcmatern(mesh,
                                 barrier.triangles = barrier.triangles,
                                 prior.range = c(range0, 0.1),
                                 prior.sigma = c(sigma, 0.01))

source("functions/barrier_model_plot.R")
# par(mfrow = c(1,2), mar = c(5, 5, 4, 5) + 0.1)
plot.bmodel(pts@coords[10,1:2], mesh = mesh, spde = b.model,
            areapol = poly.barrier, crs = CRS(proj), range = range0, msd = 3)
# plot.bmodel(pts@coords[10,1:2], mesh = mesh, spde = spde,
#             areapol = poly.barrier, spmode = T, range = range0) # Plot SPDE



# Expand raster layers to cover mesh points ----

# Get integration points
source("functions/get_weights.R")
ips <- get.weights(mesh, starea)

getd <- function(rast, ip, pts){
        rast <- extend(rast, (extent(min(mesh$loc[,1]),
                                    max(mesh$loc[,1]),
                                    min(mesh$loc[,2]),
                                    max(mesh$loc[,2])))+
                               c(-200, 200, -200, 200))

        epts <- rbind(extract(rast, ip),
                      extract(rast, pts))
        epts <- data.frame(epts, rbind(coordinates(ip),
                                       coordinates(pts)))
        
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

env.e <- getd(env, ip = ips, pts = pts)

# Convert to SpatialPixels*
env.e <- as(env.e, "SpatialPixelsDataFrame")



# Create models formulas ----

# Creates 1D spde for modeling temperature as a non-linear response
# Two models are created - one in the range of SST and Coldest Month
# and another in the range of Warmest Month

# Get values
min.sst <- data.frame(current = minValue(env[[4:6]]),
                      r12 = minValue(r12[[4:6]]),
                      r24 = minValue(r24[[4:6]]),
                      r37 = minValue(r37[[4:6]]))
round(apply(min.sst, 1, min), 2)

max.sst <- data.frame(current = maxValue(env[[4:6]]),
                      r12 = maxValue(r12[[4:6]]),
                      r24 = maxValue(r24[[4:6]]),
                      r37 = maxValue(r37[[4:6]]))
round(apply(max.sst, 1, max), 2)

# SST and Coldest Month
knots <- seq(-2.6, 1.6, length = 25)
d1mesh <- inla.mesh.1d(knots, interval = c(-2.6, 1.6), degree = 2,
                       boundary = "free")

d1spde <- inla.spde2.pcmatern(d1mesh,
                              prior.range = c(3.2, NA),
                              prior.sigma = c(1, 0.1),
                              constr = T)

# Warmest Month
knots.wm <- seq(-3.5, 1.9, length = 25)
d1mesh.wm <- inla.mesh.1d(knots.wm, interval = c(-3.5, 1.9), degree = 2,
                          boundary = "free")

d1spde.wm <- inla.spde2.pcmatern(d1mesh.wm,
                              prior.range = c(3.8, NA),
                              prior.sigma = c(1, 0.1),
                              constr = T)


##### Write formulas ----

# We will test 4 different base models:
# salinity + temperature --> at least those two variables should affect the dist.
# sal + temp + pH --> this model considers acidification
# sal + temp + chl --> this model considers productivity/light
# sal + temp + pH + chl --> full model.

# We also have three different variables representing sea temperature. To chose
# which one is most suitable for the species distribution model, we run the most
# basic model (sal + temp) with each one and evaluate the response, summary, WAIC,
# and judge which one to use.
# mean SST
# mean temperature of the Warmest Month
# mean temperature of the Coldest Month

# Finally, the models [with the selected temperature component] will be tested
# with the inclusion of a Spatial component
# (an SPDE) considering a barrier model. This yelds a total of 10 models.

forms <- list()

base.f <- coordinates ~ Intercept(1) + sal(env.e, model = "linear")
        
forms[[1]] <- update(base.f, ~ . + sst(env.e, model = d1spde,
                                       mapper = bru_mapper(d1mesh, indexed = T)))

forms[[2]] <- update(base.f, ~ . + coldm(env.e, model = d1spde,
                                         mapper = bru_mapper(d1mesh, indexed = T)))

forms[[3]] <- update(base.f, ~ . + warmm(env.e, model = d1spde.wm,
                                         mapper = bru_mapper(d1mesh.wm, indexed = T)))

# After judgement, for this species we chose: WARMEST MONTH
forms[[4]] <- update(forms[[3]], ~ . + ph(env.e, model = "linear"))

forms[[5]] <- update(forms[[3]], ~ . + chl(env.e, model = "linear"))

forms[[6]] <- update(forms[[3]], ~ . + ph(env.e, model = "linear") 
                     + chl(env.e, model = "linear"))

# Include models with SPDE
forms[[7]] <- update(forms[[3]], ~ . + spatial(coordinates, model = b.model,
                                               mapper = bru_mapper(mesh)))

forms[[8]] <- update(forms[[4]], ~ . + spatial(coordinates, model = b.model,
                                               mapper = bru_mapper(mesh)))

forms[[9]] <- update(forms[[5]], ~ . + spatial(coordinates, model = b.model,
                                               mapper = bru_mapper(mesh)))

forms[[10]] <- update(forms[[6]], ~ . + spatial(coordinates, model = b.model,
                                               mapper = bru_mapper(mesh)))



# Run models ----
# Initiate a list to hold models
m <- list()

# Creates a data.frame to hold WAIC, DIC and CPO (summary) results
m.metrics <- data.frame(waic = rep(NA, length(forms)), dic = NA, cpo = NA)
m.lambda <- list()

# Run models in loop
for (i in 1:length(forms)) {
        
        cat("Runing model", i, "\n")
        print(forms[[i]])
        
        m[[i]] <- lgcp(forms[[i]], pts,
                       ips = ips,
                       domain = list(coordinates = mesh),
                       options = list(control.inla = list(int.strategy = "eb"),
                                      control.compute(cpo = T),
                                      inla.mode = "classic"
                                      # For verbosity, uncoment those lines:
                                      # , bru_verbose = TRUE,
                                      # verbose = TRUE
                                      ))
        
        m.metrics$waic[i] <- m[[i]]$waic$waic
        m.metrics$dic[i] <- m[[i]]$dic$dic
        m.metrics$cpo[i] <- -sum(log(m[[i]]$cpo$cpo))
        
        lambda.exp <- predict(m[[i]], ips,
                              as.formula(paste0(
                                      "~ sum(weight * exp(",
                                      paste(if (is.null(names(
                                              m[[i]]$summary.random))) {
                                              m[[i]]$names.fixed
                                      } else{
                                              c(m[[i]]$names.fixed,
                                                names(m[[i]]$summary.random))
                                      }
                                      , collapse = "+"), "))"
                              )))
        cat("Number of points:", length(pts), "--- Lambda:\n")
        print(m.lambda[[i]] <- lambda.exp)
        cat("\n ====================== \n")

}


# Plot effects (change numbers according to model)
source("functions/plot_inla.R")
plot.random(m[[1]])
plot.post(m[[1]])

# Get summary
summary(m[[1]])

# Plot spatial effect (only for the ones with a spatial component!)
plot.spatial(m[[7]], mesh)

# Compare models
deltaIC(m[[3]], m[[4]], m[[5]], m[[6]],
        m[[7]], m[[8]], m[[9]], m[[10]], criterion = "WAIC")
# Now disconsidering the SPDE
deltaIC(m[[3]], m[[4]], m[[5]], m[[6]], criterion = "WAIC")

# Predict each component effect spatially
df <- pixels(mesh, mask = starea, nx = 400, ny = 400)
pred.comp <- predict(m[[4]], data = df,
                 formula = ~ list(
                         # Change those according to mode and
                         # combinations
                         warmm = warmm,
                         sal = sal,
                         ph = ph,
                         int = exp(warmm + sal + ph + Intercept)
                 ))

ggplot()+
        gg(pred.comp$ph)+ # change here according to layer
        coord_equal()

# We select the model 4 as the best one, considering WAIC and exploratory analysis
sel.model <- 4



### Cross-validation ----

# We use a spatial block cross validation, which leaves entire sections of
# the environment out in each run of the model.
# As a first step, we divide the area in blocks, using package blockCV

spat.range <- spatialAutoRange(env, showPlots = F) # get ideal range for the blocks

# Get coordinates and PA
dat <- data.frame(rbind(
    cbind(ips@coords, pa = 0),
    cbind(pts@coords, pa = 1)
))
coordinates(dat) <- ~ x + y

# Divide in blocks
blocks <- spatialBlock(dat, species = "pa",
                       theRange = spat.range$range)
# Get the ID
blocks <- blocks$foldID

b.ips <- blocks[1:length(ips)]
b.pts <- blocks[(length(ips)+1):length(blocks)]

# For which model(s) run the CV?
selmodels <- c(sel.model, (sel.model+4)) # Here we run for the chosen model and 
# it's equivalent with spatial (SPDE) component

cvm.preds <- list()

for (j in 1:length(selmodels)) {
    
    cat("Runing model", selmodels[j], "CV. \n")
    
    # Creates a vector to hold predicted intesities
    pred.intensity.ip <- c()
    pred.intensity.pts <- c()
    
    pred.rast <- stack()
    
    cv.f <- as.formula(paste0(
        "~ exp(",
        paste(if (is.null(names(
            m[[selmodels[j]]]$summary.random))) {
            m[[selmodels[j]]]$names.fixed
        } else{
            c(m[[selmodels[j]]]$names.fixed,
              names(m[[selmodels[j]]]$summary.random))
        }
        , collapse = "+"), ")"
    ))
    
    df.cv <- as(env$ph, "SpatialPixelsDataFrame")
    df.cv <- df.cv[,-1]
    
    for (i in 1:5) {
        cat("Block", i, "\n")
        
        cv.ips <- ips[which(b.ips != i),]
        cv.pts <- pts[which(b.pts != i),]
        
        cv.m <- lgcp(forms[[selmodels[j]]], cv.pts,
                     ips = cv.ips,
                     domain = list(coordinates = mesh),
                     options = list(control.inla = list(int.strategy = "eb"),
                                    inla.mode = "classic"))
        
        pred.ip <- predict(cv.m, data = ips[which(b.ips == i),], cv.f)
        pred.pts <- predict(cv.m, data = pts[which(b.pts == i),], cv.f)
        
        pred.intensity.ip[which(b.ips == i)] <- pred.ip$mean
        pred.intensity.pts[which(b.pts == i)] <- pred.pts$mean
        
        pred.r <- predict(cv.m, data = df.cv, cv.f)
        pred.r <- as(pred.r, "RasterStack")
        
        pred.rast <- stack(pred.rast, pred.r$mean)
    }
    
    full.intensity.ip <- predict(m[[selmodels[j]]], data = ips, cv.f)
    
    full.intensity.ip <- full.intensity.ip$mean
    
    full.intensity.pts <- predict(m[[selmodels[j]]], data = pts, cv.f)
    
    full.intensity.pts <- full.intensity.pts$mean
    
    original <- c(full.intensity.ip, full.intensity.pts)
    predicted <- c(pred.intensity.ip, pred.intensity.pts)
    
    #### aparentemente funcionando! testar loocv
    rmse <- sqrt(mean((original - predicted)^2))
    
    pred.rast <- mean(pred.rast)
    
    pred.rast.full <- predict(m[[selmodels[j]]], data = df.cv, cv.f)
    pred.rast.full <- as(pred.rast.full, "RasterStack")
    
    im <- dismo::nicheOverlap(pred.rast, pred.rast.full$mean, stat = "I")
    dm <- dismo::nicheOverlap(pred.rast, pred.rast.full$mean, stat = "D")
    
    cvm.preds[[j]] <- list(
        "rmse" = rmse,
        "pred.rast" = pred.rast,
        "pred.rsat.full" = pred.rast.full$mean,
        "niche.overlap" = c("I" = im, "D" = dm)
    )
}

# See results
cvm.preds[[1]]$rmse;cvm.preds[[2]]$rmse
cvm.preds[[1]]$niche.overlap;cvm.preds[[2]]$niche.overlap



# Save individual components effects for the chosen model ----
# This was already done after running the models as a exploratory tool
# Now we run again and save the results

# Create a function to save results
dir <- paste("results", species, "effects", sep = "/")
dir_create(dir)

save.rast <- function(x, model){
        r <- as(pred.comp[[x]], "RasterStack")
        writeRaster(r, paste0(dir, "/", species, "_m", model, "_",
                              x, "_", names(r), "_effect.tif"),
                    bylayer = T, overwrite = T)
        return(paste(x, "saved."))
}

# Get a base PixelsDataFrame with same reso as env layers
df <- as(env$ph, "SpatialPixelsDataFrame")
df <- df[,-1]

# Predict model without spatial component
pred.comp <- predict(m[[4]], data = df,
                     formula = ~ list(
                             warmm = warmm,
                             sal = sal,
                             ph = ph))

# save
sapply(names(pred.comp), save.rast, model = 4)

# Predict model with spatial component
pred.comp <- predict(m[[8]], data = df,
                     formula = ~ list(
                             warmm = warmm,
                             sal = sal,
                             ph = ph,
                             spatial = spatial))

# save
sapply(names(pred.comp), save.rast, model = 8)



# Predictions ----
# inlabru gets the environmental layer from the environment
# thus, to predict to different scenarios we need to change the environmental
# layers object (i.e. 'env.e') to the new object for which we want the prediction

# Get formula
pred.f <- as.formula(paste0(
        "~ exp(",
        paste(if (is.null(names(
                m[[sel.model]]$summary.random))) {
                m[[sel.model]]$names.fixed
        } else{
                c(m[[sel.model]]$names.fixed,
                  names(m[[sel.model]]$summary.random))
        }
        , collapse = "+"), ")"
))

# Predict to current scenario
env.e <- as(env, "SpatialPixelsDataFrame")

pred.cur <- predict(m[[sel.model]],
                    data = df,
                    formula = pred.f)
ggplot() +
        gg(pred.cur, aes(fill = mean)) + coord_equal()

# Predict to SSP 1
env.e <- as(r12, "SpatialPixelsDataFrame")

pred.r12 <- predict(m[[sel.model]],
                    data = df,
                    formula = pred.f)
ggplot() +
        gg(pred.r12, aes(fill = mean)) + coord_equal()

# Predict to SSP 2
env.e <- as(r24, "SpatialPixelsDataFrame")

pred.r24 <- predict(m[[sel.model]],
                    data = df,
                    formula = pred.f)
ggplot() +
        gg(pred.r24, aes(fill = mean)) + coord_equal()

# Predict to SSP 3
env.e <- as(r37, "SpatialPixelsDataFrame")

pred.r37 <- predict(m[[sel.model]],
                    data = df,
                    formula = pred.f)
ggplot() +
        gg(pred.r37, aes(fill = mean)) + coord_equal()


# Convert all to raster (easier handling)
pred.cur <- as(pred.cur, "RasterStack")
pred.r12 <- as(pred.r12, "RasterStack")
pred.r24 <- as(pred.r24, "RasterStack")
pred.r37 <- as(pred.r37, "RasterStack")

to.remove <- names(pred.cur)[!names(pred.cur) %in% c("mean", "sd", "q0.025", "q0.975")]

pred.cur <- dropLayer(pred.cur, to.remove)
pred.r12 <- dropLayer(pred.r12, to.remove)
pred.r24 <- dropLayer(pred.r24, to.remove)
pred.r37 <- dropLayer(pred.r37, to.remove)

# Save results
dir <- paste("results", species, "predictions", sep = "/")
dir_create(dir)

writeRaster(pred.cur, paste0(dir, "/", species, "_", names(pred.cur), "_current.tif"),
            bylayer = T, overwrite = T)
writeRaster(pred.r12, paste0(dir, "/", species, "_", names(pred.r12), "_ssp1.tif"),
            bylayer = T, overwrite = T)
writeRaster(pred.r24, paste0(dir, "/", species, "_", names(pred.r24), "_ssp2.tif"),
            bylayer = T, overwrite = T)
writeRaster(pred.r37, paste0(dir, "/", species, "_", names(pred.r37), "_ssp3.tif"),
            bylayer = T, overwrite = T)

# Get differences
dif.r12 <- pred.r12 - pred.cur
dif.r24 <- pred.r24 - pred.cur
dif.r37 <- pred.r37 - pred.cur

names(dif.r12) <- names(dif.r24) <- names(dif.r37) <- names(pred.cur)

# Save results
writeRaster(dif.r12, paste0(dir, "/", species, "_", names(dif.r12), "_dif_ssp1.tif"),
            bylayer = T, overwrite = T)
writeRaster(dif.r24, paste0(dir, "/", species, "_", names(dif.r24), "_dif_ssp2.tif"),
            bylayer = T, overwrite = T)
writeRaster(dif.r37, paste0(dir, "/", species, "_", names(dif.r37), "_dif_ssp3.tif"),
            bylayer = T, overwrite = T)



# Save other results ----
dir <- paste("results", species, sep = "/")
saveRDS(cvm.preds, file = paste0(dir, "/", species, "_cvresults.rds"))

# Save session info
dir.info <- paste0(dir, "/", species, "_sessioninfo.txt")
cat("\n ============= \n \n Research info \n \n
    Author: Silas C. Principe | silasprincipe@usp.br \n
    Co-authors: Tito M.C. Lotufo | André L. Acosta \n
    Modelling of reef herbivorous invertebrates (sea-urchins) \n
    This work is part of a PhD project being held at the Oceanographic Institute - USP",
    file = dir.info)
cat("\n ============= \n \n Session info \n \n", file = dir.info, append = T)
write(capture.output(sessionInfo()), file = dir.info, append = T)
cat("\n ============= \n \n INLA info \n \n", file = dir.info, append = T)
write(capture.output(inla.version()), file = dir.info, append = T)

## END OF CODE