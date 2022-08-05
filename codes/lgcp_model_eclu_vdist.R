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
# Plotting and utilities
library(ggplot2)
library(fs)
# Settings
set.seed(2932) # Replicability of sampling
spatt <- theme_classic()+theme(axis.line = element_blank()) # Theme for better ploting



# Load environmental and species data ----
# These are prepared in the code lgcp_prepare_data.R
list2env(readRDS('data/lgcp_data.rds'),globalenv())

# Load environmental layers
env <- stack(list.files("data/env/ready_layers", pattern = "_cur", full.names = T))
r12 <- stack(list.files("data/env/ready_layers", pattern = "_r12", full.names = T))
r24 <- stack(list.files("data/env/ready_layers", pattern = "_r24", full.names = T))
r37 <- stack(list.files("data/env/ready_layers", pattern = "_r37", full.names = T))

# Include the distance to coast layer
env <- stack(env, raster("data/env/ready_layers/distcoast.tif"))
r12[[8]] <- r24[[8]] <- r37[[8]] <- env[[8]]

# Change names
names(env) <- c("chl", "coldm", "ph", "pho", "sal", "sst", "warmm", "dist")
names(r12) <- names(r24) <- names(r37) <- names(env)

# Load species data
species <- "eclu"

pts <- read.csv(paste0("data/", species, "/", species, "_cell.csv"))[,1:2]
pts <- SpatialPoints(data.frame(x = pts[,1], y = pts[,2]), 
                     proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))

# Reproject species data to the Equal Area Lambers projection
proj <- "+proj=laea +lat_0=0 +lon_0=-70 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs"
pts <- spTransform(pts, CRS(proj))

# Plot INLA mesh
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
size <- diff(bbox(starea)[1,])
range0 <- as.numeric(round(size/5))
sigma <- 0.5

# Stationary model for comparison (if wanted)
# spde <- inla.spde2.pcmatern(mesh,
#                             prior.range = c(range0, 0.1),
#                             prior.sigma = c(sigma, 0.01))

# Barrier model
b.model <- inla.barrier.pcmatern(mesh,
                                 barrier.triangles = barrier.triangles,
                                 prior.range = c(range0, 0.1),
                                 prior.sigma = c(sigma, 0.01))

source("functions/barrier_model_plot.R")
# par(mfrow = c(1,2), mar = c(5, 5, 4, 5) + 0.1)
plot.bmodel(pts@coords[10,1:2], mesh = mesh, spde = b.model,
            areapol = poly.barrier, crs = CRS(proj), range = range0, msd = 0.5)
# To plot the stationary version:
# plot.bmodel(pts@coords[10,1:2], mesh = mesh, spde = spde,
#             areapol = poly.barrier, spmode = T, range = range0)



# Create models formulas ----

# Creates 1D spde for modeling temperature as a non-linear response
# Two models are created - one in the range of SST and Coldest Month
# and another in the range of Warmest Month

# Get values
min.sst <- data.frame(current = minValue(env[[c("sst", "coldm", "warmm")]]),
                      r12 = minValue(r12[[c("sst", "coldm", "warmm")]]),
                      r24 = minValue(r24[[c("sst", "coldm", "warmm")]]),
                      r37 = minValue(r37[[c("sst", "coldm", "warmm")]]))
round(apply(min.sst, 1, min), 2)

max.sst <- data.frame(current = maxValue(env[[c("sst", "coldm", "warmm")]]),
                      r12 = maxValue(r12[[c("sst", "coldm", "warmm")]]),
                      r24 = maxValue(r24[[c("sst", "coldm", "warmm")]]),
                      r37 = maxValue(r37[[c("sst", "coldm", "warmm")]]))
round(apply(max.sst, 1, max), 2)

# SST
knots.st <- seq(-2.6, 1.6, length = 40)
d1mesh.st <- inla.mesh.1d(knots.st, interval = c(-2.6, 1.6), degree = 2)

d1spde.st <- inla.spde2.pcmatern(d1mesh.st,
                              prior.range = c(diff(c(-2.6, 1.6)), NA),
                              prior.sigma = c(1, 0.1),
                              constr = T)

# Coldest Month
knots.cm <- seq(-2.5, 1.4, length = 40)
d1mesh.cm <- inla.mesh.1d(knots.cm, interval = c(-2.5, 1.4), degree = 2)

d1spde.cm <- inla.spde2.pcmatern(d1mesh.cm,
                              prior.range = c(diff(c(-2.5, 1.4)), NA),
                              prior.sigma = c(1, 0.1),
                              constr = T)

# Warmest Month
knots.wm <- seq(-3.5, 1.9, length = 40)
d1mesh.wm <- inla.mesh.1d(knots.wm, interval = c(-3.5, 1.9), degree = 2)

d1spde.wm <- inla.spde2.pcmatern(d1mesh.wm,
                              prior.range = c(diff(c(-3.5, 1.9)), NA),
                              prior.sigma = c(1, 0.1),
                              constr = T)


##### Write formulas ----

# We will test 3 different model.
# All of them contains at least temperature + salinity + distance from coast,
# what we call "base components"
# base components + pH --> this model considers acidification
# base components + chl-a + phosphate --> this model considers productivity/light
# base components + pH + chlorophyll-a + phosphate --> full model.

# We also have three different variables representing sea temperature. To chose
# which one is most suitable for the species distribution model, we run the most
# basic model (base components) with each one and evaluate the response, summary, WAIC,
# and judge which one to use.
# mean SST
# mean temperature of the Warmest Month
# mean temperature of the Coldest Month

# Finally, the models [with the selected temperature component] will be tested
# with the inclusion of a Spatial component
# (an SPDE) considering a barrier model. This yelds a total of 6 testable models.
# (NOTE: in total we run 9 models, because 3 are to chose the temperature variable)

forms <- list()

base.f <- coordinates ~ 
    Intercept(1, mean.linear = log(length(pts)/gArea(starea)), prec.linear = 1) + 
    sal(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
    dist(env.e, model = "linear", mean.linear = 0, prec.linear = 1)
        
forms[[1]] <- update(base.f, ~ . + sst(env.e, model = d1spde.st))

forms[[2]] <- update(base.f, ~ . + coldm(env.e, model = d1spde.cm))

forms[[3]] <- update(base.f, ~ . + warmm(env.e, model = d1spde.wm))

# After judgement, for this species we chose: SST
forms[[4]] <- update(forms[[1]], ~ . +
                         ph(env.e, model = "linear", mean.linear = 0, prec.linear = 1))

forms[[5]] <- update(forms[[1]], ~ . +
                         chl(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
                         pho(env.e, model = "linear", mean.linear = 0, prec.linear = 1))

forms[[6]] <- update(forms[[1]], ~ . +
                         ph(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
                         chl(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
                         pho(env.e, model = "linear", mean.linear = 0, prec.linear = 1))

# Include models with SPDE
forms[[7]] <- update(forms[[4]], ~ . + spatial(coordinates, model = b.model,
                                               mapper = bru_mapper(mesh)))

forms[[8]] <- update(forms[[5]], ~ . + spatial(coordinates, model = b.model,
                                               mapper = bru_mapper(mesh)))

forms[[9]] <- update(forms[[6]], ~ . + spatial(coordinates, model = b.model,
                                               mapper = bru_mapper(mesh)))


# Run models ----
# Initiate a list to hold models
m <- list()

# Creates a data.frame to hold WAIC, DIC and CPO (summary) results
m.metrics <- data.frame(waic = rep(NA, length(forms)), dic = NA)
m.lambda <- list()

# Run models in loop
for (i in 1:length(forms)) {
        
        cat("Runing model", i, "\n")
        print(forms[[i]])
        
        m[[i]] <- lgcp(forms[[i]], pts,
                       ips = ips,
                       #domain = list(coordinates = mesh),
                       options = list(control.inla = list(int.strategy = "auto",
                                                          strategy = "simplified.laplace"),
                                      inla.mode = "classic"
                                      # For verbosity, uncoment those lines:
                                      # , bru_verbose = TRUE,
                                      # verbose = TRUE
                                      ))
        
        m.metrics$waic[i] <- m[[i]]$waic$waic
        m.metrics$dic[i] <- m[[i]]$dic$dic

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
plot.random(m[[4]])
plot.post(m[[4]])

# Get summary
summary(m[[4]])

# Plot spatial effect (only for the ones with a spatial component!)
plot.spatial(m[[7]], mesh)

# Compare models
deltaIC(m[[1]], m[[2]], m[[3]], criterion = "WAIC")

deltaIC(m[[4]], m[[5]], m[[6]],
        m[[7]], m[[8]], m[[9]], criterion = "WAIC")


# Predict each component effect spatially
# Get a base PixelsDataFrame with same reso as env layers
df <- as(env$ph, "SpatialPixelsDataFrame")
df <- df[,-1]

pred.comp <- predict(m[[4]], data = df,
                 formula = ~ list(
                         # Change those according to mode and
                         # combinations
                         sst = sst,
                         sal = sal,
                         dist = dist,
                         ph = ph,
                         #spatial = spatial,
                         int = exp(sst + sal + dist + ph +
                                      # spatial +
                                       Intercept)
                 ))

ggplot()+
        gg(pred.comp$ph, aes(fill = mean))+ # change here according to layer
        coord_equal()

# We select the model 4 as the best one, considering WAIC and exploratory analysis
sel.model <- c(4, 7) # we also project the one with the spatial component



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

# Predict model without spatial component
pred.comp <- predict(m[[sel.model[1]]], data = df,
                     formula = ~ list(
                             sst = sst,
                             sal = sal,
                             dist = dist,
                             ph = ph,
                             combl = sst + sal + dist + ph + Intercept,
                             combl_exp = exp(sst + sal + dist + ph + Intercept)),
                     n.samples = 1000, seed = 2932)

# save
sapply(names(pred.comp), save.rast, model = sel.model[1])

# Predict model with spatial component
pred.comp <- predict(m[[sel.model[2]]], data = df,
                     formula = ~ list(
                             sst = sst,
                             sal = sal,
                             ph = ph,
                             dist = dist,
                             spatial = spatial,
                             combl = sst + sal + dist + ph + Intercept,
                             combl_exp = exp(sst + sal + dist + ph + Intercept)),
                     n.samples = 1000, seed = 2932)

# save
sapply(names(pred.comp), save.rast, model = sel.model[2])



# Predictions ----
# inlabru gets the environmental layer from the environment
# thus, to predict to different scenarios we need to change the environmental
# layers object (i.e. 'env.e') to the new object for which we want the prediction

# we do the same for the model without and with the SPDE effect
for (i in sel.model) {
    # Get formula
    pred.f <- as.formula(paste0(
        "~ exp(",
        paste(if (is.null(names(
            m[[i]]$summary.random))) {
            m[[i]]$names.fixed
        } else{
            c(m[[i]]$names.fixed,
              names(m[[i]]$summary.random))
        }
        , collapse = "+"), ")"
    ))
    
    # Predict to current scenario
    env.e <- as(env, "SpatialPixelsDataFrame")
    
    pred.cur <- predict(m[[i]],
                        data = df,
                        formula = pred.f,
                        n.samples = 1000, seed = 2932)
    ggplot() +
        gg(pred.cur, aes(fill = mean)) + coord_equal()
    
    # Predict to SSP 1
    env.e <- as(r12, "SpatialPixelsDataFrame")
    
    pred.r12 <- predict(m[[i]],
                        data = df,
                        formula = pred.f,
                        n.samples = 1000, seed = 2932)
    ggplot() +
        gg(pred.r12, aes(fill = mean)) + coord_equal()
    
    # Predict to SSP 2
    env.e <- as(r24, "SpatialPixelsDataFrame")
    
    pred.r24 <- predict(m[[i]],
                        data = df,
                        formula = pred.f,
                        n.samples = 1000, seed = 2932)
    ggplot() +
        gg(pred.r24, aes(fill = mean)) + coord_equal()
    
    # Predict to SSP 3
    env.e <- as(r37, "SpatialPixelsDataFrame")
    
    pred.r37 <- predict(m[[i]],
                        data = df,
                        formula = pred.f,
                        n.samples = 1000, seed = 2932)
    ggplot() +
        gg(pred.r37, aes(fill = mean)) + coord_equal()
    
    
    # Convert all to raster (easier handling)
    pred.cur <- as(pred.cur, "RasterStack")
    pred.r12 <- as(pred.r12, "RasterStack")
    pred.r24 <- as(pred.r24, "RasterStack")
    pred.r37 <- as(pred.r37, "RasterStack")
    
    to.remove <- names(pred.cur)[!names(pred.cur) %in% 
                                     c("mean", "sd", "q0.025", "q0.975")]
    
    pred.cur <- dropLayer(pred.cur, to.remove)
    pred.r12 <- dropLayer(pred.r12, to.remove)
    pred.r24 <- dropLayer(pred.r24, to.remove)
    pred.r37 <- dropLayer(pred.r37, to.remove)
    
    # Save results
    dir <- paste("results", species, "predictions", sep = "/")
    dir_create(dir)
    
    writeRaster(pred.cur, paste0(dir, "/", species, "_",
                                 names(pred.cur), "_m", i, "_current.tif"),
                bylayer = T, overwrite = T)
    writeRaster(pred.r12, paste0(dir, "/", species, "_",
                                 names(pred.r12), "_m", i, "_ssp1.tif"),
                bylayer = T, overwrite = T)
    writeRaster(pred.r24, paste0(dir, "/", species, "_",
                                 names(pred.r24), "_m", i, "_ssp2.tif"),
                bylayer = T, overwrite = T)
    writeRaster(pred.r37, paste0(dir, "/", species, "_",
                                 names(pred.r37), "_m", i, "_ssp3.tif"),
                bylayer = T, overwrite = T)
}



# Save metrics ----
dir <- paste0("results/", species)
write.csv(m.metrics, paste0(dir, "/", species, "_model_metrics.csv"),
          row.names = F)



# Save session info ----
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