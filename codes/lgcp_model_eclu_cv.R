#### Modelling of coral reef herbivorous invertebrates ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2022 ##

## Barrier Model under the Log-Gaussian Cox Process framework ##

## Species: Echinometra lucunter

# We run the model for each species separately. That way we can optimize the
# process accordingly.




# Load packages and define settings ----
# Modeling
library(INLA)
library(inlabru)
# Spatial data
library(raster)
library(rgeos)
# Plotting and utilities
library(ggplot2)
library(patchwork)
library(fs)

# Settings
set.seed(2932) # Replicability of sampling
spatt <- theme_classic()+theme(axis.line = element_blank()) # Theme for better ploting
# Integration strategy (change here for "eb" for a fast approximation)
intest <- "auto"
# Number of samples to calculate posteriors (decrease to a fast approx.)
nsamp <- 1000
# RGDAL: avoid saving rasters with auxiliary files
rgdal::setCPLConfigOption("GDAL_PAM_ENABLED", "FALSE")



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

pts <- read.csv(paste0("data/", species, "/", species, "_filt.csv"))
proj <- "+proj=laea +lat_0=0 +lon_0=-70 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs"
pts <- SpatialPoints(pts, proj4string = CRS(proj))

# Plot INLA mesh
ggplot()+gg(mesh)+coord_fixed()+spatt;mesh$n

# Get integration points
ips <- ipoints(starea, mesh)



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
range0 <- as.numeric(round(size/10))
sigma <- 0.1

# Stationary model for comparison (if wanted)
# spde <- inla.spde2.pcmatern(mesh,
#                             prior.range = c(range0, 0.1),
#                             prior.sigma = c(sigma, 0.01))

# Barrier model
b.model <- inla.barrier.pcmatern(mesh,
                                 barrier.triangles = barrier.triangles,
                                 prior.range = c(range0, 0.5),
                                 prior.sigma = c(sigma, 0.01))

source("functions/barrier_model_plot.R")
# par(mfrow = c(1,2), mar = c(5, 5, 4, 5) + 0.1)
plot.bmodel(pts@coords[10,1:2], mesh = mesh, spde = b.model,
            areapol = poly.barrier, crs = CRS(proj), range = range0, msd = sigma)
# To plot the stationary version:
# plot.bmodel(pts@coords[10,1:2], mesh = mesh, spde = spde,
#             areapol = poly.barrier, spmode = T, range = range0)



# Create models formulas ----

# Creates 1D spde for modeling temperature as a non-linear response

# Get values for reference
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
# Leave an "outer bound" of ~0.5 unit for each side
knots.st <- seq(-3, 2, length = 30)
d1mesh.st <- inla.mesh.1d(knots.st, degree = 2,
                          boundary = "free")

d1spde.st <- inla.spde2.pcmatern(d1mesh.st,
                              prior.range = c(diff(c(-3, 2)), NA),
                              prior.sigma = c(0.5, 0.5),
                              constr = T)

# Coldest Month can use the same as SST
# knots.cm <- seq(-2.5, 1.4, length = 40)
# d1mesh.cm <- inla.mesh.1d(knots.cm, interval = c(-2.5, 1.4), degree = 2)
# 
# d1spde.cm <- inla.spde2.pcmatern(d1mesh.cm,
#                               prior.range = c(diff(c(-2.5, 1.4)), NA),
#                               prior.sigma = c(1, 0.1),
#                               constr = T)

# Warmest Month
knots.wm <- seq(-4, 2.5, length = 30)
d1mesh.wm <- inla.mesh.1d(knots.wm, degree = 2,
                          boundary = "free")

d1spde.wm <- inla.spde2.pcmatern(d1mesh.wm,
                              prior.range = c(diff(c(-4, 2.5)), NA),
                              prior.sigma = c(0.5, 0.5),
                              constr = T)



##### Write formulas ----

# We will test 3 different model.
# All of them contains at least temperature + salinity + distance from coast,
# what we call "base components"
# The base components also include the spatial component (SPDE)
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

# This yields a total of 6 testable models

forms <- list()

base.f <- coordinates ~ 
    Intercept(1, mean.linear = log(length(pts)/gArea(starea)), prec.linear = 1) + 
    sal(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
    dist(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
    spatial(coordinates, model = b.model, mapper = bru_mapper(mesh))
        
forms[[1]] <- update(base.f, ~ . + sst(env.e, model = d1spde.st))

forms[[2]] <- update(base.f, ~ . + coldm(env.e, model = d1spde.st))

forms[[3]] <- update(base.f, ~ . + warmm(env.e, model = d1spde.wm))

# After judgement, for this species we chose: mean of coldest month
chosen.bm <- forms[[2]]

forms[[4]] <- update(chosen.bm, ~ . +
                         ph(env.e, model = "linear", mean.linear = 0, prec.linear = 1))

forms[[5]] <- update(chosen.bm, ~ . +
                         chl(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
                         pho(env.e, model = "linear", mean.linear = 0, prec.linear = 1))

forms[[6]] <- update(chosen.bm, ~ . +
                         ph(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
                         chl(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
                         pho(env.e, model = "linear", mean.linear = 0, prec.linear = 1))



# Run models ----
# Initiate a list to hold models
m <- list()

# Creates a data.frame to hold WAIC and DIC results
m.metrics <- data.frame(waic = rep(NA, 6), dic = NA)
m.lambda <- list()

# Run models in loop
for (i in 1:length(forms)) {
        
        cat("Runing model", i, "\n")
        print(forms[[i]])
        
        m[[i]] <- lgcp(forms[[i]], pts,
                       ips = ips,
                       options = list(control.inla = list(int.strategy = intest,
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

# Plot spatial effect
plot.spatial(m[[4]], mesh)

# Compare models
# Just for when chosing between the three SST components
#deltaIC(m[[1]], m[[2]], m[[3]], criterion = "WAIC")

deltaIC(m[[4]], m[[5]], m[[6]], criterion = "WAIC")


# Predict each component effect spatially
# Get a base PixelsDataFrame with same reso as env layers
df <- as(env$ph, "SpatialPixelsDataFrame")
df <- df[,-1]

pred.comp <- predict(m[[4]], data = df,
                 formula = ~ list(
                         ### Change those according to mode and combinations
                         # Individual components
                         coldm = coldm,
                         sal = sal,
                         dist = dist,
                         ph = ph,
                         spatial = spatial,
                         # Integrated (exp(lambda))
                         int = exp(coldm + sal + ph + dist + spatial + Intercept),
                         # Linear predictor scale
                         lin = coldm + sal + ph + dist + spatial + Intercept,
                         # Contrast models (excluding spatial component)
                         int_cont = exp(coldm + sal + ph + dist + Intercept),
                         lin_cont = coldm + sal + ph + dist + Intercept
                 ))

# Get a better color palete for ploting:
pc <- function(component = NULL, var = "mean", pred.ob = pred.comp){
  sca <- function(...){
    scale_fill_gradientn(
      colours = rev(c("#A84C00", "#D97D27", "#F5BD44", "#FFD561", "#FFF291",
                      "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695")),
      limits = range(...))
  }
  
  if (is.null(component)) {
    cc <- sca(pred.ob[[var]])
  } else {
    cc <- sca(pred.ob[[component]][[var]])
  }
  
  return(cc)
}

# Plot components
ggplot()+
        gg(pred.comp$int_cont, aes(fill = mean)) + # change here according to layer
        pc("int_cont", "mean") +
        #gg(pts, size = .8, alpha = .4) +
        coord_equal()

# We select the model 4 as the best one, considering WAIC and exploratory analysis
sel.model <- 4



# Save individual components effects for the chosen model ----
# This was already done after running the models as a exploratory tool
# Now we run again and save the results

# Create a function to save results
dir <- paste("results", species, "effects", sep = "/")
dir_create(dir)

save.rast <- function(x, model){
        r <- as(pred.comp[[x]], "RasterStack")
        for (i in 1:nlayers(r)) {
          writeRaster(r[[i]], paste0(dir, "/", species, "_m", model, "_",
                                x, "_", names(r)[i], "_effect.tif"),
                      overwrite = T, format="GTiff")
        }
        return(paste(x, "saved."))
}

# Predict model components
pred.comp <- predict(m[[sel.model]], data = df,
                     formula = ~ list(
                       # Individual components
                       coldm = coldm,
                       sal = sal,
                       dist = dist,
                       ph = ph,
                       spatial = spatial),
                     n.samples = nsamp, seed = 2932)

# save
sapply(names(pred.comp), save.rast, model = sel.model)



# Predictions ----
# inlabru gets the environmental layer from the R environment
# thus, to predict to different scenarios we need to change the environmental
# layers object (i.e. 'env.e') to the new object for which we want the prediction

# We generate 2 predictions: 
# 1 with the spatial component and 1 without it (contrast)

dir <- paste("results", species, "predictions", sep = "/")
dir_create(dir)

for (i in 1:2) {
    
   # Get formula
    pred.f <- list(
      as.formula("~exp(coldm + sal + dist + spatial + ph + Intercept)"),
      as.formula("~exp(coldm + sal + dist + ph + Intercept)")
    )[[i]]
    
    # Save name
    snam <- paste0("_m", sel.model, "_", c("int", "cont")[i])
    
    cat("Running", snam, "\n")
    print(pred.f)
    
    # Predict to current scenario
    cat("Predicting for current scenario... \n")
    env.e <- as(env, "SpatialPixelsDataFrame")
    
    pred.cur <- predict(m[[sel.model]],
                        data = df,
                        formula = pred.f,
                        n.samples = nsamp, seed = 2932)
    ggplot() +
        gg(pred.cur, aes(fill = mean)) + pc(pred.ob = pred.cur) + coord_equal()
    
    # Predict to SSP 1
    cat("Predicting for SSP1 scenario... \n")
    env.e <- as(r12, "SpatialPixelsDataFrame")
    
    pred.r12 <- predict(m[[sel.model]],
                        data = df,
                        formula = pred.f,
                        n.samples = nsamp, seed = 2932)
    ggplot() +
        gg(pred.r12, aes(fill = mean)) + coord_equal()
    
    # Predict to SSP 2
    cat("Predicting for SSP2 scenario... \n")
    env.e <- as(r24, "SpatialPixelsDataFrame")
    
    pred.r24 <- predict(m[[sel.model]],
                        data = df,
                        formula = pred.f,
                        n.samples = nsamp, seed = 2932)
    ggplot() +
        gg(pred.r24, aes(fill = mean)) + pc(pred.ob = pred.r24) + coord_equal()
    
    # Predict to SSP 3
    cat("Predicting for SSP3 scenario... \n")
    env.e <- as(r37, "SpatialPixelsDataFrame")
    
    pred.r37 <- predict(m[[sel.model]],
                        data = df,
                        formula = pred.f,
                        n.samples = nsamp, seed = 2932)
    ggplot() +
        gg(pred.r37, aes(fill = mean)) + coord_equal()
    
    
    # Convert all to raster (easier handling)
    cat("Saving files... \n")
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
    
    # Raster was not being saved by layers for some reason
    # Thus we just use a small function to get the work done
    save.rast <- function(x, fnames){
      for (z in 1:nlayers(x)) {
        writeRaster(x[[z]], fnames[z], overwrite = T)
      }
      return(invisible(NULL))
    }
    
    save.rast(pred.cur, paste0(dir, "/", species, "_",
                               names(pred.cur), snam,
                               "_current.tif"))
    
    save.rast(pred.r12, paste0(dir, "/", species, "_",
                               names(pred.r12), snam,
                               "_ssp1.tif"))
    
    save.rast(pred.r24, paste0(dir, "/", species, "_",
                               names(pred.r24), snam,
                               "_ssp2.tif"))
    
    save.rast(pred.r37, paste0(dir, "/", species, "_",
                               names(pred.r37), snam,
                               "_ssp3.tif"))
    
    cat("Done! \n")
}



# Model cross-validation ----
# We cross validate based on the integration points
# We need the integrated lambda (intensity * weights)
# To get the "predicted value" on the presence points
# we look for the closest integration point for the presence point

# Note that this is not a perfect way of cross validation because each sample
# of points is a realization of the point process. Even though, it should give
# a glimpse of how models are performing. We also get the COR and the
# Boyce index (just for the full data).
# For an example of CV for point process (but using another modeling approach)
# see: https://doi.org/10.1016/B978-0-12-815226-3.00003-X
# https://github.com/thopitz/landslides-point-process/blob/master/Chapter.R

dir <- paste0("results/", species)

library(blockCV)
library(pROC)

spat.range <- spatialAutoRange(env, showPlots = F) # get ideal range for the blocks

# Divide in blocks
blocks <- spatialBlock(pts,
                       theRange = spat.range$range,
                       rasterLayer = env[[1]])
# Get the ID
blocks <- blocks$foldID

# Get distance to points of the full dataset
distpts <- pointDistance(ips, pts)

# Get the closest points
clpts <- apply(distpts, 2, function(x){which(x == min(x))})

# Get cross validate values
pa.pts <- ips
pa.pts$pa <- 0
pa.pts[clpts,"pa"] <- 1

cv.roc <- list()

for (i in 1:5) {
  
  cat("Block", i, "\n")
  
  cv.pts <- pts[blocks != i,]
  
  cv.m <- lgcp(forms[[sel.model]], cv.pts,
               ips = ips,
               options = list(control.inla = list(int.strategy = intest,
                                                  strategy = "simplified.laplace"),
                              control.mode=list(restart=T, theta=m[[sel.model]]$mode$theta),
                              inla.mode = "classic"))
  
  pred.ips <- predict(cv.m, data = ips,
                      formula = ~ weight * exp(coldm + sal + dist + spatial + ph + Intercept))
  
  pred.ips$cv.pred <- 1 - exp(-pred.ips$mean)
  
  cv.roc[[i]] <- roc(pa.pts$pa ~ pred.ips$cv.pred)
  
}

# Predict on all integration points
pred.ips.full <- predict(m[[sel.model]], data = ips,
                         formula = ~ weight * exp(coldm + sal + dist + spatial + ph + Intercept))


# Get ROC full and means of ROC
pa.pts$pred <- 1 - exp(-pred.ips.full$mean)
model.roc <- roc(pa.pts$pa~pa.pts$pred)

model.roc.cv <- list()
model.roc.cv$auc <- mean(unlist(lapply(cv.roc, function(x){as.numeric(x$auc)})))
model.roc.cv$auc.sd <- sd(unlist(lapply(cv.roc, function(x){as.numeric(x$auc)})))
model.roc.cv$spec <- apply(sapply(cv.roc, function(x){as.numeric(x$specificities)}),
                           1, mean)
model.roc.cv$sens <- apply(sapply(cv.roc, function(x){as.numeric(x$sensitivities)}),
                           1, mean)

# Save ROC plot
jpeg(paste0(dir, "/cv_results.jpg"), quality = 100)
plot(1 - model.roc$specificities, model.roc$sensitivities,
     xlim = c(0, 1), ylim = c(0, 1), col = "blue", type = "l",
     xlab = "1-Specificity", ylab = "Sensitivity", main = "ROC curve")
for (j in 1:5) {lines(1 - cv.roc[[j]]$specificities, cv.roc[[j]]$sensitivities, col = "gray")}
lines(1 - model.roc.cv$spec, model.roc.cv$sens)
legend("bottomright", legend = c("Fit", "CV"), col = c("blue", "black"),
       lty = c(1, 1))
dev.off()

# Save results
write.csv(data.frame(auc = model.roc$auc, auc_cv = model.roc.cv$auc,
                     auc_cv_sd = model.roc.cv$auc.sd),
          paste0(dir, "/", species, "_cv_results.csv"),
          row.names = F)


# Save metrics and summaries ----
dir <- paste0("results/", species)
write.csv(m.metrics, paste0(dir, "/", species, "_model_metrics.csv"),
          row.names = F)

# Extract summaries of the chosen model
model.summ <- m[[sel.model]]$summary.fixed

sigma.mar <- inla.tmarginal(function(x) exp(x),
                            m[[sel.model]]$marginals.hyperpar[[1]])
sigma.mar <- inla.zmarginal(sigma.mar)

range.mar <- inla.tmarginal(function(x) exp(x),
                            m[[sel.model]]$marginals.hyperpar[[2]])
range.mar <- inla.zmarginal(range.mar)

sst.sd.mar <- inla.zmarginal(m[[sel.model]]$marginals.hyperpar[[3]])

hyper <- rbind(cbind(as.data.frame(sigma.mar), mode = NA, kld = NA),
               cbind(as.data.frame(range.mar), mode = NA, kld = NA),
               cbind(as.data.frame(sst.sd.mar), mode = NA, kld = NA))
row.names(hyper) <- c("sigma", "range", "temp_sd")
hyper <- hyper[,-c(4, 6)]
names(hyper) <- names(model.summ)

model.summ <- rbind(model.summ, hyper)

# Save summaries
write.csv(model.summ, paste0(dir, "/", species, "_model_summary.csv"))

# Get plots of probability density
plots <- list()

# x label titles
xl <- c("Intercept", "Salinity", "Distance from coast", "pH")

for (i in 1:length(m[[sel.model]]$names.fixed)) {
  plots[[i]] <- plot(m[[sel.model]], m[[sel.model]]$names.fixed[i])+
    ylab(ifelse(i == 1, "Density", ""))+xlab(xl[i])+theme_classic()
}

pl <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]]
ggsave(paste0(dir, "/effects_density.jpg"), pl)

# Get plots of probability density for the SPDE parameters
sigp <- ggplot(data.frame(
  inla.smarginal(
    inla.tmarginal(function(x) exp(x),
                   m[[sel.model]]$marginals.hyperpar$`Theta1 for spatial`))),
  aes(x, y)) + geom_line() + xlab("Sigma") + ylab("Density") + theme_classic()

rangep <- ggplot(data.frame(
  inla.smarginal(
    inla.tmarginal(function(x) exp(x),
                   m[[sel.model]]$marginals.hyperpar$`Theta2 for spatial`))),
  aes(x, y)) + geom_line() + xlab("Range") + ylab("") + theme_classic()

sstp <- ggplot(data.frame(
  inla.smarginal(m[[sel.model]]$marginals.hyperpar[[3]])),
  aes(x, y)) + geom_line() + xlab("SD(SST)") + ylab("") + theme_classic()

sigp + rangep + sstp
ggsave(paste0(dir, "/hyperpar_density.jpg"), height = 3)




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