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
# Settings
set.seed(2932) # Replicability of sampling
spatt <- theme_classic() # Theme for better ploting


# Load environmental and species data ----
source("functions/Varload.r")
lays <- c("salinitymean", "tempmean", "ph", "chlomean")

env <- var.load(layers = lays)

r12 <- var.load(folder = "proj_layers/ssp126/", layers = lays)
r24 <- var.load(folder = "proj_layers/ssp245/", layers = lays)
r37 <- var.load(folder = "proj_layers/ssp370/", layers = lays)

# Load study area shapefile
starea <- shapefile("gis/starea.shp")

# Mask environmental layers
env <- mask(env, starea)
r12 <- mask(r12, starea)
r24 <- mask(r24, starea)
r37 <- mask(r37, starea)

# Change names for easier handling
names(env) <- c("ph", "chl", "sal", "sst")
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

# Matérn models
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
                                 prior.range = c(range0, 0.01),
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
env.e$sstg <- inla.group(env.e$sst, n = 20)

# Create models formulas ----

# Creates a PC prior for random walks
pcprior <- list(theta = list(prior='pc.prec', param=c(1,0.01), initial = log(25)))

# Establish model formulas
# We will test 7 different models
forms <- list()

forms[[1]] <- coordinates ~ 
        sal(env.e, model = "linear")+
        sst(env.e, model = "linear") +
        Intercept(1)

forms[[2]] <- coordinates ~ 
        sal(env.e, model = "linear")+
        sst(env.e, model = "linear") +
        spatial(coordinates, model = spde)+
        Intercept(1)

forms[[3]] <- coordinates ~ 
        sal(env.e, model = "linear")+
        sst(env.e, model = "linear") +
        spatial(coordinates, model = b.model, mapper = bru_mapper(mesh)) +
        Intercept(1)

forms[[4]] <- coordinates ~ 
        sal(env.e, model = "linear")+
        sstg(env.e, model = "rw2", hyper = pcprior) +
        Intercept(1)

forms[[5]] <- coordinates ~ 
        sal(env.e, model = "linear")+
        sstg(env.e, model = "rw2", hyper = pcprior) +
        spatial(coordinates, model = spde)+
        Intercept(1)

forms[[6]] <- coordinates ~ 
        sal(env.e, model = "linear")+
        sstg(env.e, model = "rw2", hyper = pcprior) +
        spatial(coordinates, model = b.model, mapper = bru_mapper(mesh)) +
        Intercept(1)


knots <- seq(-3, 2, length = 25)
d1mesh <- inla.mesh.1d(knots, interval = c(-3, 2), degree = 2)
d1spde <- inla.spde2.pcmatern(d1mesh,
                              prior.range = c(2, NA),
                              prior.sigma = c(1, 0.01))


forms[[7]] <- coordinates ~ 
        sal(env.e, model = "linear")+
        sst(env.e, model = d1spde, mapper = bru_mapper(d1mesh, indexed = T)) +
        Intercept(1)

forms[[8]] <- coordinates ~ 
        sal(env.e, model = "linear")+
        sst(env.e, model = d1spde, mapper = bru_mapper(d1mesh, indexed = T)) +
        spatial(coordinates, model = spde)+
        Intercept(1)

forms[[9]] <- coordinates ~ 
        sal(env.e, model = "linear")+
        sst(env.e, model = d1spde, mapper = bru_mapper(d1mesh, indexed = T)) +
        spatial(coordinates, model = b.model, mapper = bru_mapper(mesh)) +
        Intercept(1)


# Run models ----
# Initiate a list to hold models
m <- list()

# Creates a data.frame to hold WAIC, DIC and CPO (summary) results
m.metrics <- data.frame(waic = rep(NA, length(forms)), dic = NA, cpo = NA)
m.lambda <- list()

# # Get mean occurrences per area
# avg.global <- sum(pa)/gArea(starea)
# 
# ips <- ipoints(starea, mesh)

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
source("functions/plot_random_effects.R")
plot.random(m[[9]])
plot.post(m[[1]])

# Plot spatial effect
plot.spatial(m[[1]], mesh)

# Compare models
deltaIC(m[[1]], m[[2]], m[[3]], m[[4]], m[[5]], m[[6]], m[[7]])

# We proceed cross-validation with model 1
sel.model <- 4




### Cross-validation ----

# To cross-validate the models we consider that it's possible to convert
# the intensity to a proability value 
# (see Lombardo et al. 2017; https://arxiv.org/abs/1708.03156).
# In that case, values are converted to probability as 1-exp(-intensity)

# We use a spatial block cross validation, which leaves entire sections of
# the environment out in each run of the model.
# As a first step, we divide the area in blocks, using package blockCV

spat.range <- spatialAutoRange(env, showPlots = F) # get ideal range for the blocks

# Get coordinates and PA
dat <- data.frame(cbind(rbind(mesh$loc[,1:2], pts@coords), pa = pa))
coordinates(dat) <- ~ x + y

# Divide in blocks
blocks <- spatialBlock(dat, species = "pa",
                       theRange = spat.range$range)
# Get the ID
blocks <- blocks$foldID

# For which model(s) run the CV?
selmodels <- c(which(m.metrics$dic %in% m.metrics$dic[order(m.metrics$dic)][1:2]))

cvm.preds <- list()

for (j in 1:length(selmodels)) {
        
        cat("Runing model", selmodels[j], "CV. \n")
        
        # Creates a vector to hold predicted intesities
        pred.intensity <- c()
        
        for (i in 1:5) {
                cat("Block", i, "\n")
                
                # Set vector of response for prediction
                pa.cv <- pa
                pa.cv[which(blocks == i)] <- NA
                
                # Set CV stack
                stk.cv <- inla.stack(
                        data = list(y = pa.cv, e = ev), 
                        A = list(A.f, 1), 
                        effects = list(list(spatial = 1:nv),
                                       list(b0 = rep(1, nrow(sp.data)), 
                                            sp.data)),
                        tag = 'est')
                
                # Run model for the selected model
                cvm <- inla(forms[[selmodels[j]]],
                                 family = 'poisson', data = inla.stack.data(stk.cv), 
                                 control.predictor=list(A=inla.stack.A(stk.cv),
                                                        compute=TRUE, link=1), 
                                 control.fixed = list(mean = list(b0 = log(avg.global)),
                                                      prec = list(b0 = 1, sal = 2)),
                                 control.mode=list(theta=m[[selmodels[j]]]$mode$theta,
                                                   restart=FALSE),
                                 E = inla.stack.data(stk.cv)$e,
                                 control.inla=list(int.strategy = "eb"),
                                 num.threads = 3)
                
                # Get predicted intensity
                idx <- inla.stack.index(stk.cv, "est")$data
                
                pred.intensity[which(blocks == i)] <-
                        cvm$summary.fitted.values$mean[idx[which(blocks == i)]]
        }
        
        rm(cvm)
        
        # Fit values for the full model
        idx <- inla.stack.index(stk.i, "est")$data
        full.intensity <- m[[selmodels[j]]]$summary.fitted.values$mean[idx]
        
        #### Calculate ROC curves and AUC values
        
        # Convert to probability
        probability <- 1-exp(-full.intensity) 
        probability.cv <- 1-exp(-pred.intensity)
        
        # Calculate ROC and AUC
        m.roc <- roc(pa ~ probability)
        m.roc.cv <- roc(pa ~ probability.cv)
        m.auc <- as.numeric(m.roc$auc)
        m.auc.cv <- as.numeric(m.roc.cv$auc)
        
        cvm.preds[[j]] <- list(
                "fint" = full.intensity,
                "pint" = pred.intensity,
                "mroc" = m.roc,
                "mrocv" = m.roc.cv,
                "fauc" = m.auc,
                "auccv" = m.auc.cv
        )
}

# See results
cvm.preds[[1]]$fauc;cvm.preds[[2]]$fauc
cvm.preds[[1]]$auccv;cvm.preds[[2]]$auccv

# Plot ROC
par(mfrow = c(1,2))
for (z in 1:2) {
        plot(1 - cvm.preds[[z]]$mroc$specificities, cvm.preds[[z]]$mroc$sensitivities,
             xlim = c(0, 1), ylim = c(0, 1), col = "blue", type = "l",
             xlab = "1-Specificity", ylab = "Sensitivity", main = "ROC curve")
        lines(1 - cvm.preds[[z]]$mrocv$specificities, cvm.preds[[z]]$mrocv$sensitivities)
        legend("bottomright", legend = c("Fit", "CV"), col = c("blue", "black"),
               lty = c(1, 1))
}



### Prediction ----
# Creates a function to extract data for predictions
extract.pred.data <- function(rast){
        pred.data <- data.frame(rasterToPoints(aggregate(rast, 5)))
        pred.data <- na.omit(pred.data)
        coordinates(pred.data) <- ~ x + y
        pred.data
        
}

# Extract data for prediction
pdat <- lapply(list(env, r12, r24, r37), extract.pred.data)
names(pdat) <- c("current", "ssp1", "ssp2", "ssp3")

# Creates a new matrix and stack
A.pred <- lapply(pdat, function(x){inla.spde.make.A(mesh, loc = x)})

stk.list <- list()
for (i in 1:length(pdat)) {
        
        stk.list[[i]] <- inla.stack(
                data = list(y = NA), 
                A = list(A.pred[[i]], 1), 
                effects = list(list(spatial = 1:nv),
                               list(b0 = 1, data.frame(pdat[[i]]@data))),
                tag = names(pdat)[i])
        
}

# Create full stack
stk.full <- inla.stack(stk.i, stk.list[[1]], stk.list[[2]], stk.list[[3]], stk.list[[4]])

pred.list <- list()

for (j in 1:length(forms)) {
        # Run prediction model
        model.pred <- inla(forms[[j]],
                           family = 'poisson', data = inla.stack.data(stk.full), 
                           control.predictor=list(A=inla.stack.A(stk.full),
                                                  compute=TRUE, link = 1), 
                           E = inla.stack.data(stk.full)$e,
                           control.mode=list(theta=m[[j]]$mode$theta,
                                             restart=FALSE),
                           control.fixed = list(mean = list(b0 = log(avg.global)),
                                                prec = list(b0 = 1, sal = 2)),
                           control.compute=list(return.marginals.predictor=TRUE),
                           control.inla=list(int.strategy = "eb"),
                           num.threads = 3, inla.mode = "classic")
        
        # Creates a function to extract different types of prediction
        extract.pred <- function(name, type = "mean"){
                # Get index of predictions
                index <- inla.stack.index(stk.full, name)$data
                
                # Get predictions
                fit.pred <- data.frame(
                        "m" = model.pred$summary.fitted.values[[type]][index])
                
                # Convert to raster 
                pred.raster <- rasterFromXYZ(cbind(pdat[[1]]@coords, fit.pred$m))
                names(pred.raster) <- name
                pred.raster
        }
        
        # Extract predictions for each scenario and convert to raster
        pred.mean <- stack(lapply(names(pdat), extract.pred, type = "mean"))
        
        #plot(pred.mean)
        
        pred.list[[j]] <- pred.mean
}

par(mfrow = c(1,2))
plot(pred.list[[sel.model]][[1]]);plot(pred.list[[sel.model]][[4]])
#lapply(pred.list, function(x){plot(x[[1]]);plot(x[[2]])})


### Convert to exceedance probability ----
# Get exceddance probability

minv <- extract(pred.list[[sel.model]][[1]], pts)
prob <- min(minv)


# Creates a function to extract different types of prediction
extract.excprob <- function(name, prob){
        # Get index of predictions
        index <- inla.stack.index(stk.full, name)$data
        
        excprob <- sapply(model.pred$marginals.fitted.values,
                          function(marg){
                                  1-inla.pmarginal(q = prob, marginal = marg)})
        
        # Get predictions
        fit.pred <- data.frame(
                "m" = model.pred$summary.fitted.values[[type]][index])
        
        # Convert to raster 
        pred.raster <- rasterFromXYZ(cbind(pdat[[1]]@coords, fit.pred$m))
        names(pred.raster) <- name
        pred.raster
}

# Extract predictions for each scenario and convert to raster
pred.mean <- stack(lapply(names(pdat), extract.excprob, prob = prob))





