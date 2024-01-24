#### Modelling of coral reef herbivorous invertebrates ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2022 ##

## Barrier Model under the Log-Gaussian Cox Process framework ##

## Get response curves and variable importance
# This step was done in 2024, prior to publication


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
source("functions/response_curves.R")

# Settings
set.seed(2932) # Replicability of sampling
spatt <- theme_classic()+theme(axis.line = element_blank()) # Theme for better ploting
# Integration strategy (change here for "eb" for a fast approximation)
intest <- "auto"
# Number of samples to calculate posteriors (decrease to a fast approx.)
nsamp <- 1000



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

# Warmest Month
knots.wm <- seq(-4, 2.5, length = 30)
d1mesh.wm <- inla.mesh.1d(knots.wm, degree = 2,
                          boundary = "free")

d1spde.wm <- inla.spde2.pcmatern(d1mesh.wm,
                                 prior.range = c(diff(c(-4, 2.5)), NA),
                                 prior.sigma = c(0.5, 0.5),
                                 constr = T)



##### Write formulas ----

forms <- list()

base.f <- coordinates ~ 
  Intercept(1, mean.linear = log(length(pts)/gArea(starea)), prec.linear = 1) + 
  sal(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
  dist(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
  spatial(coordinates, model = b.model, mapper = bru_mapper(mesh))

forms[[1]] <- update(base.f, ~ . + sst(env.e, model = d1spde.st))
forms[[2]] <- update(base.f, ~ . + coldm(env.e, model = d1spde.st))
forms[[3]] <- update(base.f, ~ . + coldm(env.e, model = d1spde.wm))

forms[[1]] <- update(forms[[1]], ~ . +
                       ph(env.e, model = "linear", mean.linear = 0, prec.linear = 1))
forms[[2]] <- update(forms[[2]], ~ . +
                       ph(env.e, model = "linear", mean.linear = 0, prec.linear = 1))
forms[[3]] <- update(forms[[3]], ~ . +
                       ph(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
                       chl(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
                       pho(env.e, model = "linear", mean.linear = 0, prec.linear = 1))

vars.resp <- list(c("sal", "dist", "sst", "ph"),
                  c("sal", "dist", "coldm", "ph"),
                  c("sal", "dist", "coldm", "ph", "chl", "pho"))

# Run models to get response curves ----
species.list <- c("lyva", "eclu", "trve")
m <- list()
respcurves <- list()
respcurves.exp <- list()
var.importance <- list()

for (i in 1:3) {
  
  species <- species.list[i]
  
  cat("Runing model for species", species, "\n")
  print(forms[[i]])
  
  pts <- read.csv(paste0("data/", species, "/", species, "_filt.csv"))
  proj <- "+proj=laea +lat_0=0 +lon_0=-70 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs"
  pts <- SpatialPoints(pts, proj4string = CRS(proj))
  
  m[[i]] <- lgcp(forms[[i]], pts,
                 ips = ips,
                 options = list(control.inla = list(int.strategy = intest,
                                                    strategy = "simplified.laplace"),
                                inla.mode = "classic"
                                # For verbosity, uncoment those lines:
                                # , bru_verbose = TRUE,
                                # verbose = TRUE
                 ))
  
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
  print(lambda.exp)
  cat("\n ====================== \n")
  
  cat("==== Getting response curves ====\n")
  
  respcurves[[i]] <- get.resp.curves(m[[i]], vars.resp[[i]], mode = NULL, samp = 1000)
  respcurves.exp[[i]] <- get.resp.curves(m[[i]], vars.resp[[i]], mode = "exp", samp = 1000)
  
  cat("==== Getting variable importance ====\n")
  old.env.e <- env.e
  iter <- 10
  vars <- vars.resp[[i]]
  fit <- m[[i]]
  
  start.pred <- predict(fit, rbind(pts, ips),
                        as.formula(paste0(
                          "~ exp(",
                          paste(if (is.null(names(
                            fit$summary.random))) {
                            fit$names.fixed
                          } else{
                            c(fit$names.fixed,
                              names(fit$summary.random))
                          }
                          , collapse = "+"), ")"
                        )))
  
  start.pred <- start.pred$mean
  
  results <- data.frame(variables = vars,
                        importance = NA,
                        importance_sd = NA)
  
  original.values <- env.e@data
  
  for (k in 1:length(vars)) {
    cat("Running variable", vars[k], "\n")
    
    cor.val <- rep(NA, iter)
    
    cli::cli_progress_bar(total = iter)
    for (it in 1:iter) {
      
      env.e[[vars[k]]] <-  env.e[[vars[k]]][sample(1:nrow(env.e))]
      
      new.pred <- predict(fit, rbind(pts, ips),
                          as.formula(paste0(
                            "~ exp(",
                            paste(if (is.null(names(
                              fit$summary.random))) {
                              fit$names.fixed
                            } else{
                              c(fit$names.fixed,
                                names(fit$summary.random)#[names(fit$summary.random) != "spatial"]
                                )
                            }
                            , collapse = "+"), ")"
                          )))
      
      new.pred <- new.pred$mean
      
      cor.val[it] <- 1-cor(start.pred, new.pred)
      
      env.e[[vars[k]]] <- original.values[,vars[k]]
      
      cli::cli_progress_update()
      
    }
    cli::cli_progress_done()
    
    results[results$variables == vars[k], "importance"] <- mean(cor.val, na.rm = T)
    results[results$variables == vars[k], "importance_sd"] <- sd(cor.val, na.rm = T)
  }
  
  var.importance[[i]] <- results
  
  env.e <- old.env.e
  
}

# Save importance of variables ----
lapply(1:3, function(x){
  sp <- species.list[x]
  write.csv(var.importance[[x]], paste0("results/", sp, "/", sp, "_varimportance.csv"), row.names = F)
  return(invisible(NULL))
})

# Save response curves ----
for (x in 1:3) {
  sp <- species.list[x]
  p.lin <- respcurves[[x]]
  p.exp <- respcurves.exp[[x]]
  plot(p.lin)
  ggsave(filename = paste0("results/", sp, "/", sp, "_respcurveslin.png"),
         width = 12, height = 8)
  plot(p.exp)
  ggsave(filename = paste0("results/", sp, "/", sp, "_respcurvesexp.png"),
         width = 12, height = 8)
}

#### END