#### Modelling of coral reef herbivorous invertebrates ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2022 ##

## Tests with common SDM methods adjusted for point process ##

# Load packages ----
library(terra)
library(SDMtune)
library(blockCV)
library(sp)

# Create folder to hold results ----
fs::dir_create("supplement/maxent")

# Load environmental data ----
env <- rast(list.files("data/env/ready_layers", pattern = "_cur", full.names = T))
r12 <- rast(list.files("data/env/ready_layers", pattern = "_r12", full.names = T))
r24 <- rast(list.files("data/env/ready_layers", pattern = "_r24", full.names = T))
r37 <- rast(list.files("data/env/ready_layers", pattern = "_r37", full.names = T))

names(env) <- c("chl", "coldm", "ph", "pho", "sal", "sst", "warmm")
names(r12) <- names(r24) <- names(r37) <- names(env)

env <- c(env, rast("data/env/ready_layers/distcoast.tif"))
names(env)[8] <- "dist"
r12$dist <- r24$dist <- r37$dist <- env[[8]]

pred.l <- list(env, r12, r24, r37)
names(pred.l) <- c("current", "ssp1", "ssp2", "ssp3")

# Define species ----
species <- c("lyva", "eclu", "trve")

# Define variables for each species according to results ----
forms <- list()
forms[[1]] <- c("sal", "dist", "sst", "ph")
forms[[2]] <- c("sal", "dist", "coldm", "ph")
forms[[3]] <- c("sal", "dist", "coldm", "ph", "chl", "pho")

# See number of background/quadrature points necessary ----
n.quad <- c(5000, 10000, 20000, 30000, 37000)
area.size <- expanse(env[[1]], unit = "km")
area.size <- area.size[,2]
proj <- "+proj=laea +lat_0=0 +lon_0=-70 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs"
par(mfrow = c(1,3))

for (z in 1:3) {
  lik.list <- list()
  
  sp <- species[z]
  
  # Load species data
  pts <- read.csv(paste0("data/", sp, "/", sp, "_filt.csv"))[,1:2]
  pts <- SpatialPoints(data.frame(x = pts[,1], y = pts[,2]), 
                       proj4string = CRS(proj))
  
  # Get environmental layers for the species
  env.sel <- subset(env, forms[[z]])
  
  for (j in 1:50) {
    
    lik <- rep(NA, length(n.quad))
    
    for (i in 1:length(n.quad)) {
      quad.test <- data.frame(spatSample(env.sel, n.quad[i], values = F,
                                         xy = T, na.rm = T))
      
      pts.xy <- coordinates(pts)
      colnames(pts.xy) <- c("x", "y")
      dat <- rbind(cbind(pts.xy, sp = 1, terra::extract(env.sel, pts.xy)),
                   cbind(quad.test, sp = 0, terra::extract(env.sel, quad.test, ID = F)))
      dat <- dat[!is.na(dat$sal),]
      
      p.wt <- rep(1.e-8, nrow(dat))
      p.wt[dat$sp == 0] <- (area.size/1e6)/n.quad[i]
      
      resp <- dat$sp/p.wt
      
      dat$resp <- resp
      
      dwpr <- glm(as.formula(paste("resp ~", paste(forms[[z]], collapse = "+"))),
                  family = poisson(), weights = p.wt, data = dat)
      
      mu <- dwpr$fitted
      
      lik[i] <- sum(p.wt*(resp*log(mu) - mu))
    }
    
    lik.list[[j]] <- lik
    
    cat(j, "|")
  }
  
  plot(y = lik.list[[1]], x = n.quad, type = "o", log = "x",
       ylim = range(unlist(lik.list)), main = sp)
  lapply(2:length(lik.list), function(x){
    lines(y = lik.list[[x]], x = n.quad, type = "o", col = sample(colors(), 1))
  })
}

# We will use 30000, almost all available points

# Run in sequence for all species ----
for (i in 1:3) {
        
        sp <- species[i]
        
        # Load species data
        pts <- read.csv(paste0("data/", sp, "/", sp, "_filt.csv"))[,1:2]
        pts <- SpatialPoints(data.frame(x = pts[,1], y = pts[,2]), 
                             proj4string = CRS(proj))
        
        # Get environmental layers for the species
        env.sel <- subset(env, forms[[i]])
        
        # Get background points
        bck.pts <- spatSample(env.sel, size = 30000, na.rm = T, xy = T)
        bck.pts <- SpatialPointsDataFrame(coords = bck.pts[,1:2],
                                          data = bck.pts[,3:length(bck.pts)],
                                          proj4string = CRS(proj))
        
        # Create the SDM tune object
        pts.coords <- data.frame(coordinates(pts))
        bck.coords <- data.frame(coordinates(bck.pts))
        names(pts.coords) <- names(bck.coords) <- c("x", "y")
        
        swd.obj <- prepareSWD(species = sp,
                              env = env.sel,
                              p = pts.coords,
                              a = bck.coords)
        
        # Get blocks for cross validation
        blocks.size <- cv_spatial_autocor(
          x = bck.pts[sample(1:nrow(bck.pts), 5000),],
          column = names(bck.pts)
        )
        
        full.pts <- SpatialPointsDataFrame(
          coords = swd.obj@coords,
          data = data.frame(presence = swd.obj@pa),
          proj4string = crs(pts)
        )
        
        blocks <- cv_spatial(full.pts,
                             column = "presence",
                             r = env.sel,
                             k = 5,
                             size = blocks.size$range)
        
        # Fit a first model
        max.m <- train("Maxnet", swd.obj, fc = "lq")
        
        # Define the hyperparameters to test
        hyper <- list(reg = seq(0.5, 5, 0.5),
                      fc = c("l", "lq"))
        
        # Test all the possible combinations
        tuning <- gridSearch(max.m, hypers = hyper,
                             metric = "aicc",
                             env = env)
        
        # See results
        head(tuning@results[order(tuning@results$delta_AICc), ])
        
        # Get the best configuration
        best.t <- tuning@results[order(tuning@results$delta_AICc), 1:2][1,]
        
        
        # Run final model
        # We run one with the blocks for cross-validation and one
        # with the full data.
        max.m <- train("Maxnet", swd.obj, fc = best.t$fc, folds = blocks,
                       reg = best.t$reg)
        
        auc(max.m) # Print AUC
        
        max.m.f <- train("Maxnet", swd.obj, fc = best.t$fc,
                         reg = best.t$reg)
        
        for (pl in 1:length(pred.l)) {
          env.pred <- subset(pred.l[[pl]], forms[[i]])
          pred.exp <- predict(max.m.f, env.pred, type = "exponential")
          pred.clog <- predict(max.m.f, env.pred, type = "cloglog")
          writeRaster(pred.exp, paste0("supplement/maxent/", sp, "_max_exp_",
                                       names(pred.l)[pl], ".tif"),
                      overwrite = TRUE)
          writeRaster(pred.clog, paste0("supplement/maxent/", sp, "_max_clog_",
                                       names(pred.l)[pl], ".tif"),
                      overwrite = TRUE)
        }
}

# Explore results ----
# Define species
sp <- "lyva"

# Define type
type <- "clog"

# Load raster
current <- rast(paste0("supplement/maxent/", sp, "_max_" , type, "_current.tif"))
ssp1 <- rast(paste0("supplement/maxent/", sp, "_max_" , type, "_ssp1.tif"))
ssp2 <- rast(paste0("supplement/maxent/", sp, "_max_" , type, "_ssp2.tif"))
ssp3 <- rast(paste0("supplement/maxent/", sp, "_max_" , type, "_ssp3.tif"))

if (type == "exp") {
  max.curr <- global(current, max, na.rm = T)[,1]
  ssp1[ssp1 > max.curr] <- max.curr
  ssp2[ssp2 > max.curr] <- max.curr
  ssp3[ssp3 > max.curr] <- max.curr
}

par(mfrow = c(1,1))
plot(current)
par(mfrow = c(1,3))
plot(ssp1);plot(ssp2);plot(ssp3)

# Deltas
ssp1.delta <- ssp1 - current
ssp2.delta <- ssp2 - current
ssp3.delta <- ssp3 - current

par(mfrow = c(1,3))
plot(ssp1.delta);plot(ssp2.delta);plot(ssp3.delta)

# Compare direct
current.max <- rast(paste0("supplement/maxent/lyva_max_exp_current.tif"))
current.lgcp <- rast(paste0("results/lyva/predictions/lyva_mean_m4_cont_current.tif"))

ssp3.max <- rast(paste0("supplement/maxent/lyva_max_exp_ssp3.tif"))
ssp3.lgcp <- rast(paste0("results/lyva/predictions/lyva_mean_m4_cont_ssp3.tif"))

par(mfrow = c(1,2))
plot(current.max)
plot(current.lgcp)

current.max <- rast(paste0("supplement/maxent/eclu_max_exp_current.tif"))
current.lgcp <- rast(paste0("results/eclu/predictions/eclu_mean_m4_cont_current.tif"))

ssp3.max <- rast(paste0("supplement/maxent/eclu_max_exp_ssp3.tif"))
ssp3.lgcp <- rast(paste0("results/eclu/predictions/eclu_mean_m4_cont_ssp3.tif"))

par(mfrow = c(1,2))
plot(current.max)
plot(current.lgcp)

plot(ssp3.max)
plot(ssp3.lgcp)
