library(maxnet)
library(raster)
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
species <- "lyva"

pts <- read.csv(paste0("data/", species, "/", species, "_filt.csv"))
proj <- "+proj=laea +lat_0=0 +lon_0=-70 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs"
pts <- SpatialPoints(pts, proj4string = CRS(proj))



abs <- sampleRandom(env[[c("sal", "sst", "dist")]], 25000, xy = T)
pre <- extract(env[[c("sal", "sst", "dist")]], pts)

dat <- rbind(cbind(abs, dps = 0), cbind(coordinates(pts), pre, dps = 1))

max <- maxnet(dat[,6], data.frame(dat[,3:5]), addsamplestobackground = F, classes = "lq")
plot(max)

pred.max <- predict(env, max, type = "exponential")

pred.max.fut <- predict(r37, max, type = "exponential")

plot(pred.max)

pred.max.gg <- as(pred.max, "SpatialPixelsDataFrame")
pred.max.fut.gg <- as(pred.max, "SpatialPixelsDataFrame")

a <- ggplot()+gg(pred.max.fut.gg)+scale_fill_viridis()+coord_equal()
b <- ggplot()+gg(fut.pred.original)+scale_fill_viridis()+coord_equal()

a+b

pred.max.clog <- predict(env, max, type = "cloglog")
pred.max.clog2 <- 1-exp(-exp(max$entropy+predict(env, max, type = "link")))


r <- raster(
  "C:/pesquisa/herbinv/results/lyva/predictions/lyva_mean_m4_int_current.tif"
)

r2 <- raster(
  "C:/pesquisa/herbinv/results/lyva/predictions/lyva_mean_m4_cont_current.tif"
)


r

rclog <- 1-exp(-r*85469.8)
rclog


rclog2 <- 1-exp(-r2*85469.8)

plot(rclog);plot(rclog2)



library(inlabru)
library(INLA)


# Load environmental and species data ----
# These are prepared in the code lgcp_prepare_data.R
list2env(readRDS('data/lgcp_data.rds'),globalenv())

cmp <- coordinates ~ 
  Intercept(1) + 
  sal(env.e, model = "linear") +
  dist(env.e, model = "linear") +
  sst(env.e, model = "linear") +
  sst_sq(env.e, model = "linear") +
  spatial(coordinates, model = spde)

ips <- ipoints(starea, mesh)


size <- diff(bbox(starea)[1,])
range0 <- as.numeric(round(size/10))
sigma <- 0.1

# Stationary model for comparison (if wanted)
spde <- inla.spde2.pcmatern(mesh,
                            prior.range = c(range0, 0.1),
                            prior.sigma = c(sigma, 0.01))


fit <- lgcp(cmp,
            pts,
            #formula = coordinates ~ Intercept + sal + dist + spatial + poly(coldm,degree = 2),
            ips = ips,
            options = list(control.inla = list(int.strategy = "eb",
                                               strategy = "simplified.laplace"),
                           inla.mode = "classic"
                           # For verbosity, uncoment those lines:
                           # , bru_verbose = TRUE,
                           # verbose = TRUE
            ))


df <- as(env$ph, "SpatialPixelsDataFrame")
df <- df[,-1]

pred.comp <- predict(fit, data = df,
                     formula = ~ exp(sst + sst_sq + sal + dist + spatial + Intercept))

pred.comp <- as(pred.comp, "RasterStack")

pred.comp$clog <- pred.comp$mean
pred.comp$clog <- 1-exp(-pred.comp$clog*85469.8)

plot(pred.comp$mean)
plot(pred.comp$clog)

plot(pred.max.clog)
plot(pred.comp$clog)
