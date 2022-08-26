
abs <- sampleRandom(env[[c("sal", "coldm", "dist", "ph")]], 25000, xy = F)
pre <- extract(env[[c("sal", "coldm", "dist", "ph")]], pts)

#dat <- rbind(cbind(abs, y = 0), cbind(coordinates(pts), pre, y = 1))
dat <- data.frame(rbind(cbind(abs, y = 0), cbind(pre, y = 1)))


boosted.ipp <- gbm(y ~ sal + coldm + I(coldm^2) + dist + ph,
                   distribution = "bernoulli",
                   data=dat, weights=1E3^(1-y))
summary(boosted.ipp)
plot(boosted.ipp)

brt.pred <- predict(env, boosted.ipp)
plot(exp(brt.pred))
plot(exp(brt.pred), xlim = c(-4000,0), ylim = c(1000,4500))

gam.ipp <- gam(y ~ sal + s(coldm) + dist + ph, family=binomial, data=dat,
               weights=1E3^(1-y))

summary(gam.ipp)

gam.pred <- predict(env, gam.ipp)
plot(exp(gam.pred))
plot(exp(gam.pred), xlim = c(-4000,0), ylim = c(1000,4500))

st <- gArea(starea)
p.wt  <- rep(1.e-6, length(dat$y))
p.wt[dat$y == 0] <- st/sum(dat$y == 0)
dwglm.ipp <- glm(y/p.wt ~ sal + coldm + I(coldm^2) + dist + ph,
                 data = dat,
                 family = poisson(), weights = p.wt)

dwglm.pred <- predict(env, dwglm.ipp)
plot(exp(dwglm.pred))
plot(exp(dwglm.pred), xlim = c(-4000,0), ylim = c(1000,4500))

up.wt <- (10^6)^(1 - dat$y)
iwlr.ipp <- glm(y ~ sal + coldm + I(coldm^2) + dist + ph,
                data = dat,
                family = binomial(), weights = up.wt)

iwlr.pred <- predict(env, iwlr.ipp)
plot(exp(iwlr.pred))
plot(exp(iwlr.pred), xlim = c(-4000,0), ylim = c(1000,4500))

###
library(spatstat)

abs <- sampleRandom(env[[c("sal", "coldm", "dist", "ph")]], 25000, xy = T)
dat <- rbind(cbind(abs, dps = 0), cbind(coordinates(pts), pre, dps = 1))

env.im <- as(env, "SpatialPixelsDataFrame")
env.im <- as.data.frame(env.im)
env.im <- as.im(env.im)
env.im <- lapply(1:nlayers(env),
                 function(x){
                   as.im.SpatialGridDataFrame(as(env[[x]], "SpatialGridDataFrame"))
                 })
names(env.im) <- names(env)
env.im[[9]] <-  as.im.SpatialGridDataFrame(as((env$coldm^2), "SpatialGridDataFrame"))

owin.ppp <- as.owin.SpatialPolygons(starea)

owin.ppp <- owin(poly = starea)
pts.ppp <- ppp(coordinates(pts)[,1], coordinates(pts)[,2],
               window = owin.ppp)

plot(env.im)

Q <- quadscheme(pts.ppp)

spat.mod <- ppm(Q,  trend = ~sal+coldm+coldm_sq+dist+ph, covariates = env.im)

spat.pred <- predict(spat.mod)
plot(spat.pred)

spat.pred.rast <- as(spat.pred, "RasterLayer")
plot(spat.pred.rast)

lgcp.bru.pred <- raster("results/eclu/predictions/eclu_mean_m4_cont_current.tif")
lgcp.bru.pred.int <- raster("results/eclu/predictions/eclu_mean_m4_int_current.tif")

####
par(mfrow = c(3,3))
plot(exp(brt.pred), main = "brt")
plot(exp(gam.pred), main = "gam")
