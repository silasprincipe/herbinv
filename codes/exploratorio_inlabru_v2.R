ext.f <- function(x, y, rast) {
        #print(length(x))
        # turn coordinates into SpatialPoints object:
        # with the appropriate coordinate reference system (CRS)
        spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = fm_sp_get_crs(rast))
        v <- extract(rast, spp)
        #v <- inla.group(v, n = 20)
        v
}

ext.f2 <- function(x, y, rast) {
        #print(length(x))
        # turn coordinates into SpatialPoints object:
        # with the appropriate coordinate reference system (CRS)
        spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = fm_sp_get_crs(rast))
        v <- extract(rast, spp)
        l <- length(v)
        v <- c(v, -3, 2.5)
        v <- inla.group(v, n = 20)
        v <- v[1:l]
        v
}


####
sp.env <- as(env.e, "SpatialPixelsDataFrame")

sp.env$sstg <- inla.group(sp.env$sst, n = 20)



pcprior <- list(theta = list(prior='pc.prec', param=c(1,0.01), initial = log(25)))

cmp.bm <- coordinates ~ 
        sst(sp.env, model = "linear") +
        #spatial(coordinates, model = b.model, mapper = bru_mapper(mesh)) +
        Intercept(1)

ipts <- ipoints(starea, mesh)


fit.2 <- lgcp(cmp.bm, pts,
              ips = ipts,
              domain = list(coordinates = mesh),
              options = list(control.inla = list(int.strategy = "eb"),
                             inla.mode = "classic",
                             bru_verbose = TRUE,
                             verbose = TRUE))

cmp.nl <- coordinates ~ 
        sstg(sp.env, model = "rw1", hyper = pcprior) +
        #spatial(coordinates, model = b.model, mapper = bru_mapper(mesh)) +
        Intercept(1)

fit.3 <- lgcp(cmp.nl, pts,
              ips = ipts,
              domain = list(coordinates = mesh),
              options = list(control.inla = list(int.strategy = "eb"),
                             inla.mode = "classic",
                             bru_verbose = TRUE,
                             verbose = TRUE))

plot.random(fit.3)



sp.env$sstc <- seq(-2.8, 1.8, length.out = length(sp.env$sst))


elev.pred <- predict(
        fit.3,
        data = sp.env,
        formula = ~ sstg_eval(sstc)
)

elev.pred <- elev.pred@data

ggplot(elev.pred) +
        geom_line(aes(sstc, mean)) +
        geom_ribbon(aes(sstc,
                        ymin = q0.025,
                        ymax = q0.975),
                    alpha = 0.2) +
        geom_ribbon(aes(sstc,
                        ymin = mean - 1 * sd,
                        ymax = mean + 1 * sd),
                    alpha = 0.2)+
        abline()







cmp.sp <- coordinates ~ spatial(coordinates, model = spde) +
        Intercept(1)

cmp.bm_t <- coordinates ~ spatial(coordinates, model = b.model, mapper = bru_mapper(mesh)) +
        sst(ext.f(x, y, env.e$sst), model = "rw2", hyper = pcprior) +
        Intercept(1)

cmp.sp_t <- coordinates ~ spatial(coordinates, model = spde) +
        sst(ext.f(x, y, env.e$sst), model = "rw2", hyper = pcprior) +
        Intercept(1)

cmp.sst <- coordinates ~
        sst(ext.f(x, y, env.e$sst), model = "rw2", hyper = pcprior, constr = T) +
        Intercept(1)

cmp.sst_lin <- coordinates ~
        sst(ext.f(x, y, env.e$sst), model = "linear") +
        Intercept(1)

cmp.sst2 <- coordinates ~
        sst(ext.f2(x, y, env.e$sst), model = "rw2", hyper = pcprior, constr = T) +
        spatial(coordinates, model = b.model, mapper = bru_mapper(mesh)) +
        Intercept(1)


fit.1 <- lgcp(cmp.bm, pts, samplers = starea,
              ips = ips,
                domain = list(coordinates = mesh),
                options = list(control.inla = list(int.strategy = "eb"),
                               inla.mode = "classic",
                               bru_verbose = TRUE,
                               verbose = TRUE))

fit.2 <- lgcp(cmp.sp, pts, samplers = starea,
              ips = ips,
              domain = list(coordinates = mesh),
              options = list(control.inla = list(int.strategy = "eb"),
                             inla.mode = "classic",
                             bru_verbose = TRUE,
                             verbose = TRUE))

fit.3 <- lgcp(cmp.bm_t, pts, samplers = starea,
              ips = ips,
              domain = list(coordinates = mesh),
              options = list(control.inla = list(int.strategy = "eb"),
                             inla.mode = "classic",
                             bru_verbose = TRUE,
                             verbose = TRUE))

fit.4 <- lgcp(cmp.sp_t, pts, samplers = starea,
              ips = ips,
              domain = list(coordinates = mesh),
              options = list(control.inla = list(int.strategy = "eb"),
                             inla.mode = "classic",
                             bru_verbose = TRUE,
                             verbose = TRUE))

fit.5 <- lgcp(cmp.sst, pts, samplers = starea,
              ips = ips,
              domain = list(coordinates = mesh),
              options = list(control.inla = list(int.strategy = "eb"),
                             inla.mode = "classic",
                             bru_verbose = TRUE,
                             verbose = TRUE))

fit.6 <- lgcp(cmp.sst_lin, pts, samplers = starea,
              ips = ips,
              domain = list(coordinates = mesh),
              options = list(control.inla = list(int.strategy = "eb"),
                             inla.mode = "classic",
                             bru_verbose = TRUE,
                             verbose = TRUE))

fit.8 <- lgcp(cmp.sst2, pts, samplers = starea,
              ips = ips,
              domain = list(coordinates = mesh),
              options = list(control.inla = list(int.strategy = "eb"),
                             inla.mode = "classic",
                             bru_verbose = TRUE,
                             verbose = TRUE))

### sst spde
knots <- seq(-3, 2, length = 25)
d1mesh <- inla.mesh.1d(knots, interval = c(-3, 2), degree = 2)

d1spde <- inla.spde2.pcmatern(d1mesh,
                              prior.range = c(1, NA),
                              prior.sigma = c(1, 0.1))

cmp.sst <- coordinates ~
        sst(env.e, model = d1spde) +
        Intercept(1)

fit.7 <- lgcp(cmp.sst, pts,
              ips = ips,
              domain = list(coordinates = mesh),
              options = list(control.inla = list(int.strategy = "eb"),
                             inla.mode = "classic",
                             bru_verbose = TRUE,
                             verbose = TRUE))


cmp.sst_spde <- coordinates ~
        sst(sp.env, model = d1spde, mapper = bru_mapper(d1mesh, indexed = T)) +
        spatial(coordinates, model = b.model, mapper = bru_mapper(mesh)) +
        Intercept(1)

fit.7_s <- lgcp(cmp.sst_spde, pts,
              ips = ipts,
              domain = list(coordinates = mesh),
              options = list(control.inla = list(int.strategy = "eb"),
                             inla.mode = "classic",
                             bru_verbose = TRUE,
                             verbose = TRUE))


cmp.sst_spde2 <- coordinates ~
        sst(sp.env, model = d1spde, mapper = bru_mapper(d1mesh, indexed = T)) +
        spatial(coordinates, model = spde) +
        Intercept(1)

fit.7_s2 <- lgcp(cmp.sst_spde2, pts,
                ips = ipts,
                domain = list(coordinates = mesh),
                options = list(control.inla = list(int.strategy = "eb"),
                               inla.mode = "classic",
                               bru_verbose = TRUE,
                               verbose = TRUE))



df <- pixels(mesh, mask = starea, nx = 400, ny = 400)
int1 <- predict(fit.7, data = df, ~ exp(sst + Intercept))
int2 <- predict(fit.7_s, data = df, ~ exp(sst + spatial + Intercept))
int3 <- predict(fit.7_s2, data = df, ~ exp(sst + spatial + Intercept))

ggplot()+gg(int1)+coord_equal()
ggplot()+gg(int2)+coord_equal()
ggplot()+gg(int3)+coord_equal()

# Non pc matern
d1spde <- inla.spde2.matern(d1mesh,
                            # prior.range = c(2, NA),
                            # prior.sigma = c(0.1, 0.01)
                            theta.prior.prec = 1e-4,
                            constr = T)
### funcionando

sp.env <- env.e

sp.env$sstc <- seq(-2.6, 1.8, length.out = length(sp.env$sst))

elev.pred <- predict(
        m[[7]],
        data = sp.env,
        formula = ~ sst_eval(sstc)
)

elev.pred <- elev.pred@data

ggplot(elev.pred) +
        geom_line(aes(sstc, mean)) +
        geom_ribbon(aes(sstc,
                        ymin = q0.025,
                        ymax = q0.975),
                    alpha = 0.2) +
        geom_ribbon(aes(sstc,
                        ymin = mean - 1 * sd,
                        ymax = mean + 1 * sd),
                    alpha = 0.2)+
        abline()


df <- pixels(mesh, mask = starea, nx = 400, ny = 400)
int1 <- predict(m[[7]], data = df, ~ exp(sst + sal + Intercept))

int.pred <- int1["mean"]

int1$mean <- int1$mean * 1230.775

int2 <- int1

int2$mean <- 1-exp(-int1$mean)

ggplot() +
        gg(int2) +
        #gg(starea, alpha = 0, lwd = 1) +
        #gg(pts, color = "DarkGreen") +
        coord_equal()


23.90894 * 30.02841

#
exposure <- c(ipts$weight, rep(0, length(pts)))
pa <- c(rep(0, length(ipts)), rep(1, length(pts)))

sst <- c(extract(env.e$sst, ipts), extract(env.e$sst, pts))

fit.gam <- gam(pa ~ s(sst, k = 10),
               family = poisson())


pred <- data.frame(pa = NA, sst =  seq(-2.8, 1.9, length = 100))

pred2.gam <- predict(fit.gam, newdata = pred, type = "response")


plot(pred2.gam)


# Predict abundance >>> should be something close to 311
Lambda.exp <- predict(
        fit.8,
        ipoints(starea, mesh),
        ~ sum(weight * exp(spatial + sst + Intercept))
);Lambda.exp

sst.pred <- predict(
        fit.8,
        data = data.frame(sst2 = seq(-3, 2, length.out = 100)),
        formula = ~ sst_eval(sst2)
)

ggplot(sst.pred) +
        geom_line(aes(sst2, mean)) +
        geom_ribbon(aes(sst2,
                        ymin = q0.025,
                        ymax = q0.975),
                    alpha = 0.2) +
        geom_vline(xintercept = max(sp.data$sst))+
        geom_ribbon(aes(sst2,
                        ymin = mean - 1 * sd,
                        ymax = mean + 1 * sd),
                    alpha = 0.2)

sst.pred <- predict(
        fit.5,
        data = data.frame(sst = seq(-3, 2, length.out = 100)),
        formula = ~ exp(sst + Intercept)
)

ggplot(sst.pred) +
        geom_line(aes(sst2, mean)) +
        geom_ribbon(aes(sst2,
                        ymin = q0.025,
                        ymax = q0.975),
                    alpha = 0.2) +
        geom_vline(xintercept = max(sp.data$sst))+
        geom_ribbon(aes(sst2,
                        ymin = mean - 1 * sd,
                        ymax = mean + 1 * sd),
                    alpha = 0.2)



### MESH 1D FUNCIONA
cmp.sst2 <- coordinates ~
        #sst(ext.f2(x, y, env.e$sst), model = "rw2", hyper = pcprior, constr = T) +
        spatial(coordinates, model = b.model, mapper = bru_mapper(mesh)) +
        sst_spde(ext.f(x, y, env.e$sst), model = spde1d, mapper = bru_mapper(mesh1d, indexed = T)) +
        Intercept(1)

fit.8 <- lgcp(cmp.sst2, pts, samplers = starea,
              ips = ips,
              domain = list(coordinates = mesh),
              options = list(control.inla = list(int.strategy = "eb"),
                             inla.mode = "classic",
                             bru_verbose = TRUE,
                             verbose = TRUE))


#### IID
cmp.sst2 <- coordinates ~
        sst(ext.f2(x, y, env.e$sst), model = "rw2", hyper = pcprior, constr = T) +
        spatial(coordinates, model = b.model, mapper = bru_mapper(mesh)) +
        site(1:3542, model = "iid", n = 3542)+
        #sst_spde(ext.f(x, y, env.e$sst), model = spde1d, mapper = bru_mapper(mesh1d, indexed = T)) +
        Intercept(1)

fit.8 <- lgcp(cmp.sst2, pts, samplers = starea,
              ips = ips,
              domain = list(coordinates = mesh),
              options = list(control.inla = list(int.strategy = "eb"),
                             inla.mode = "classic",
                             bru_verbose = TRUE,
                             verbose = TRUE))
