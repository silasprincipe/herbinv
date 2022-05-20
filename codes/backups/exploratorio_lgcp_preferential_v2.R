pts.c <- read.csv("data/lyva/lyva_final.csv")[,2:3]
colnames(pts.c) <- c("x", "y")
coordinates(pts.c) <- ~x+y
crs(pts.c) <- CRS(projargs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
pts.c <- spTransform(pts.c, CRS(proj))

test <- env[[1]]
test[] <- NA
for (z in 1:length(pts.c)) {
        test[cellFromXY(test, pts.c@coords[z,])] <-
                ifelse(is.na(test[cellFromXY(test, pts.c@coords[z,])]), 1,
                       (test[cellFromXY(test, pts.c@coords[z,])] + 1))
}

test <- mask(test, env[[1]])

test <- rasterToPoints(test)

colnames(test) <- c("x", "y", "count")
test <- data.frame(test)
coordinates(test) <- ~x+y
crs(test) <- crs(env)

pts.c <- test


#gaus.prior <- list(prior = 'gaussian', param = c(0, 2))
cmp = ~ spde(coordinates, model = spde) +
        spdeCopy(coordinates, copy = "spde", fixed = FALSE) +
        sal(env.e, model = "linear")+
        sst(env.e, model = d1spde,
            mapper = bru_mapper(d1mesh, indexed = T))+
        sstCopy(env.e, copy = "sst", fixed = F)+
        #sst(env.e, model = "linear")+
        # salCopy(copy = "sal", fixed = FALSE)+
        # sstCopy(copy = "sst", fixed = FALSE)+
        lgcpIntercept(1) +
        Intercept(1)

lik1 <- like(data =  pts.c,
             family = "poisson",
             formula = count ~ spde + sst + sal + Intercept)

lik2 <- like(data =  pts.c,
             family = "cp",
             ips = ips,
             #domain = list(coordinates = mesh),
             formula = coordinates ~ spdeCopy + sstCopy + sal + lgcpIntercept)




#gaus.prior <- list(prior = 'gaussian', param = c(0, 2))
cmp = ~ spde(coordinates, model = spde) +
        spdeCopy(coordinates, copy = "spde", fixed = FALSE) +
        sal(env.e, model = "linear")+
        sst(env.e, model = d1spde,
            mapper = bru_mapper(d1mesh, indexed = T))+
        #sst(env.e, model = "linear")+
        # salCopy(copy = "sal", fixed = FALSE)+
        # sstCopy(copy = "sst", fixed = FALSE)+
        lgcpIntercept(1) +
        Intercept(1)

lik1 <- like(data =  pts.c,
             family = "poisson",
             formula = count ~ spde + Intercept)

lik2 <- like(data =  pts.c,
             family = "cp",
             ips = ips,
             #domain = list(coordinates = mesh),
             formula = coordinates ~ spdeCopy + sst + sal + lgcpIntercept)

# Even though both likelihoods are log-linear, INLA needs some
# assistance in finding the modes; bru_max_iter = 3 seems sufficient here, but
# allowing more iterations is safer; this iteration stops at 4.
fit <- bru(components = cmp,
           lik1,
           lik2,
           options = list(bru_max_iter = 2,
                          bru_verbose = 4,
                          control.inla = list(int.strategy = "eb")))

predict(fit, ips, ~sum(weight * exp(spdeCopy + sst + sal + lgcpIntercept)))


pxl <- pixels(mesh, mask = starea)
fish.intensity <- predict(fit, pxl, ~ list(
        spde = spde,
        spde_lg = spdeCopy,
        effects = exp(lgcpIntercept + sst + sal),
        count = exp(spde + Intercept),
        point = exp(spdeCopy + lgcpIntercept + sst + sal)
))

fish.intensity_b <- predict(fitb, pxl, ~ list(
        spde = spde,
        #spde_lg = spdeCopy,
        effects = exp(lgcpIntercept + sst + sal),
        count = exp(spde + Intercept),
        point = exp(spde_b + lgcpIntercept + sst + sal)
))

colsc <- function(...) {
        scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"RdYlBu")),
                             limits = range(...))
}


p1 <- ggplot() +
        gg(fish.intensity$point) +
        #gg(pts.c, color = "white", size = .8, alpha = .5)+
        #gg(pts.c, size = .5, color = "white")+
        # scale_fill_stepsn(#breaks = seq(0.1,0.9,length.out = 9),
        #                   n.breaks = 9,
        #                   # limits = c(min(fish.intensity$point$mean),
        #                   #            max(fish.intensity$point$mean)),
        #                   colors = RColorBrewer::
        #                           brewer.pal(9, "GnBu"),
        #                   na.value = NA)+
        colsc(fish.intensity$point$mean)+
       # scale_fill_viridis()+
        gg(pts.c, size = .8, color = "black", alpha = .3)+
        theme_classic()+
        ggtitle("Preferencial")+
        coord_equal()

p2 <- ggplot() +
        gg(fish.intensity$effects) +
        #gg(pts.c, color = "white", size = .8, alpha = .5)+
        #gg(pts.c, size = .5, color = "white")+
        # scale_fill_stepsn(#breaks = seq(0.1,0.9,length.out = 9),
        #         n.breaks = 9,
        #         # limits = c(min(fish.intensity$point$mean),
        #         #            max(fish.intensity$point$mean)),
        #         colors = RColorBrewer::
        #                 brewer.pal(9, "GnBu"),
        #         na.value = NA)+
        colsc(fish.intensity$effects$mean)+
        # scale_fill_viridis()+
        gg(pts.c, size = .8, color = "black", alpha = .3)+
        theme_classic()+
        ggtitle("Pref. - Apenas efeitos ambientais")+
        coord_equal()



p3 <- ggplot() +
        gg(fish.intensity_b$point) +
        #gg(pts.c, color = "white", size = .8, alpha = .5)+
        #gg(pts.c, size = .5, color = "white")+
        # scale_fill_stepsn(#breaks = seq(0.1,0.9,length.out = 9),
        #         n.breaks = 9,
        #         # limits = c(min(fish.intensity$point$mean),
        #         #            max(fish.intensity$point$mean)),
        #         colors = RColorBrewer::
        #                 brewer.pal(9, "GnBu"),
        #         na.value = NA)+
        colsc(fish.intensity_b$point$mean)+
        gg(pts.c, size = .8, color = "black", alpha = .3)+
        # scale_fill_viridis()+
        theme_classic()+
        ggtitle("Normal")+
        coord_equal()

p4 <- ggplot() +
        gg(fish.intensity_b$effects) +
        #gg(pts.c, color = "white", size = .8, alpha = .5)+
        #gg(pts.c, size = .5, color = "white")+
        # scale_fill_stepsn(#breaks = seq(0.1,0.9,length.out = 9),
        #         n.breaks = 9,
        #         # limits = c(min(fish.intensity$point$mean),
        #         #            max(fish.intensity$point$mean)),
        #         colors = RColorBrewer::
        #                 brewer.pal(9, "GnBu"),
        #         na.value = NA)+
        colsc(fish.intensity_b$effects$mean)+
        gg(pts.c, size = .8, color = "black", alpha = .3)+
        # scale_fill_viridis()+
        theme_classic()+
        ggtitle("Normal - Apenas efeitos ambientais")+
        coord_equal()

library(patchwork)
(p1+p3)
(p2+p4)


cmp.lg = coordinates ~ spde(coordinates, model = spde) +
        sal(env.e, model = "linear")+
        sst(env.e, model = d1spde,
            mapper = bru_mapper(d1mesh, indexed = T))+
        #spdeCopy(coordinates, copy = "spde", fixed = FALSE) +
        #lgcpIntercept(1) +
        Intercept(1)

fit.lg <- lgcp(cmp.lg, data = pts.c, 
              ips = ips,
           options = list(bru_max_iter = 2,
                          bru_verbose = 4,
                          control.inla = list(int.strategy = "eb")))

fish.intensity.lg <- predict(fit.lg, pxl, ~ exp(spde + Intercept))
ggplot() +
        gg(fish.intensity.lg) +
        gg(pts.c, size = .5, color = "white")+
        #gg(pts, size = .5, color = "orange")+
        scale_fill_viridis()+
        coord_equal()



predict(fit, ipoints(domain = shrimp$mesh), ~ sum(weight * exp(spde + lgcpIntercept)))




lik1 <- like(data =  pts.c,
             family = "poisson",
             formula = count ~ spde + Intercept)

lik2 <- like(data =  pts.c,
             family = "cp",
             ips = ips,
             #domain = list(coordinates = mesh),
             formula = coordinates ~ spdeCopy + sst + sal + lgcpIntercept)

# Even though both likelihoods are log-linear, INLA needs some
# assistance in finding the modes; bru_max_iter = 3 seems sufficient here, but
# allowing more iterations is safer; this iteration stops at 4.
fit2 <- bru(components = cmp,
           lik1,
           lik2,
           options = list(bru_max_iter = 2,
                          bru_verbose = 4,
                          control.inla = list(int.strategy = "eb")))



lik1 <- like(data =  pts.c,
             family = "poisson",
             formula = count ~ spde + sst + sal  + Intercept)

lik2 <- like(data =  pts.c,
             family = "cp",
             ips = ips,
             #domain = list(coordinates = mesh),
             formula = coordinates ~ spdeCopy + lgcpIntercept)

# Even though both likelihoods are log-linear, INLA needs some
# assistance in finding the modes; bru_max_iter = 3 seems sufficient here, but
# allowing more iterations is safer; this iteration stops at 4.
fit3 <- bru(components = cmp,
            lik1,
            lik2,
            options = list(bru_max_iter = 5,
                           bru_verbose = 4,
                           control.inla = list(int.strategy = "eb")))




lik1 <- like(data =  pts.c,
             family = "poisson",
             formula = count ~ spde + Intercept)

lik2 <- like(data =  pts.c,
             family = "cp",
             ips = ips,
             #domain = list(coordinates = mesh),
             formula = coordinates ~ spdeCopy + lgcpIntercept)

# Even though both likelihoods are log-linear, INLA needs some
# assistance in finding the modes; bru_max_iter = 3 seems sufficient here, but
# allowing more iterations is safer; this iteration stops at 4.
fit4 <- bru(components = cmp,
            lik1,
            lik2,
            options = list(bru_max_iter = 10,
                           bru_verbose = 4,
                           control.inla = list(int.strategy = "eb")))



lik1 <- like(data =  pts.c,
             family = "poisson",
             formula = count ~ spde + Intercept)

lik2 <- like(data =  pts.c,
             family = "cp",
             ips = ips,
             #domain = list(coordinates = mesh),
             formula = coordinates ~ spdeCopy + lgcpIntercept)

# Even though both likelihoods are log-linear, INLA needs some
# assistance in finding the modes; bru_max_iter = 3 seems sufficient here, but
# allowing more iterations is safer; this iteration stops at 4.
fit5 <- bru(components = cmp,
            lik2,
            options = list(bru_max_iter = 5,
                           bru_verbose = 4,
                           control.inla = list(int.strategy = "eb")))




deltaIC(fit, fit2, fit3, fit4, fit5, criterion = "WAIC")









form_lgcp <- coordinates ~ b0.lgcp + spat1

form_col <- count ~ b0.count + spat2


lik_lgcp <- like(family = 'cp', data = pts.c,
                 ips = ips,
                 domain = list(coordinates = mesh),
                 formula = form_lgcp)

lik_col <- like(family = 'gaussian', data = pts.c,
                ips = ips,
                domain = list(coordinates = mesh),
                formula = form_col)

fit <- bru(components = cmp, lik_lgcp, lik_col,
           options = list(control.inla=list(int.strategy = 'eb')))

summary(fit)

pfit <- predict(fit, data = df,
                     formula = ~ list(
                             spat1 = spat1,
                             spat2 = spat2,
                             lgcp = exp(b0.lgcp + spat1),
                             col = exp(b0.count + spat2),
                             total = exp(b0.lgcp + spat1 + b0.count)
                     ))

predict(fit, ips,
        ~sum(weight * exp(b0.lgcp + spat1)))

ggplot()+
        gg(pfit$lgcp)+ # change here according to layer
        scale_fill_viridis()+
        #gg(pts, color = 'white')+
        coord_equal()

p1




###### SHRIMP EXAMPLE
data(shrimp)

matern <- inla.spde2.pcmatern(shrimp$mesh, alpha = 2,
                              prior.sigma = c(0.1, 0.01),
                              prior.range = c(1, 0.01))

plot.bmodel(shrimp$hauls@coords[10,1:2], mesh = shrimp$mesh, spde = matern,
            areapol = test, spmode = T, range = 1)

cmp = ~ spde(coordinates, model = matern) +
        spdeCopy(coordinates, copy = "spde", fixed = FALSE) +
        lgcpIntercept(1) +
        Intercept(1)

# The catch data has some non-integer value, so Poisson is inappropriate
# Using the pseudo-model "xpoisson" that allows non-integer observations instead
lik1 <- like(data =  shrimp$hauls,
             family = "xpoisson",
             formula = catch ~ spde + Intercept)

lik2 <- like(data =  shrimp$hauls,
             family = "cp",
             domain = list(coordinates = shrimp$mesh),
             formula = coordinates ~ spdeCopy + lgcpIntercept)

# Even though both likelihoods are log-linear, INLA needs some
# assistance in finding the modes; bru_max_iter = 3 seems sufficient here, but
# allowing more iterations is safer; this iteration stops at 4.
fit <- bru(components = cmp,
           lik1,
           lik2,
           options = list(bru_max_iter = 20,
                          bru_verbose = 4))

predict(fit, ipoints(domain = shrimp$mesh), ~ sum(weight * exp(spde + lgcpIntercept)))

pxl <- pixels(shrimp$mesh)
fish.intensity <- predict(fit, pxl, ~ exp(spde + Intercept))
ggplot() +
        gg(fish.intensity) +
        gg(shrimp$hauls, size = .5, color = "white")+
        coord_equal(xlim = c(-0.6, 0), ylim = c(37.5, 38))

predict(fit, formula = ~ exp(Intercept + spde), integrate = "coordinates")

