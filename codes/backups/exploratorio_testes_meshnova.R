range0 <- as.numeric(round(size/5))
range0 = 2500
sigma <- 0.1

# Stationary model for comparison (if wanted)
spde <- inla.spde2.pcmatern(mesh, constr = T,
                            prior.range = c(range0, 0.5),
                            prior.sigma = c(sigma, 0.01))

# Barrier model
b.model <- inla.barrier.pcmatern(mesh,
                                 barrier.triangles = barrier.triangles,
                                 prior.range = c(range0, 0.5),
                                 prior.sigma = c(sigma, 0.01))


forms <- list()

forms[[1]] <- coordinates ~ Intercept(1) + spatial(main = coordinates, model = spde)
forms[[1]] <- coordinates ~ Intercept(1) + spatial(main = coordinates, model = b.model, mapper = bru_mapper(mesh))

forms[[2]] <- update(forms[[1]], ~ . + sal(env.e, model = "linear"))

forms[[3]] <- update(forms[[2]], ~ . + sal(env.e, model = "linear") + sst(env.e, model = d1spde,
                                                                      mapper = bru_mapper(d1mesh, indexed = T)))

forms[[4]] <- update(forms[[2]], ~ . + sal(env.e, model = "linear") + sst(env.e, model = "linear"))

forms[[5]] <- coordinates ~ Intercept(1) + sal(env.e, model = "linear") + sst(env.e, model = "linear")

m <- list()

# Creates a data.frame to hold WAIC, DIC and CPO (summary) results
m.metrics <- data.frame(waic = rep(NA, length(forms)), dic = NA, cpo = NA)
m.lambda <- list()

# Run models in loop
for (i in 1:length(forms)) {
        
        cat("Runing model", i, "\n")
        print(forms[[i]])
        
        m[[i]] <- lgcp(forms[[i]], pts.c,
                       ips = ips,
                       domain = list(coordinates = mesh),
                       options = list(control.inla = list(int.strategy = "eb",
                                                          strategy="simplified.laplace"),
                                      control.compute(cpo = T),
                                      # control.fixed(mean = list(Intercept = log(length(pts)/gArea(starea))),
                                      #               prec = list(Intercept = 1)),
                                      #control.fixed(prec = list(spatial = 1)),
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


# Predict each component effect spatially
df <- pixels(mesh, mask = starea, nx = 150, ny = 150)
p1 <- predict(m[[1]], data = df,
                     formula = ~ list(
                             spatial = spatial,
                             int = exp(spatial + Intercept)))
p2 <- predict(m[[2]], data = df,
              formula = ~ list(
                      spatial = spatial,
                      sal = sal,
                      int = exp(sal + spatial + Intercept)))
p3 <- predict(m[[3]], data = df,
              formula = ~ list(
                      spatial = spatial,
                      sal = sal,
                      sst = sst,
                      int = exp(sal + sst + spatial + Intercept)))
p4 <- predict(m[[4]], data = df,
              formula = ~ list(
                      spatial = spatial,
                      sal = sal,
                      sst = sst,
                      int = exp(sal + sst + spatial + Intercept)))
p5 <- predict(m[[5]], data = df,
              formula = ~ list(
                      sal = sal,
                      sst = sst,
                      int = exp(sal + sst + Intercept)))

(pp1 <- ggplot()+gg(p1$int)+ scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = median(p1$int$mean))+ coord_equal())
(pp2 <- ggplot()+gg(p2$int)+ scale_fill_viridis()+ coord_equal()+#coord_equal(ylim = c(0,4800)))+
        colsc(p2$int$mean))
(pp3 <- ggplot()+gg(p3$teste)+ scale_fill_viridis()+ coord_equal()+colsc(p3$teste$mean))
(pp4 <- ggplot()+gg(p4$int)+ scale_fill_viridis()+ coord_equal())
(pp5 <- ggplot()+gg(p5$int)+ scale_fill_viridis()+ coord_equal())

p3$teste <- p3$int
p3$teste$mean <- 1-exp(-(p3$int$mean))

library(patchwork)

(pp1 + pp2) / (pp3 + pp4)

fit <- m[[1]]

int.plot <- plot(fit, "Intercept")
spde.range <- spde.posterior(fit, "spatial", what = "range")
spde.logvar <- spde.posterior(fit, "spatial", what = "variance")
range.plot <- plot(spde.range)
var.plot <- plot(spde.logvar)

multiplot(range.plot, var.plot, int.plot)

corplot <- plot(spde.posterior(fit, "spatial", what = "matern.correlation"))
covplot <- plot(spde.posterior(fit, "spatial", what = "matern.covariance"))
multiplot(covplot, corplot)

plot.spatial(m[[2]], mesh, addb = F,# pts = pts
             )











D <- dist(coordinates(pts))
par(mfrow = c(1,2), mar = c(5,5,2,2), cex.lab = 1.5)
hist(D, 
     freq = TRUE,
     main = "", 
     xlab = "Dist. between sites (km)",
     ylab = "Frequency")

plot(x = sort(D), 
     y = (1:length(D))/length(D), 
     type = "l",
     xlab = "Dist. sites (km)",
     ylab = "Cumulative proportion")
