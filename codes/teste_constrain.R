forms[[3]] <-
        coordinates ~ Intercept(1, mean.linear = log(length(pts)/gArea(starea)), prec.linear = 1) +
        sal(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        sst(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        dist(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        spatial(coordinates, model = b.model,
                mapper = bru_mapper(mesh), extraconstr = list(A = A, e = e))



forms[[4]] <-
        coordinates ~ Intercept(1, mean.linear = log(length(pts)/gArea(starea)), prec.linear = 1) +
        sal(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        sst(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        dist(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        spatial(coordinates, model = b.model, mapper = bru_mapper(mesh))





forms[[1]] <-
        coordinates ~ Intercept(1, mean.linear = log(length(pts)/gArea(starea)), prec.linear = 1) +
        sal(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        sst(env.e, model = d1spde, mapper = bru_mapper(d1mesh, indexed = T)) +
        dist(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        ph(env.e, model = "linear", mean.linear = 0, prec.linear = 1) 

forms[[2]] <-
        coordinates ~ Intercept(1, mean.linear = log(length(pts)/gArea(starea)), prec.linear = 1) +
        sal(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        sst(env.e, model = d1spde, mapper = bru_mapper(d1mesh, indexed = T)) +
        dist(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        chl(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        pho(env.e, model = "linear", mean.linear = 0, prec.linear = 1)

forms[[3]] <-
        coordinates ~ Intercept(1, mean.linear = log(length(pts)/gArea(starea)), prec.linear = 1) +
        sal(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        sst(env.e, model = d1spde, mapper = bru_mapper(d1mesh, indexed = T)) +
        dist(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        chl(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        pho(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        ph(env.e, model = "linear", mean.linear = 0, prec.linear = 1) 
        
forms[[4]] <- update(forms[[1]], ~. + spatial(coordinates, model = b.model, mapper = bru_mapper(mesh)))
forms[[5]] <- update(forms[[2]], ~. + spatial(coordinates, model = b.model, mapper = bru_mapper(mesh)))
forms[[6]] <- update(forms[[3]], ~. + spatial(coordinates, model = b.model, mapper = bru_mapper(mesh)))




forms[[7]] <-
        coordinates ~ Intercept(1, mean.linear = log(length(pts)/gArea(starea)), prec.linear = 1) +
        sal(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        sst(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        dist(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        chl(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        pho(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        nit(env.e, model = "linear", mean.linear = 0, prec.linear = 1)

forms[[8]] <-
        coordinates ~ Intercept(1, mean.linear = log(length(pts)/gArea(starea)), prec.linear = 1) +
        sal(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        sst(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        dist(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        chl(env.e, model = "linear", mean.linear = 0, prec.linear = 1)

forms[[8]] <-
        coordinates ~ Intercept(1, mean.linear = log(length(pts)/gArea(starea)), prec.linear = 1) +
        sal(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        sst(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        dist(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        chl(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        pho(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        sil(env.e, model = "linear", mean.linear = 0, prec.linear = 1)

forms[[9]] <-
        coordinates ~ Intercept(1, mean.linear = log(length(pts)/gArea(starea)), prec.linear = 1) +
        sal(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        sst(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        dist(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        chl(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        nit(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        sil(env.e, model = "linear", mean.linear = 0, prec.linear = 1)

forms[[10]] <-
        coordinates ~ Intercept(1, mean.linear = log(length(pts)/gArea(starea)), prec.linear = 1) +
        sal(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        sst(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        dist(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        chl(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        sil(env.e, model = "linear", mean.linear = 0, prec.linear = 1)


        
# Matérn models ----
# Get range and sigma
size <- diff(bbox(starea)[1,])
range0 <- as.numeric(round(size/5))
sigma <- 0.5

# Stationary model for comparison (if wanted)
b.model <- inla.barrier.pcmatern(mesh,
                                 barrier.triangles = barrier.triangles,
                                 prior.range = c(range0, 0.5),
                                 prior.sigma = c(sigma, 0.01))

spde <- inla.spde2.pcmatern(mesh,
                            prior.range = c(range0, 0.5),
                            prior.sigma = c(sigma, 0.01))

env.tot <- getd(env, ip = mesh$loc[,1:2])
nodes.d <- extract(env.tot, mesh$loc[,1:2])
#nodes.d[is.na(nodes.d[,8]),7:9] <- apply(nodes.d[,7:9], 2, mean, na.rm = T)

A <- matrix(1, ncol = mesh$n, nrow = 1)
e <- matrix(0, ncol = 1)

A[1,] <- nodes.d[,"sst"]


sstfun <- function(x, y){
        p <- SpatialPoints(cbind(x, y), proj4string = inlabru::fm_sp_get_crs(env.e))
        vals <- over(p, env.e[,"sst"])
        gvals <- inla.group(c(vals[,1], (max(vals)+1)), n = 20)
        return(gvals[1:(length(gvals)-1)])
}

forms[[1]] <-
        coordinates ~ Intercept(1, mean.linear = log(length(pts)/gArea(starea)), prec.linear = 1) +
        sal(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        sst(sstfun(x, y), model = "rw2", hyper = list(theta = list(prior='pc.prec', param=c(1, 0.01), initial = log(25)))) +
        dist(env.e, model = "linear", mean.linear = 0, prec.linear = 1) +
        ph(env.e, model = "linear", mean.linear = 0, prec.linear = 1) 

forms[[1]] <- update(forms[[1]], ~. + spatial(coordinates, model = b.model, mapper = bru_mapper(mesh)))

# precision of  0.001 = SD of 31
# 1/sqrt(0.001)

# SD to precision >>> 1/sd^2

# Run models in loop
for (i in 1) {
        
        cat("Runing model", i, "\n")
        print(forms[[i]])
        
        m[[i]] <- lgcp(forms[[i]], pts,
                       ips = ips,
                      # domain = list(coordinates = mesh),
                       options = list(control.inla = list(int.strategy = "eb"),
                                      inla.mode = "classic"
                                      # For verbosity, uncoment those lines:
                                      , bru_verbose = TRUE,
                                      verbose = TRUE
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


sca <- function(...){
        scale_fill_gradientn(
                colours = rev(c("#A84C00", "#D97D27", "#F5BD44", "#FFD561", "#FFF291",
                                "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695")),
                limits = range(...))
}


df <- as(env$ph, "SpatialPixelsDataFrame")
df <- df[,-1]

pred.comp <- lapply(1:6, function(x){
        predict(m[[x]], data = df,
                formula = ~ list(
                        # Change those according to mode and
                        # combinations
                        sst = sst,
                        sal = sal,
                        dist = dist,
                        #spatial = spatial,
                        int = as.formula(paste0(
                                "~exp(",
                                paste(if (is.null(names(
                                        m[[x]]$summary.random))) {
                                        m[[x]]$names.fixed
                                } else{
                                        c(m[[x]]$names.fixed,
                                          names(m[[x]]$summary.random))
                                }
                                , collapse = "+"), ")"
                        )))
                )
})

pred.comp[[6]] <- predict(m[[6]], data = df,
                          formula = ~ list(
                                  sst = sst,
                                  dist = dist,
                                  sal = sal,
                                  int = exp(sst + sal + dist + chl + pho + ph + spatial + Intercept)
                          ))

pred.comp[[2]] <- predict(m[[8]], data = df,
                          formula = ~ list(
                                  chl = chl,
                                  sal = sal,
                                  int = exp(sst + sal + dist + chl + Intercept)
                          ))


p5 <- ggplot()+
        gg(pred.comp[[5]]$int)+ # change here according to layer
        sca(pred.comp[[5]]$int$mean)+
        coord_equal()

ggplot()+
        gg(pred.comp[[6]]$int)+ # change here according to layer
        sca(pred.comp[[6]]$int$mean)+
        coord_equal()

library(patchwork)
p1+p2+p3+p4+p5


2-2/((qnorm(0.975, 0, 1)-qnorm(0.5, 0, 1)))




####
#the default priors specify a standard deviation of:
(SD<-(1/sqrt(0.25) * qnorm(0.975,0,1)))

prior <- data.frame(x=seq(-3*SD,3*SD,len=1000))
prior$density <- dnorm(prior$x,0,SD)
post <- data.frame(m[[1]]$marginals.fixed[[2]])
head(post)


ggplot(prior, aes(y=density, x=x)) + geom_line(aes(color='Prior')) +
        geom_line(data=post, aes(y=y, x=x, color='Posterior')) +
        scale_color_discrete('')+
        theme_classic()+theme(legend.position=c(1,1), legend.justification=c(1,1))
