#### Modelling of coral reef herbivorous invertebrates ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2022 ##

## Tests with common SDM methods adjusted for point process ##
env <- stack(list.files("data/env/ready_layers", pattern = "_cur", full.names = T))
r12 <- stack(list.files("data/env/ready_layers", pattern = "_r12", full.names = T))
r24 <- stack(list.files("data/env/ready_layers", pattern = "_r24", full.names = T))
r37 <- stack(list.files("data/env/ready_layers", pattern = "_r37", full.names = T))

names(env) <- c("chl", "coldm", "ph", "sal", "sst", "warmm")
names(r12) <- names(r24) <- names(r37) <- names(env)

pred.l <- list(env, r12, r24, r37)


species <- c("lyva", "eclu", "trve")

# Define formulas for each species according to results
forms <- list()
forms[[1]] <- occ ~ sal + sst + I(sst^2)

for (i in 1:3) {
        
        sp <- species[i]
        
        pts <- read.csv(paste0("data/", sp, "/", sp, "_cell.csv"))[,1:2]
        pts <- SpatialPoints(data.frame(x = pts[,1], y = pts[,2]), 
                             proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
        
        # Reproject species data to the Equal Area Lambers projection
        proj <- "+proj=laea +lat_0=0 +lon_0=-70 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs"
        pts <- spTransform(pts, CRS(proj))
        
        ### RandomForest (Down Weighted) ----
        
        # generating background points
        bg <- data.frame(randomPoints(env$sst, n = 20000))
        
        training <- data.frame(
                rbind(cbind(bg, occ = 0),
                      cbind(pts@coords, occ = 1))
        )
        
        # print the first few rows and columns
        head(training)
        
        
        # calculate sub-samples
        prNum <- as.numeric(table(training$occ)["1"]) # number of presence records
        spsize <- c("0" = prNum, "1" = prNum) # sample size for both classes
        
        training$occ <- as.factor(training$occ)
        
        training <- data.frame(cbind(training, raster::extract(env, training[,1:2])))
        
        # RF with down-sampling
        m.rf <- randomForest(occ ~ sal + sst + I(sst^2) + ph,
                               data = training,
                               ntree = 1000,
                               sampsize = spsize,
                               replace = TRUE) # make sure samples are with replacement (default)
        
        # predict to raster layers
        pred.rf <-  lapply(pred.l, predict, m.rf, type = "prob", index = 2)

        
        ### GLM (Down Weighted) ----
        p.wt <- rep(1.e-6, length(occ))
        p.wt[occ == 0] = gArea(starea)/sum(occ == 0)
        m.glm <- glm(Pres/p.wt ~ sal + sst + I(sst^2) + ph, family = poisson(), weights = p.wt)
        
        
        ### Maxent (point-process mode) ----
        m.max <- maxnet(p = occ,
                        data = sp.data,
                        f = maxnet.formula(
                                occ, sp.data, classes = "lqpht"
                        ))
        
        # ver se inclui ENMeval
        
        ### BRT (not point process) ----
        m.dismo <- dismo::gbm.step()
        
        
        ### GLM (not point process) ----
        m.glm.npp <- glm(occ ~ sal + sst + I(sst^2) + ph, family = binomial())
        
}