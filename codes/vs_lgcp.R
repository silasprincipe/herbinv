#### Modelling of Herbivorous invertebrates of coral reefs ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2021

### Generating virtual species for testing methods ###

# For establishing the best methods for the modeling exercise, we first created
# two virtual species with caracteristics that are similar to the sea-urchin
# species being studied here.
#
# This enables us to test four things:
# 1) The method chosen to reduce sampling bias and spatial autocorrelation
# 2) The method chosen to generate pseudo-absences
# 3) The ability of the chosen algorithms in capturing:
#       a. the expected response curves
#       b. the expected direction of change in future climates

# Load libraries ----
library(virtualspecies)
library(raster)
source("functions/Varload.r")

# Set seed for replicability
set.seed(2932)

# Load environmental layers ----
env <- var.load(folder = "crop_layers/", "env_layers.txt", "2_100")

r12 <- var.load(folder = "proj_layers/ssp126/", "env_layers.txt", "2_100")

r24 <- var.load(folder = "proj_layers/ssp245/", "env_layers.txt", "2_100")

r37 <- var.load(folder = "proj_layers/ssp370/", "env_layers.txt", "2_100")

r85 <- var.load(folder = "proj_layers/ssp585/", "env_layers.txt", "2_100")

# Set parameters for the virtual species ----

# The sp1 represents a species with thermal limits closest to L. variegatus
# and E. lucunter
vsp.param1 <- formatFunctions(
        BO_ph = c(fun = 'dnorm',
                  mean = 8.1,
                  sd = 0.2),
        BO21_chlomean_ss = c(fun = 'logisticFun',
                             beta = 0.8,
                             alpha = 0.1),
        # BO21_silicatemean_ss = c(fun = 'logisticFun',
        #                          beta = 25,
        #                          alpha = 1),
        BO21_tempmean_ss = c(fun = 'custnorm',
                             mean = 25,
                             diff = 10,
                             prob = 0.99),
        BO21_salinitymean_ss = c(fun = 'dnorm',
                                 mean = 35,
                                 sd = 10),
        bath_2_100 = c(fun = 'dnorm',
                       mean = -20,
                       sd = 20))

# The sp2 represents a species with thermal limits closest to T. ventricosus
vsp.param2 <- formatFunctions(
        BO_ph = c(fun = 'dnorm',
                  mean = 8.1,
                  sd = 0.2),
        BO21_chlomean_ss = c(fun = 'logisticFun',
                             beta = 0.8,
                             alpha = 0.1),
        # BO21_silicatemean_ss = c(fun = 'logisticFun',
        #                          beta = 25,
        #                          alpha = 1),
        BO21_tempmean_ss = c(fun = 'custnorm',
                             mean = 26.5,
                             diff = 7.5,
                             prob = 0.99),
        BO21_salinitymean_ss = c(fun = 'dnorm',
                                 mean = 35,
                                 sd = 10),
        bath_2_100 = c(fun = 'dnorm',
                       mean = -20,
                       sd = 20))

# Generate virtual species ----
vsp1 <- generateSpFromFun(env, vsp.param1, plot = T)
vsp2 <- generateSpFromFun(env, vsp.param2, plot = T)

# Plot response curves to verify
plotResponse(vsp1)
plotResponse(vsp2)

# Convert to presence absence
vsp1.pa <- convertToPA(vsp1,
                       PA.method = "probability",
                       prob.method = "logistic",
                       beta = 0.5, alpha = -0.05,
                       plot = TRUE)

vsp2.pa <- convertToPA(vsp2,
                       PA.method = "probability",
                       prob.method = "logistic",
                       beta = 0.5, alpha = -0.05,
                       plot = TRUE)

# Plot to verify
par(mfrow = c(1, 2))
plot(vsp1.pa$pa.raster, legend = F, main = "VSP 1")
plot(vsp2.pa$pa.raster, legend = F, main = "VSP 2")


# Sample occurrence points ----

# Load bias polygon
bias.polygon <- shapefile("gis/vsp_bias_polygon.shp")

# SP 1 sampling
# High number of points sampling:
vsp1.pns.hn <- sampleOccurrences(vsp1.pa, n = 300)

# Low number of points sampling
vsp1.pns.ln <- sampleOccurrences(vsp1.pa, n = 150)

# High number of points sampling with polygon bias:
vsp1.pbs.hn <- sampleOccurrences(vsp1.pa, n = 300,
                                 bias = "polygon",
                                 bias.strength = 10,
                                 bias.area = bias.polygon)

# Low number of points sampling with polygon bias:
vsp1.pbs.ln <- sampleOccurrences(vsp1.pa, n = 150,
                                 bias = "polygon",
                                 bias.strength = 10,
                                 bias.area = bias.polygon)


# SP 2 sampling - in that case, only low number (80), the same as T. ventricosus
vsp2.pns <- sampleOccurrences(vsp2.pa, n = 80)

vsp2.pbs <- sampleOccurrences(vsp2.pa, n = 80,
                              bias = "polygon",
                              bias.strength = 10,
                              bias.area = bias.polygon)

vsp2.psbs <- sampleOccurrences(vsp2.pa, n = 80,
                               bias = "manual",
                               bias.strength = 10,
                               weights = vsp1.pa$suitab.raster)



# Simulate future conditions ----
vsp1.r12 <- generateSpFromFun(r12, vsp.param1, plot = T)
vsp1.r24 <- generateSpFromFun(r24, vsp.param1, plot = T)
vsp1.r37 <- generateSpFromFun(r37, vsp.param1, plot = T)
vsp1.r85 <- generateSpFromFun(r85, vsp.param1, plot = T)

vsp2.r12 <- generateSpFromFun(r12, vsp.param2, plot = T)
vsp2.r24 <- generateSpFromFun(r24, vsp.param2, plot = T)
vsp2.r37 <- generateSpFromFun(r37, vsp.param2, plot = T)
vsp2.r85 <- generateSpFromFun(r85, vsp.param2, plot = T)

# Plot response curves to verify (just two, to get an idea)
plotResponse(vsp1.r12)
plotResponse(vsp1.r85)

# Convert to PA
vsp1.fut.pa <- lapply(list(vsp1.r12, vsp1.r24, vsp1.r37, vsp1.r85),
                      convertToPA, PA.method = "probability",
                      prob.method = "logistic",
                      beta = 0.5, alpha = -0.05)

vsp2.fut.pa <- lapply(list(vsp2.r12, vsp2.r24, vsp2.r37, vsp2.r85),
                      convertToPA, PA.method = "probability",
                      prob.method = "logistic",
                      beta = 0.5, alpha = -0.05)

names(vsp2.fut.pa) <- names(vsp1.fut.pa) <- c("r12", "r24", "r37", "r85")



# Save data ----
if(!dir.exists("data/vsp/")){dir.create("data/vsp/")}

# Save vsp 1 data
vsp1.data <- list(vsp1, vsp1.pa,
                  list("normal_hn" = vsp1.pns.hn, "normal_ln" = vsp1.pns.ln,
                       "bias_hn" = vsp1.pbs.hn, "bias_ln" = vsp1.pbs.ln),
                  list("r12" = vsp1.r12, "r24" = vsp1.r24,
                       "r37" = vsp1.r37, "r85" = vsp1.r85),
                  vsp1.fut.pa)

names(vsp1.data) <- c("vsp_models", "vsp_pa", "vsp_points_sampling", 
                      "vsp_fut_models", "vsp_fut_pa")

saveRDS(vsp1.data, file = "data/vsp/vsp1.rds")


# Save vsp 2 data
vsp2.data <- list(vsp2, vsp2.pa,
                  list("normal" = vsp2.pns, "bias" = vsp2.pbs),
                  list("r12" = vsp2.r12, "r24" = vsp2.r24,
                       "r37" = vsp2.r37, "r85" = vsp2.r85),
                  vsp2.fut.pa)

names(vsp2.data) <- names(vsp1.data)

saveRDS(vsp2.data, file = "data/vsp/vsp2.rds")


# Save essential data in other formats

# PA projections rasters
writeRaster(vsp1.pa$pa.raster, "data/vsp/vsp1_pa.tif")
writeRaster(vsp2.pa$pa.raster, "data/vsp/vsp2_pa.tif")

lapply(c("r12", "r24", "r37", "r85"), function(x){
        writeRaster(vsp1.fut.pa[[x]]$pa.raster,
                    paste0("data/vsp/vsp1_",x,"_pa.tif"))
})

lapply(c("r12", "r24", "r37", "r85"), function(x){
        writeRaster(vsp2.fut.pa[[x]]$pa.raster,
                    paste0("data/vsp/vsp2_",x,"_pa.tif"))
})

# Point sampling csv
write.csv(vsp1.pns.hn$sample.points, "data/vsp/vsp1_ptnormalsamp_hn.csv",
          row.names = F)
write.csv(vsp1.pns.ln$sample.points, "data/vsp/vsp1_ptnormalsamp_ln.csv",
          row.names = F)
write.csv(vsp1.pbs.hn$sample.points, "data/vsp/vsp1_ptpolbiassamp_hn.csv",
          row.names = F)
write.csv(vsp1.pbs.ln$sample.points, "data/vsp/vsp1_ptpolbiassamp_ln.csv",
          row.names = F)

write.csv(vsp2.pns$sample.points, "data/vsp/vsp2_ptnormalsamp.csv",
          row.names = F)
write.csv(vsp2.pbs$sample.points, "data/vsp/vsp2_ptpolbiassamp.csv",
          row.names = F)

#rm(list = ls())
### END
