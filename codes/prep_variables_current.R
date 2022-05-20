#### Modelling of coral reef herbivorous invertebrates ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2021

#### Preparing environmental layers #####

library(sdmpredictors)
library(dplyr)
library(raster)
library(sp)
library(usdm)

#Define study extent
ext <- extent(-99, -29, -42.5, 42.5)

# Shapefile of the selected area
# This was created based masking the SST Bio-oRACLE Layer with
# the bathymetry layer + the crop_shape shapefile (see the prep_bath file).
# Then, some "islands" (areas up to 100m) on the borders of the study area
# (e.g.-40 W-39 N) were removed to restrict the study area and improve analysis.
sel.area <- shapefile("gis/starea.shp")

# Load bathymetry layer
bath <- raster("data/env/bath_layers/bath_2_100.tif")

#Load layer codes
codes <- c("BO21_tempmean_ss",
           "BO21_tempmax_ss",
           "BO21_salinitymean_ss",
           "BO21_salinitymin_ss",
           "BO_ph",
           "BO21_chlomean_ss",
           "BO21_chlomax_ss",
           "BO21_silicatemax_ss",
           "BO21_silicatemean_ss",
           "BO21_dissoxmean_ss",
           "BO21_dissoxmin_ss")

# Load environmental layers
env <- load_layers(
  layercodes = codes ,
  equalarea = FALSE,
  rasterstack = FALSE,
  datadir = "data/env/env_layers"
)

# Crop environmental layers to study extent
for (i in 1:length(env)) {
  env[[i]] <- crop(env[[i]], ext)
}

env <- stack(env)

# Mask layers to bathymetry and to selected area
env.b <- mask(env, bath)
env <- mask(env.b, sel.area)

plot(env$BO21_tempmean_ss_lonlat)
plot(bath)

names(env)

# Load bioclimatic layers
bioclim <- stack("data/env/bioclim_layers/mtemp_warmm_current_hr.tif",
                 "data/env/bioclim_layers/mtemp_coldm_current_hr.tif")

names(bioclim) <- c("warmest_month", "coldest_month")

#### Colinearity verification
vifstep.env <- vifstep(env, th = 5)

vifstep.env

# NEW COLINEARITY VERIFICATION - INCLUDING LAYERS THAT WE DECIDED ARE IMPORTANT
#Exclude based on vifstep
env.2 <- exclude(env, vifstep.env)
env.2 <- dropLayer(env.2, c("BO21_silicatemean_ss_lonlat"))
env.2 <- stack(env.2, env$BO21_salinitymean_ss_lonlat)

vifstep(env.2, th = 5)
cor(sampleRandom(env.2, 10000))

# Now consider the occurrence points and a sample of points
lyva <- read.csv("data/lyva/lyva_cell.csv")
eclu <- read.csv("data/eclu/eclu_cell.csv")
trve <- read.csv("data/trve/trve_cell.csv")

lyva.env <- rbind(extract(env.2, lyva[,1:2]),
                  extract(env.2, sampleRandom(env.2, 1000, xy = T)[,1:2]))
eclu.env <- rbind(extract(env.2, eclu[,1:2]),
                  extract(env.2, sampleRandom(env.2, 1000, xy = T)[,1:2]))
trve.env <- rbind(extract(env.2, trve[,1:2]),
                  extract(env.2, sampleRandom(env.2, 1000, xy = T)[,1:2]))

vif(data.frame(lyva.env))
vif(data.frame(eclu.env))
vif(data.frame(trve.env))

cor(data.frame(lyva.env))
cor(data.frame(eclu.env))
cor(data.frame(trve.env))

# Remove dissox min, that is causing a high VIF
env.3 <- dropLayer(env.2, "BO21_dissoxmin_ss_lonlat")

# Run final VIF verification
vif(env.3)

# Now test changing the "tempmean" layer for one of the bioclimatics
bc.test <- env.3
bc.test$BO21_tempmean_ss_lonlat <- bioclim$warmest_month
vif(bc.test)

bc.test$BO21_tempmean_ss_lonlat <- bioclim$coldest_month
vif(bc.test)

#Write final list of layers
write.table(names(env.2), "data/env/env_layers.txt", col.names = F)

#Save Vifstep outputs
capture.output(vifstep.env, file = "data/env/vifstep_result.txt")
capture.output(vifstep.env@corMatrix, file = "data/env/vifstep_matrix.txt")
capture.output(vif(env.3), file = "data/env/vif_result_final.txt")

# Separate rasters and save
names(env) # We save the full files, so if needed it's ready.

if (dir.exists("data/env/crop_layers") == F) {
  dir.create("data/env/crop_layers", recursive = T)
}

writeRaster(
  env,
  filename = paste("data/env/crop_layers/", names(env), sep = ""),
  format = "GTiff",
  bylayer = TRUE,
  overwrite = T
)

write.table(names(env.3), "data/env/selected_env_layers.txt", col.names = F)

#### END
