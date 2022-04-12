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

sel.area <- shapefile("data/env/crop_shape.shp")

# Load bathymetry layer

bath <- raster("data/env/bath_layers/bath_2_100.tif")

#Load layer codes

#codes <- read.table("data/env/layercodes.txt")
codes <- c("BO21_tempmean_ss",
           "BO21_tempmax_ss",
           "BO21_salinitymean_ss",
           "BO21_salinitymin_ss",
           #"BO21_lightbotmean_bdmean",
           #"BO21_lightbotmax_bdmean",
           "BO_ph",
           "BO21_chlomean_ss",
           "BO21_chlomax_ss",
           "BO21_silicatemax_ss",
           "BO21_silicatemean_ss",
           "BO21_dissoxmean_ss",
           "BO21_dissoxmin_ss")

# Load environmental layers

#layer.codes <- as.vector(codes$V1)

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

env.f <- mask(env, bath)
env <- mask(env.f, sel.area)

plot(env$BO21_tempmean_ss)
plot(bath)

names(env)

env <- stack(env, bath)

### Load bioclimatic layers and stack
bioclim <- stack("data/env/bioclim_layers/mtemp_warmq_current_hr.tif",
                 "data/env/bioclim_layers/mtemp_coldq_current_hr.tif")

names(bioclim) <- c("warmest_quarter", "coldest_quarter")

#env <- stack(env, bioclim)

#### Colinearity verification
vifstep.env <- vifstep(env, th = 5)

vifstep.env


### NEW COLINEARITY VERIFICATION - INCLUDING LAYERS THAT WE DECIDED ARE IMPORTANT

#Exclude based on vifstep
env.2 <- exclude(env, vifstep.env)

names(env.2)

#
env.2 <- stack(env.2, env$BO21_salinitymean_ss)

env.2 <- dropLayer(env.2, c("BO21_silicatemean_ss"))

vifstep(env.2, th = 5)
cor(sampleRandom(env.2, 10000))

lyva <- read.csv("data/lyva/lyva_cell.csv")
eclu <- read.csv("data/eclu/eclu_cell.csv")
trve <- read.csv("data/trve/trve_cell.csv")

lyva.env <- extract(env.2, lyva[,1:2])
eclu.env <- extract(env.2, eclu[,1:2])
trve.env <- extract(env.2, trve[,1:2])

vif(data.frame(lyva.env))
vif(data.frame(eclu.env))
vif(data.frame(trve.env))

cor(data.frame(lyva.env))
cor(data.frame(eclu.env))
cor(data.frame(trve.env))



env.3 <- dropLayer(env.2, "BO21_tempmean_ss")

env.3 <- stack(env.3, bioclim)

vif(env.3)
cor(sampleRandom(env.3, 10000))

lyva <- read.csv("data/lyva/lyva_cell.csv")
eclu <- read.csv("data/eclu/eclu_cell.csv")
trve <- read.csv("data/trve/trve_cell.csv")

lyva.env <- extract(env.3, lyva[,1:2])
eclu.env <- extract(env.3, eclu[,1:2])
trve.env <- extract(env.3, trve[,1:2])

vif(data.frame(lyva.env))
vif(data.frame(eclu.env))
vif(data.frame(trve.env))

cor(data.frame(lyva.env))
cor(data.frame(eclu.env))
cor(data.frame(trve.env))

cor(sampleRandom(env.2, 1000))

###Now exclude variations of the same variable based on the ones
### that have stronger biological connection

env.2 <- dropLayer(env.2, "BO21_dissoxmin_ss")

#We can make another vifstep verification
vifstep.env2 <- vifstep(env.2, th = 5)
vifstep.env2

env.2 <- dropLayer(env.2, "bath_2_100")

#Write final list of layers
write.table(names(env.2), "data/env/env_layers.txt", col.names = F)

#Save Vifstep outputs
capture.output(vifstep.env, file = "data/env/vifstep_result.txt")
capture.output(vifstep.env@corMatrix, file = "data/env/vifstep_matrix.txt")
capture.output(vifstep.env2, file = "data/env/vif_result_afterexcluding.txt")
capture.output(vif(env.2), file = "data/env/vif_result_final.txt")

##### Separate rasters
names(env)

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

write.table(names(env), "data/env/layers_names_full.txt", col.names = F)

#### END
