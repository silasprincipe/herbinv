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
           "BO21_lightbotmean_bdmax",
           "BO21_lightbotmax_bdmax",
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

#### Colinearity verification
vifstep.env <- vifstep(env, th = 10)

vifstep.env


### NEW COLINEARITY VERIFICATION - INCLUDING LAYERS THAT WE DECIDED ARE IMPORTANT

#Exclude based on vifstep
env.2 <- exclude(env, vifstep.env)

names(env.2)

###Now exclude variations of the same variable based on the ones
### that have stronger biological connection

env.3 <- dropLayer(env.2,
                   c(
                     "BO21_salinitymean_ss"
                     ))

names(env.3)

#We can make another vifstep verification
vifstep.env3 <- vifstep(env.3, th = 10)
vifstep.env3

#Write final list of layers
write.table(names(env.3), "data/env/env_layers.txt", col.names = F)

#Save Vifstep outputs
capture.output(vifstep.env, file = "data/env/vifstep_result.txt")
capture.output(vifstep.env@corMatrix, file = "data/env/vifstep_matrix.txt")
capture.output(vifstep.env3, file = "data/env/vif_result_afterexcluding.txt")


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
