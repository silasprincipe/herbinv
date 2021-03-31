#### Modelling of Herbivorous invertebrates of coral reefs ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2020

#### Bathymetry Layer preparation #####

### Important Note: to run this code you first need to download the
# GEBCO (General Bathymetric Chart of the Oceans) from www.gebco.net
# and place it in the data/env/bath_layers/ folder.
# Also note that you don't need to run it, as we already provide
# the croped layer.

# Load libraries ----
library(raster)

# Load data and define study parameters ----
# Define study extent
ext <- extent(-99, -29, -45, 45)

# Load GEBCO file
bath <- raster("data/env/bath_layers/gebco_america_africa.tif")

# Load shapefile
sel.area <- shapefile("data/env/crop_shape.shp")

# Prepare the bath layer according to parameters ----
# Crop for the extension
bath.crop <- crop(bath, ext)

plot(bath.crop)

# Restrict according to depth
bath.2_100 <-
  calc(
    bath.crop,
    fun = function(x) {
      x[x > 2 | x < -100] <- NA
      return(x)
    }
  )

# Adjust resolution
bath.2_100 <- aggregate(bath.2_100, fact = 20, fun = mean)

# Mask in the selected area
m2_100 <- mask(bath.2_100, sel.area)

# Plot to verify
plot(m2_100)

# Write final rasters

writeRaster(m2_100, filename = "data/env/bath_layers/bath_2_100.tif",
            overwrite = T)


###END