
library(sf) 
library(ncdf4)
library(raster)
library(rasterVis)
library(RColorBrewer)

# rasterVis plot parameters
mapTheme <- rasterTheme(region = rev(brewer.pal(10, "RdBu")))
cutpts <- c(-30,-20,-10,0,10,20,30,40,50)

# OISST data
noaa.oi <- list.files("data/temperature/noaa_oisst/", full.names = T)

oisst <- raster(noaa.oi, varname = 'sst')

oisst <- rotate(oisst)

(plt <- levelplot(subset(oisst, 1), margin = F, at=cutpts, cuts=11, pretty=TRUE, par.settings = mapTheme,
                 main="OI SST"))


# OISST data
noaa.whoi <- list.files("data/temperature/noaa_whoi/", full.names = T)

whoi <- raster(noaa.whoi, varname = 'sea_surface_temperature')

whoi <- rotate(whoi)

(plt <- levelplot(subset(whoi, 1), margin = F, at=cutpts, cuts=11, pretty=TRUE, par.settings = mapTheme,
                  main="WHOI SST"))


# OISST data
noaa.pathf <- list.files("data/temperature/noaa_pathfinder/", full.names = T)

pathf <- raster(noaa.pathf[1], varname = 'sea_surface_temperature')

whoi <- rotate(whoi)

(plt <- levelplot(subset(whoi, 1), margin = F, at=cutpts, cuts=11, pretty=TRUE, par.settings = mapTheme,
                  main="WHOI SST"))



# Explore data from GEE - Pathfinder dataset and Lytechinus variegatus

lyva.gee <- read.csv("data/temperature/gee_data/lyva_datatable_oisst.csv")

lyva.gee$area <- lyva.gee$area * 0.000001

plot(lyva.gee$system.index, lyva.gee$area)

(ggplot(data=lyva.gee, aes(x=system.index, y=area, group=1)) +
                geom_line()+
                geom_smooth(method = 'lm')+
                theme_classic())
