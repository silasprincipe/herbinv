#Exploratory variables + points

source("functions/Varload.r")
env <- var.load(folder = "crop_layers/", layers = "env_layers.txt",
                bath = "2_100")


lyva <- read.csv("data/lyva/lyva_cell.csv")
coordinates(lyva) <- ~decimalLongitude+decimalLatitude

library(rasterVis)
library(RColorBrewer)

# rasterVis plot parameters
mapTheme <- rasterTheme(region = rev(brewer.pal(10, "RdBu")))

plot.raster <- env$bath_2_100

cutpts <- c(seq(from = minValue(plot.raster),
                to = maxValue(plot.raster), by = 20))
(plt <- levelplot(subset(plot.raster, 1), margin = F, at=cutpts, cuts=11,
                  pretty=TRUE, par.settings = mapTheme,
                  main = names(plot.raster))+
                layer(sp.points(lyva, pch=20, cex=0.5, col='green')))


lyva.env.data <- data.frame(extract(env, lyva))

env.data <- data.frame(extract(env, xyFromCell(env, 1:ncell(env))))
env.data <- env.data[!is.na(env.data$BO2_tempmax_ss),]

lyva.env.data$group <- "lyva"
env.data$group <- "env"

full.data <- rbind(lyva.env.data, env.data) 

(ggplot(full.data, aes(color = group))+
                geom_density(aes(x = bath_2_100), size = 1))


