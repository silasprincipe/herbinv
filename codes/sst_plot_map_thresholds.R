#### Modelling of coral reef herbivorous invertebrates ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2021

### Create maps of temperature limits for species of sea-urchins ###
# Based on data extracted from OISST (current) and GFDL-NOAA (future)
# Methods are based on the article of Castro et al. (2020)
# DOI: 10.1111/ddi.13073

# Load libraries ----
library(raster)
library(tidyverse)
library(patchwork)
library(sf)

# Load threshold data ----
load("data/sst_limits/allspecies_oisst_thvalues.RData")

# Establish species and plot details ----

# List species
species <- c("lyva", "eclu", "trve")

## PLOT DETAILS
# Establish ggplot theme details
main.theme <-
        theme(
                panel.background = element_rect(fill = "white", 
                                                colour = "lightgrey"),
                panel.grid.major = element_line(linetype = 'solid', 
                                                colour = "lightgrey"),
                #panel.grid.major = element_blank(), #If wanted to remove line grid
                axis.ticks.length = unit(0.8, "mm"),
                axis.text = element_text(size = 10),
                legend.justification = c(1, 0),
                legend.position = c(0.18, 0.02),
                legend.background = element_blank(),
                # legend.title = element_text(size = 10),
                # legend.text = element_text(size = 10),
                # legend.key.height = unit(7, 'mm'),
                # legend.key.width = unit(7, 'mm'),
                # legend.key = element_rect(colour = "#696969", size =
                #                                   2.5),
                plot.margin = unit(c(0.3, 0.35, 0.3, 0.35), 'mm')
        )

# Load base maps
# Pacific ocean
base2 <- shapefile("gis/basemaps/ne_110m_land_no_oc.shp")
# America
base1 <- shapefile("gis/basemaps/ne_110m_land_edited.shp")
# Convert to SF
base2 <- st_as_sf(base2)
base1 <- st_as_sf(base1)

#Base maps colors
c.sea <- "grey90"
c.land <- "grey70"

# Maps scales
main.s.x <-
        scale_x_continuous(breaks = seq(-35,-100,-10),
                           limits = c(-100,-29))
main.s.y <-
        scale_y_continuous(breaks = seq(-40, 40, 10),
                           limits = c(-42.5, 42.5))

step.guide <- guide_coloursteps(title = "% time",
                                show.limits = TRUE,
                                barheight = unit(2.4, "in"),
                                barwidth = unit(0.18, "in"),
                                ticks = T,
                                ticks.colour = "grey20",
                                frame.colour = "grey20",
                                title.vjust = unit(0.1, "in"))
###

# Run in loop ----
for (i in 1:length(species)){

    # Load rasters generated on GEE
    curr <- raster(paste0("data/sst_limits/", species[i], "_current_thresh.tif"))
    fut <- raster(paste0("data/sst_limits/", species[i], "_", ssp, "_thresh.tif"))

    # Get the percentage of time to use as threshold (mean of min and max point)
    lval <- round(((thresholds[[species[i]]]$time_inrange_hottest_point +
                        thresholds[[species[i]]]$time_inrange_coolest_point)/2),
                  2) # round to 2 digits

    
    # Get the polygons of the areas that are suitable
    curr.p <- curr

    curr.p <- calc(curr.p, function(x){
        x[x < lval] <- NA
        x[x >= lval] <- 1
        x
    })

    curr.p <- st_as_sf(rasterToPolygons(curr.p, dissolve = T))
    
    fut.p <- fut

    fut.p <- calc(fut.p, function(x){
        x[x < lval] <- NA
        x[x >= lval] <- 1
        x
    })

    fut.p <- st_as_sf(rasterToPolygons(fut.p, dissolve = T))
    
    fut.p <- st_set_crs(fut.p, crs(curr.p))
    

    # Convert to data-frame
    curr <- data.frame(cbind(xyFromCell(curr, 1:ncell(curr)), values(curr)))
    colnames(curr) <- c("x", "y", "val")
    curr$val <- curr$val * 100

    fut <- data.frame(cbind(xyFromCell(fut, 1:ncell(fut)), values(fut)))
    colnames(fut) <- c("x", "y", "val")
    fut$val <- fut$val * 100

    # Generate current map
    current.pl <-   ggplot()+
                    # Base maps
                    geom_sf(data = base2,
                            color = c.sea,
                            fill = c.sea) +
                    geom_sf(data = base1,
                            color = c.land,
                            fill = c.land) +
                    # Get the raster
                    geom_raster(data = curr, aes(x = x, y = y, fill = val))+
                    scale_fill_stepsn(breaks = seq(10,90,10),
                                      limits = c(0, 100),
                                      colors = RColorBrewer::
                                          brewer.pal(10, "Spectral"),
                                      na.value = NA,
                                      guide = step.guide)+
                    # Add the suitable area
                    geom_sf(data = curr.p,
                            color = "blue",
                            fill = NA)+
                    # Establish area
                    coord_sf(
                            xlim = c(-99, -29),
                            ylim = c(-42.5, 42.5),
                            expand = FALSE,
                    ) +
                    # Add plot details
                    xlab(NULL) + ylab(NULL)+
                    main.theme + main.s.x + main.s.y

    # Generate future map
    future.pl <-   ggplot()+
                    # Base maps
                    geom_sf(data = base2,
                            color = c.sea,
                            fill = c.sea) +
                    geom_sf(data = base1,
                            color = c.land,
                            fill = c.land) +
                    # Get the raster
                    geom_raster(data = fut, aes(x = x, y = y, fill = val))+
                    scale_fill_stepsn(breaks = seq(10,90,10),
                                      limits = c(0, 100),
                                      colors = RColorBrewer::
                                          brewer.pal(10, "Spectral"),
                                      na.value = NA,
                                      guide = step.guide)+
                    # Add the suitable area
                    geom_sf(data = fut.p,
                            color = "blue",
                            fill = NA)+
                    # Establish area
                    coord_sf(
                            xlim = c(-99, -29),
                            ylim = c(-42.5, 42.5),
                            expand = FALSE,
                            label_axes = list(
                                   bottom = "E",
                                   right = "N",
                                   top = ""
                           )
                    ) +
                    # Add plot details
                    xlab(NULL) + ylab(NULL)+
                    main.theme + main.s.x + main.s.y +
                    theme(legend.position = 'none')

    # Save the plot
    map <- current.pl + future.pl

    ggsave(paste0("figures/", species[i], "_temp_thresh.tiff"),
        map,
        width = 25, height = 15, dpi = 300, units = 'cm')

}
### END