#### Modelling of coral reef herbivorous invertebrates ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2022 ##

# Plot of LGCP results - all species

# Load needed packages ----
library(ggplot2)
library(sf)
library(raster)
library(patchwork)

# Load base shapefiles ----
base <- shapefile("gis/basemaps/ne_110m_land_edited.shp")
# Crop to the extent
base <- buffer(base, 0)
base <- crop(base, extent(-120, -10, -55, 65))
# Reproject
proj <- "+proj=laea +lat_0=0 +lon_0=-70 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs"
base <- spTransform(base, CRS(proj))
# Convert to SF
base <- st_as_sf(base)

# Load results ----
sp <- "lyva"

curr <- raster(paste0("results/", sp, "/predictions/", sp, "_mean_current.tif"))

ssp1 <- raster(paste0("results/", sp, "/predictions/", sp, "_mean_ssp1.tif"))
ssp2 <- raster(paste0("results/", sp, "/predictions/", sp, "_mean_ssp2.tif"))
ssp3 <- raster(paste0("results/", sp, "/predictions/", sp, "_mean_ssp3.tif"))

# Get differences
d.ssp1 <- (ssp1 - curr) / curr
d.ssp2 <- (ssp2 - curr) / curr
d.ssp3 <- (ssp3 - curr) / curr

# Convert to data.frame
get.val <- function(x){
        temp <- as(x, "SpatialPixelsDataFrame")
        temp <- as.data.frame(temp)
        colnames(temp) <- c("val", "x", "y")
        temp
}

ssp1.v <- get.val(d.ssp1)
ssp2.v <- get.val(d.ssp2)
ssp3.v <- get.val(d.ssp3)

# Prepare theme/ploting stuff ----
# Guide for legend
step.guide <- guide_coloursteps(title = "Change in intensity",
                                show.limits = TRUE,
                                barheight = unit(0.12, "in"),
                                barwidth = unit(3.5, "in"),
                                ticks = T,
                                ticks.colour = "grey20",
                                frame.colour = "grey20",
                                title.position = "top")

# Dvision of multiple maps
mlines <- data.frame(x = c(1256.5727, 1256.5727, 4890.4723, 4890.4723),
                    xend = c(1256.5727, 4890.4723,4890.4723, 2325.4219),
                    y = c(4470.2438, 1996.0985, -373.9687, -993.9687),
                    yend = c(1996.0985, -373.9687,-993.9687, -4848.4605))

# Themes
nlt <- theme_classic()+
        theme(axis.line = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              legend.position="none",
              legend.title.align=0.5,
              panel.background = element_blank(),
              plot.background = element_blank()
        )

wlt <- theme_classic()+
        theme(axis.line = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              legend.position="bottom",
              legend.title.align=0.5,
              legend.spacing.x = unit(0, 'cm'),
              panel.background = element_blank(),
              plot.background = element_blank()
        )

# Scale of values
sca <- scale_fill_stepsn(breaks = c(-10, -5, -2.5, -1, 0,
                                    1, 2.5, 5, 10),
                         limits = c(-20, 20),
                         values = scales::rescale(c(-20,-10, -5, -2.5, -1, 0,
                                                    1, 2.5, 5, 10, 20)),
                         colors = RColorBrewer::
                                 brewer.pal(11, "BrBG"),
                         na.value = NA,
                         guide = step.guide)

sca <- scale_fill_manual(
        breaks = c(-1, 0, 1, 5),
        values = c("blue", "grey", "green", "yellow"),
        guide = guide_legend(
                title = "Change in intensity",
                keyheight = unit(0.12, "in"),
                keywidth = unit(3.5, "in"),
                label.position = "bottom",
                ticks = T,
                ticks.colour = "grey20",
                frame.colour = "grey20",
                title.position = "top"
        )
        
)



# Generate plots ----
# SSP1
(p1 <- ggplot()+
                # Base maps
                geom_sf(data = base,
                        color = NA,
                        size = 0.6,
                        fill = NA) +
                # Get the raster
                geom_tile(data = ssp1.v, aes(x = x, y = y, fill = val)) +
                sca + # Add color scale
                # Establish area
                coord_sf(xlim = c(-3239.984, 5200.636),
                         ylim = c(-4614.575, 4596.965),
                         datum = st_crs(base),
                         expand = F) +
                # Add division lines
                geom_segment(data = mlines,
                             aes(x = x, xend = xend, y = y, yend = yend),
                             size = 2, color = "grey90",
                             lineend = "round")+
                # Add title
                geom_text(aes(label = "SSP1", x = -2700, y = 4400),
                  size = 10, fontface = "bold", color = "grey30")+
                # Remove axis labels and add theme
                xlab(NULL) + ylab(NULL) + nlt
                )

# SSP2
(p2 <- ggplot()+
                # Base maps
                geom_sf(data = base,
                        color = NA,
                        size = 0.6,
                        fill = NA) +
                # Get the raster
                geom_tile(data = ssp2.v, aes(x = x, y = y, fill = val)) +
                sca + # Add color scale
                # Establish area
                coord_sf(xlim = c(-3239.984, 5200.636),
                         ylim = c(-4614.575, 4596.965),
                         datum = st_crs(base),
                         expand = F) +
                # Add division lines
                geom_segment(data = mlines,
                             aes(x = x, xend = xend, y = y, yend = yend),
                             size = 2, color = "grey90",
                             lineend = "round")+
                # Add title
                geom_text(aes(label = "SSP2", x = -2700, y = 4400),
                          size = 10, fontface = "bold", color = "grey30")+
                # Remove axis labels and add theme
                xlab(NULL) + ylab(NULL) + wlt
)

# SSP3
(p3 <- ggplot()+
                # Base maps
                geom_sf(data = base,
                        color = NA,
                        size = 0.6,
                        fill = NA) +
                # Get the raster
                geom_tile(data = ssp3.v, aes(x = x, y = y, fill = val)) +
                sca + # Add color scale
                # Establish area
                coord_sf(xlim = c(-3239.984, 5200.636),
                         ylim = c(-4614.575, 4596.965),
                         datum = st_crs(base),
                         expand = F) +
                # Add title
                geom_text(aes(label = "SSP3", x = -2700, y = 4400),
                          size = 10, fontface = "bold", color = "grey30")+
                # Remove axis labels and add theme
                xlab(NULL) + ylab(NULL) + nlt
)

# Get final plot ----

# Creates a layout
layout <- c(
        patchwork::area(t = 1.1, l = 1  , b = 5.5, r = 4),
        patchwork::area(t = 1  , l = 3  , b = 5  , r = 6),
        patchwork::area(t = 1.1, l = 5  , b = 5.5  , r = 8)
)

# Get final plot
(fp <- p1 + p2 + p3 + plot_layout(design = layout))

# save
