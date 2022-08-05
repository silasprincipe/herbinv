#### Modelling of coral reef herbivorous invertebrates ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2022 ##

# Plot of SST analysis - all species

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

# Load threshold data ----
load("data/sst_limits/allspecies_oisst_thvalues.RData")

# Load results ----
sp <- "lyva"

# Load rasters generated before
curr <- raster(paste0("data/sst_limits/", sp, "_current_thresh.tif"))
ssp1 <- raster(paste0("data/sst_limits/", sp, "_", "ssp126", "_thresh.tif"))
ssp2 <- raster(paste0("data/sst_limits/", sp, "_", "ssp245", "_thresh.tif"))
ssp3 <- raster(paste0("data/sst_limits/", sp, "_", "ssp370", "_thresh.tif"))

# Get the percentage of time to use as threshold (mean of min and max point)
lval <- round(((thresholds[[sp]]$time_inrange_hottest_point +
                        thresholds[[sp]]$time_inrange_coolest_point)/2),
              2) # round to 2 digits

# Reproject to equal-area lambers
proj <- "+proj=laea +lat_0=0 +lon_0=-70 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs"

crs(ssp1) <- crs(ssp2) <- crs(ssp3) <- crs(curr)

curr <- projectRaster(curr, crs = proj)
ssp1 <- projectRaster(ssp1, crs = proj)
ssp2 <- projectRaster(ssp2, crs = proj)
ssp3 <- projectRaster(ssp3, crs = proj)

# Get the polygons of the areas that are suitable
get.pol <- function(x){
        temp <- calc(x, function(x){
                x[x < lval] <- NA
                x[x >= lval] <- 1
                x
        })
        temp <- st_as_sf(rasterToPolygons(temp, dissolve = T))
        temp <- st_set_crs(temp, crs(curr))
        temp
}

curr.p <- get.pol(curr)
ssp1.p <- get.pol(ssp1)
ssp2.p <- get.pol(ssp2)
ssp3.p <- get.pol(ssp3)

# Convert to data.frame
get.val <- function(x){
        temp <- as(x, "SpatialPixelsDataFrame")
        temp <- as.data.frame(temp)
        colnames(temp) <- c("val", "x", "y")
        temp$val <- temp$val * 100
        temp
}

curr.v <- get.val(curr)
ssp1.v <- get.val(ssp1)
ssp2.v <- get.val(ssp2)
ssp3.v <- get.val(ssp3)

# Prepare theme/ploting stuff ----
# Guide for legend
step.guide <- guide_coloursteps(title = "% time",
                                show.limits = TRUE,
                                barheight = unit(2.4, "in"),
                                barwidth = unit(0.18, "in"),
                                ticks = T,
                                ticks.colour = "grey20",
                                frame.colour = "grey20",
                                title.vjust = unit(0.1, "in"))

# Themes
nlt <- theme_classic()+
        theme(axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_text(colour = "grey60"),
              axis.title.x.top = element_text(vjust = 2),
              legend.position="none",
              panel.grid.major = element_line(linetype = 'dashed', 
                                              colour = "grey70",
                                              size = .05),
              panel.background = element_blank(),
              plot.background = element_rect(color = "grey30", fill = 'white')
        )

wlt <- theme_classic()+
        theme(axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_text(colour = "grey60"),
              axis.title.x.top = element_text(vjust = 2),
              panel.background = element_blank(),
              panel.grid.major = element_line(linetype = 'dashed', 
                                              colour = "grey70",
                                              size = .1),
              legend.justification = c(1, 0),
              legend.position = c(1.2, 0),
              legend.background = element_rect(fill = "white")
        )

# Scale of values
sca <- scale_fill_stepsn(breaks = seq(10,90,10),
                         limits = c(0, 100),
                         colors = RColorBrewer::
                                 brewer.pal(10, "Spectral"),
                         na.value = NA,
                         guide = step.guide)



# Generate plots ----
# Current
(pc <- ggplot()+
         # Base maps
         geom_sf(data = base,
                 color = c("#CFCFCF"),
                 size = 0.6,
                 fill = "white") +
         # Get the raster
         geom_raster(data = curr.v, aes(x = x, y = y, fill = val)) +
         sca + # Add color scale
         # Add the suitable area polygon
         geom_sf(data = curr.p, color = "blue", fill = NA) +
         # Establish area
         coord_sf(xlim = c(-3239.984, 4510.636),
                  ylim = c(-4614.575, 4596.965),
                  datum = st_crs(base),
                  expand = F,
                  label_axes = list(
                          top = "E",
                          left = "N",
                          top = ""
                  )) +
         # Add title
         geom_label(aes(label = "A", x = 4000, y = 4100),
                    size = 10, fontface = "bold", color = "grey30",
                    label.size = 0)+
         # Remove axis labels and add theme
         xlab("Easting") + ylab("Northing") + wlt +
         # Add grid
         scale_x_discrete(position = "top") 
)

# SSP1
(p1 <- ggplot()+
                # Base maps
                geom_sf(data = base,
                        color = c("#CFCFCF"),
                        size = 0.6,
                        fill = "white") +
                # Get the raster
                geom_raster(data = ssp1.v, aes(x = x, y = y, fill = val)) +
                sca + # Add color scale
                # Add the suitable area polygon
                geom_sf(data = ssp1.p, color = "blue", fill = NA) +
                # Establish area
                coord_sf(xlim = c(-3239.984, 4510.636),
                         ylim = c(-4614.575, 4596.965),
                         datum = st_crs(base),
                         expand = F) +
                # Add title
                geom_label(aes(label = "B", x = 4100, y = 4200),
                  size = 10, fontface = "bold", color = "grey30",
                  label.size = 0)+
                # Remove axis labels and add theme
                xlab(NULL) + ylab(NULL) + nlt +
                theme(plot.margin = margin(1, 1, 0, 0, "pt"),
                      plot.background = 
                              element_rect(color = "grey30", fill = 'white'),
                      axis.text = element_blank(),
                      panel.grid.major = element_blank()
                      )
)

# SSP2
(p2 <- ggplot()+
                # Base maps
                geom_sf(data = base,
                        color = c("#CFCFCF"),
                        size = 0.6,
                        fill = "white") +
                # Get the raster
                geom_raster(data = ssp2.v, aes(x = x, y = y, fill = val)) +
                sca + # Add color scale
                # Add the suitable area polygon
                geom_sf(data = ssp2.p, color = "blue", fill = NA) +
                # Establish area
                coord_sf(xlim = c(-3239.984, 4510.636),
                         ylim = c(-4614.575, 4596.965),
                         datum = st_crs(base),
                         expand = F) +
                # Add title
                geom_label(aes(label = "C", x = 4100, y = 4200),
                           size = 10, fontface = "bold", color = "grey30",
                           label.size = 0)+
                # Remove axis labels and add theme
                xlab(NULL) + ylab(NULL) + nlt +
                theme(plot.margin = margin(1, 1, 0, 0, "pt"),
                      plot.background = 
                              element_rect(color = "grey30", fill = 'white'),
                      axis.text = element_blank(),
                      panel.grid.major = element_blank()
                )
)

# SSP3
(p3 <- ggplot()+
                # Base maps
                geom_sf(data = base,
                        color = c("#CFCFCF"),
                        size = 0.6,
                        fill = "white") +
                # Get the raster
                geom_raster(data = ssp3.v, aes(x = x, y = y, fill = val)) +
                sca + # Add color scale
                # Add the suitable area polygon
                geom_sf(data = ssp3.p, color = "blue", fill = NA) +
                # Establish area
                coord_sf(xlim = c(-3239.984, 4510.636),
                         ylim = c(-4614.575, 4596.965),
                         datum = st_crs(base),
                         expand = F) +
                # Add title
                geom_label(aes(label = "D", x = 4100, y = 4200),
                           size = 10, fontface = "bold", color = "grey30",
                           label.size = 0)+
                # Remove axis labels and add theme
                xlab(NULL) + ylab(NULL) + nlt +
                theme(plot.margin = margin(1, 1, 0, 0, "pt"),
                      plot.background = 
                              element_rect(color = "grey30", fill = 'white'),
                      axis.text = element_blank(),
                      panel.grid.major = element_blank()
                )
)

# Get final plot ----
(fp <- p1 + p2 + p3)

ggsave(paste0("figures/", sp, "_sst_composite.jpg"), fp,
       width = 48, height = 18, units = "cm")

ggsave(paste0("figures/", sp, "_sst_mainmap.jpg"), pc,
       width = 21, height = 18, units = "cm")
