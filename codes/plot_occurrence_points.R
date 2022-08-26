#### Modelling of coral reef herbivorous invertebrates ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2022 ##

# Plot of the occurrence points

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

# Load data ----
species <- c("lyva", "eclu", "trve")

pts <- lapply(species, function(x){
  read.csv(paste0("data/", x, "/", x, "_filt.csv"))[,1:2]
})

# Load SST
sst <- raster("data/env/crop_layers/BO21_tempmean_ss.tif")
baseproj <- raster("data/env/ready_layers/sst_cur.tif")
sst <- projectRaster(sst, crs = crs(baseproj))

# Convert to data.frame
get.val <- function(x){
  temp <- as(x, "SpatialPixelsDataFrame")
  temp <- as.data.frame(temp)
  colnames(temp) <- c("val", "x", "y")
  return(temp)
}

sst.v <- get.val(sst)

# Prepare theme/ploting stuff ----
# Guide for legend
step.guide <- guide_colorbar(title = "Mean SST (°C)",
                             show.limits = TRUE,
                             barheight = unit(0.12, "in"),
                             barwidth = unit(3.5, "in"),
                             ticks = F,
                             ticks.colour = "grey20",
                             frame.colour = "grey20",
                             title.position = "top")

nlt <- theme_classic()+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x.top = element_text(vjust = 2, size = 16),
        axis.title.y.left = element_text(size = 16),
        legend.position="none",
        panel.grid.major = element_line(linetype = 'dashed', 
                                        colour = "grey70",
                                        size = .05),
        panel.background = element_blank()#,
        #plot.background = element_rect(color = "grey30", fill = 'white')
  )

wlt <- theme_classic()+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(colour = "grey60", size = 11),
        axis.title.x.top = element_text(vjust = 2, size = 11),
        axis.title.y.left = element_text(size = 11),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey60"),
        panel.grid.major = element_line(linetype = 'dashed', 
                                        colour = "grey70",
                                        size = .1),
        legend.position="bottom",
        legend.title.align=0.5,
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.background = element_rect(fill = "white"),
        
  )

# Theme for the density plot
theme_dplot <- theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.ticks.y = element_blank(),
        plot.margin = margin(0,0,0,0),
        plot.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey60"))

int <- theme_classic()+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position="none",
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(color = "grey30", fill = 'white'),
        plot.margin = margin(0.5, 0.5, 0, 0, "pt")
  )

sca <- scale_fill_stepsn(breaks = seq(8,32,2),
                         limits = c(8, 32),
                         colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
                         na.value = NA,
                         guide = step.guide)



# Generate plots ----
### Lyva ----
# Rectangles to highlight
rects <- list(data.frame(x1 = -1600, x2 = -400, y1 = 2600, y2= 3600, label = "a"),
              data.frame(x1 = 2400, x2 = 3200, y1 = -2800, y2= -2200, label = "b"))

(ptsl <- ggplot()+
    # Base maps
    geom_sf(data = base,
            color = c("#CFCFCF"),
            size = 0.6,
            fill = "white") +
    # Get the raster
    geom_raster(data = sst.v, aes(x = x, y = y, fill = val)) +
    sca + # Add color scale
    # Add occurrence points
    geom_point(data = pts[[1]],
               aes(x = x, y = y), size = 1,
               alpha = .5, shape = 16)+
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
    # Draw rectangles
    geom_rect(data = rects[[1]],
              aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
              fill = NA, color = "grey20", size = .3)+
    # geom_text(data = rects[[1]],
    #            aes(x = (x2 - 100), y = (y2 - 100), label = label))+
    geom_rect(data = rects[[2]],
              aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
              fill = NA, color = "grey20", size = .3)+
    # geom_text(data = rects[[1]],
    #            aes(x = (x2 - 100), y = (y2 - 100), label = label))+
    # Add title
    geom_label(aes(label = "A", x = 4000, y = 4100),
               size = 8, fontface = "bold", color = "grey30",
               label.size = 0)+
    # Remove axis labels and add theme
    xlab("Easting") + ylab("Northing") + wlt +
    # Add grid
    scale_x_discrete(position = "top", breaks = c(-2000, 0, 2000, 4000)) 
)

# Prepare insets
(ptsl.i1 <- ggplot()+
    # Base maps
    geom_sf(data = base,
            color = c("#CFCFCF"),
            size = 0.6,
            fill = "white") +
    # Get the raster
    geom_raster(data = sst.v, aes(x = x, y = y, fill = val)) +
    sca + # Add color scale
    geom_point(data = pts[[1]],
               aes(x = x, y = y), size = 1,
               alpha = .7, shape = 16)+
    # Establish area
    coord_sf(xlim = c(rects[[1]]$x1, rects[[1]]$x2),
             ylim = c(rects[[1]]$y1, rects[[1]]$y2),
             datum = st_crs(base),
             expand = F) + int)

(ptsl.i2 <- ggplot()+
    # Base maps
    geom_sf(data = base,
            color = c("#CFCFCF"),
            size = 0.6,
            fill = "white") +
    # Get the raster
    geom_raster(data = sst.v, aes(x = x, y = y, fill = val)) +
    sca + # Add color scale
    geom_point(data = pts[[1]],
               aes(x = x, y = y), size = 1,
               alpha = .7, shape = 16)+
    # Establish area
    coord_sf(xlim = c(rects[[2]]$x1, rects[[2]]$x2),
             ylim = c(rects[[2]]$y1, rects[[2]]$y2),
             datum = st_crs(base),
             expand = F) + int)

# Density plot of temperature
vals <- data.frame(sst = raster::extract(sst, pts[[1]]))
(dplot.ptsl <- ggplot(vals)+
    ggdist::stat_halfeye(aes(x = sst),
                         fill = "grey90",#"#5F7591",
                         interval_size_range = c(0.5,1),
                         color = "#192331"
                         )+
    geom_point(
      aes(x = sst, y = -0.1),
      ## draw horizontal lines instead of points
      shape = "|",
      size = 3,
      alpha = .2,
      color = "#192331"
    ) +
    scale_y_continuous(limits = c(-0.12, 1), expand = c(0,0))+
    scale_x_continuous(breaks = seq(20, 30, by = 2),
                       limits = c(21, 29), position = "top")+
    xlab("")+ylab("")+
    theme_dplot)

ptsl.f <- ptsl +
  annotate(
    "segment",
    x = c(rects[[1]]$x2, (rects[[2]]$x1+rects[[2]]$x2)/2),
    y = c((rects[[1]]$y1+rects[[1]]$y2)/2,
          rects[[2]]$y1),
    xend = c(1600, (rects[[2]]$x1+rects[[2]]$x2)/2),
    yend = c((rects[[1]]$y1+rects[[1]]$y2)/2,
             -4000),
    lineend = "round",
    colour = "grey20",
    size = 0.3
  ) +
  annotation_custom(
    grob = ggplotGrob(ptsl.i1),
    xmin = 900,
    xmax = 3700,
    ymin = 4000,
    ymax = 2000) +
  annotation_custom(
    grob = ggplotGrob(ptsl.i2),
    xmin = 2400,
    xmax = 4400,
    ymin = -4400,
    ymax = -3100) +
  annotation_custom(
    grob = ggplotGrob(dplot.ptsl),
    xmin = -3620,
    xmax = 0,
    ymin = -4615,
    ymax = -2100)

#### Eclu ----

rects <- list(data.frame(x1 = -2400, x2 = -1500, y1 = 2000, y2= 2700, label = "a"),
              data.frame(x1 = 2400, x2 = 3200, y1 = -2800, y2= -2200, label = "b"))

(ptse <- ggplot()+
    # Base maps
    geom_sf(data = base,
            color = c("#CFCFCF"),
            size = 0.6,
            fill = "white") +
    # Get the raster
    geom_raster(data = sst.v, aes(x = x, y = y, fill = val)) +
    sca + # Add color scale
    # Add occurrence points
    geom_point(data = pts[[2]],
               aes(x = x, y = y), size = 1,
               alpha = .5, shape = 16)+
    # Establish area
    coord_sf(xlim = c(-3239.984, 4510.636),
             ylim = c(-4614.575, 4596.965),
             datum = st_crs(base),
             expand = F,
             label_axes = list(
               top = "E",
               left = "",
               top = ""
             )) +
    # Draw rectangles
    geom_rect(data = rects[[1]],
              aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
              fill = NA, color = "grey20", size = .3)+
    # geom_text(data = rects[[1]],
    #            aes(x = (x2 - 100), y = (y2 - 100), label = label))+
    geom_rect(data = rects[[2]],
              aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
              fill = NA, color = "grey20", size = .3)+
    # geom_text(data = rects[[1]],
    #            aes(x = (x2 - 100), y = (y2 - 100), label = label))+
    # Add title
    geom_label(aes(label = "B", x = 4000, y = 4100),
               size = 8, fontface = "bold", color = "grey30",
               label.size = 0)+
    # Remove axis labels and add theme
    xlab(NULL) + ylab(NULL) + wlt +
    # Add grid
    scale_x_discrete(position = "top", breaks = c(-2000, 0, 2000, 4000)) 
)

# Prepare insets
(ptse.i1 <- ggplot()+
    # Base maps
    geom_sf(data = base,
            color = c("#CFCFCF"),
            size = 0.6,
            fill = "white") +
    # Get the raster
    geom_raster(data = sst.v, aes(x = x, y = y, fill = val)) +
    sca + # Add color scale
    geom_point(data = pts[[2]],
               aes(x = x, y = y), size = 1,
               alpha = .7, shape = 16)+
    # Establish area
    coord_sf(xlim = c(rects[[1]]$x1, rects[[1]]$x2),
             ylim = c(rects[[1]]$y1, rects[[1]]$y2),
             datum = st_crs(base),
             expand = F) + int)

(ptse.i2 <- ggplot()+
    # Base maps
    geom_sf(data = base,
            color = c("#CFCFCF"),
            size = 0.6,
            fill = "white") +
    # Get the raster
    geom_raster(data = sst.v, aes(x = x, y = y, fill = val)) +
    sca + # Add color scale
    geom_point(data = pts[[2]],
               aes(x = x, y = y), size = 1,
               alpha = .7, shape = 16)+
    # Establish area
    coord_sf(xlim = c(rects[[2]]$x1, rects[[2]]$x2),
             ylim = c(rects[[2]]$y1, rects[[2]]$y2),
             datum = st_crs(base),
             expand = F) + int)

# Density plot of temperature
vals <- data.frame(sst = raster::extract(sst, pts[[2]]))
(dplot.ptse <- ggplot(vals)+
    ggdist::stat_halfeye(aes(x = sst),
                         fill = "grey90",#"#5F7591",
                         interval_size_range = c(0.5,1),
                         color = "#192331"
    )+
    geom_point(
      aes(x = sst, y = -0.1),
      ## draw horizontal lines instead of points
      shape = "|",
      size = 3,
      alpha = .2,
      color = "#192331"
    ) +
    scale_y_continuous(limits = c(-0.12, 1), expand = c(0,0))+
    scale_x_continuous(breaks = seq(20, 30, by = 2),
                       limits = c(21, 29), position = "top")+
    xlab("")+ylab("")+
    theme_dplot)

ptse.f <- ptse +
  annotate(
    "segment",
    x = c(rects[[1]]$x2, (rects[[2]]$x1+rects[[2]]$x2)/2),
    y = c((rects[[1]]$y1+rects[[1]]$y2)/2,
          rects[[2]]$y1),
    xend = c(1600, (rects[[2]]$x1+rects[[2]]$x2)/2),
    yend = c((rects[[1]]$y1+rects[[1]]$y2)/2,
             -4000),
    lineend = "round",
    colour = "grey20",
    size = 0.3
  ) +
  annotation_custom(
    grob = ggplotGrob(ptse.i1),
    xmin = 900,
    xmax = 3700,
    ymin = 4000,
    ymax = 2000) +
  annotation_custom(
    grob = ggplotGrob(ptse.i2),
    xmin = 2400,
    xmax = 4400,
    ymin = -4400,
    ymax = -3100) +
  annotation_custom(
    grob = ggplotGrob(dplot.ptse),
    xmin = -3620,
    xmax = 0,
    ymin = -4615,
    ymax = -2100)



#### Trve ----
rects <- list(data.frame(x1 = -2400, x2 = -1500, y1 = 2000, y2= 2700, label = "a"),
              data.frame(x1 = 2400, x2 = 3200, y1 = -2800, y2= -2200, label = "b"))

(ptst <- ggplot()+
    # Base maps
    geom_sf(data = base,
            color = c("#CFCFCF"),
            size = 0.6,
            fill = "white") +
    # Get the raster
    geom_raster(data = sst.v, aes(x = x, y = y, fill = val)) +
    sca + # Add color scale
    # Add occurrence points
    geom_point(data = pts[[3]],
               aes(x = x, y = y), size = 1,
               alpha = .5, shape = 16)+
    # Establish area
    coord_sf(xlim = c(-3239.984, 4510.636),
             ylim = c(-4614.575, 4596.965),
             datum = st_crs(base),
             expand = F,
             label_axes = list(
               top = "E",
               left = "",
               top = ""
             )) +
    # Draw rectangles
    geom_rect(data = rects[[1]],
              aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
              fill = NA, color = "grey20", size = .3)+
    # geom_text(data = rects[[1]],
    #            aes(x = (x2 - 100), y = (y2 - 100), label = label))+
    geom_rect(data = rects[[2]],
              aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
              fill = NA, color = "grey20", size = .3)+
    # geom_text(data = rects[[1]],
    #            aes(x = (x2 - 100), y = (y2 - 100), label = label))+
    # Add title
    geom_label(aes(label = "C", x = 4000, y = 4100),
               size = 8, fontface = "bold", color = "grey30",
               label.size = 0)+
    # Remove axis labels and add theme
    xlab(NULL) + ylab(NULL) + wlt +
    # Add grid
    scale_x_discrete(position = "top", breaks = c(-2000, 0, 2000, 4000)) 
)

# Prepare insets
(ptst.i1 <- ggplot()+
    # Base maps
    geom_sf(data = base,
            color = c("#CFCFCF"),
            size = 0.6,
            fill = "white") +
    # Get the raster
    geom_raster(data = sst.v, aes(x = x, y = y, fill = val)) +
    sca + # Add color scale
    geom_point(data = pts[[3]],
               aes(x = x, y = y), size = 1,
               alpha = .7, shape = 16)+
    # Establish area
    coord_sf(xlim = c(rects[[1]]$x1, rects[[1]]$x2),
             ylim = c(rects[[1]]$y1, rects[[1]]$y2),
             datum = st_crs(base),
             expand = F) + int)

(ptst.i2 <- ggplot()+
    # Base maps
    geom_sf(data = base,
            color = c("#CFCFCF"),
            size = 0.6,
            fill = "white") +
    # Get the raster
    geom_raster(data = sst.v, aes(x = x, y = y, fill = val)) +
    sca + # Add color scale
    geom_point(data = pts[[3]],
               aes(x = x, y = y), size = 1,
               alpha = .7, shape = 16)+
    # Establish area
    coord_sf(xlim = c(rects[[2]]$x1, rects[[2]]$x2),
             ylim = c(rects[[2]]$y1, rects[[2]]$y2),
             datum = st_crs(base),
             expand = F) + int)

# Density plot of temperature
vals <- data.frame(sst = raster::extract(sst, pts[[3]]))
(dplot.ptst <- ggplot(vals)+
    ggdist::stat_halfeye(aes(x = sst),
                         fill = "grey90",#"#5F7591",
                         interval_size_range = c(0.5,1),
                         color = "#192331"
    )+
    geom_point(
      aes(x = sst, y = -0.1),
      ## draw horizontal lines instead of points
      shape = "|",
      size = 3,
      alpha = .2,
      color = "#192331"
    ) +
    scale_y_continuous(limits = c(-0.12, 1), expand = c(0,0))+
    scale_x_continuous(breaks = seq(20, 30, by = 2),
                       limits = c(21, 29), position = "top")+
    xlab("")+ylab("")+
    theme_dplot)

ptst.f <- ptst +
  annotate(
    "segment",
    x = c(rects[[1]]$x2, (rects[[2]]$x1+rects[[2]]$x2)/2),
    y = c((rects[[1]]$y1+rects[[1]]$y2)/2,
          rects[[2]]$y1),
    xend = c(1600, (rects[[2]]$x1+rects[[2]]$x2)/2),
    yend = c((rects[[1]]$y1+rects[[1]]$y2)/2,
             -4000),
    lineend = "round",
    colour = "grey20",
    size = 0.3
  ) +
  annotation_custom(
    grob = ggplotGrob(ptst.i1),
    xmin = 900,
    xmax = 3700,
    ymin = 4000,
    ymax = 2000) +
  annotation_custom(
    grob = ggplotGrob(ptst.i2),
    xmin = 2400,
    xmax = 4400,
    ymin = -4400,
    ymax = -3100) +
  annotation_custom(
    grob = ggplotGrob(dplot.ptst),
    xmin = -3620,
    xmax = 0,
    ymin = -4615,
    ymax = -2100)

# Save composite ----
final <- ptsl.f + ptse.f + ptst.f + plot_layout(nrow = 1, guides = "collect") &
  theme(legend.position='bottom')

ggsave(paste0("figures/species_pts_sst.jpg"), final,
       width = 38, height = 18, units = "cm")

### END