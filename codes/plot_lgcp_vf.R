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

# Load results ----
# Define species
sp <- "lyva"

# Load species occurrence data
occ <- SpatialPoints(read.csv(paste0("data/", sp, "/", sp, "_cell.csv"))[,1:2],
                     proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
occ <- spTransform(occ, crs(base))


# Define model
m <- 4

# Load rasters generated before
curr <- raster(paste0("results/", sp, "/predictions/", sp, "_mean_m", m, "_current.tif"))
ssp1 <- raster(paste0("results/", sp, "/predictions/", sp, "_mean_m", m, "_ssp1.tif"))
ssp2 <- raster(paste0("results/", sp, "/predictions/", sp, "_mean_m", m, "_ssp2.tif"))
ssp3 <- raster(paste0("results/", sp, "/predictions/", sp, "_mean_m", m, "_ssp3.tif"))

# Normalize for easier handling
miv <- min(c(minValue(curr), minValue(ssp1), minValue(ssp2), minValue(ssp3)))
mav <- max(c(maxValue(curr), maxValue(ssp1), maxValue(ssp2), maxValue(ssp3)))

curr <- (curr - miv)/(mav - miv)
ssp1 <- (ssp1 - miv)/(mav - miv)
ssp2 <- (ssp2 - miv)/(mav - miv)
ssp3 <- (ssp3 - miv)/(mav - miv)

# Convert to data.frame
get.val <- function(x, mode = "normal"){
        temp <- as(x, "SpatialPixelsDataFrame")
        temp <- as.data.frame(temp)
        colnames(temp) <- c("val", "x", "y")
        return(temp)
}

curr.v <- get.val(curr)
ssp1.v <- get.val(ssp1)
ssp2.v <- get.val(ssp2)
ssp3.v <- get.val(ssp3)



# Prepare theme/ploting stuff ----
# Guide for legend
step.guide <- guide_colorbar(title = "Intensity",
                             show.limits = TRUE,
                             barheight = unit(0.12, "in"),
                             barwidth = unit(3.5, "in"),
                             ticks = F,
                             ticks.colour = "grey20",
                             frame.colour = "grey20",
                             title.position = "top")


# Themes
nlt <- theme_classic()+
        theme(axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.title.x.top = element_text(vjust = 2),
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
              axis.text = element_text(colour = "grey60"),
              axis.title.x.top = element_text(vjust = 2),
              panel.background = element_blank(),
              panel.border = element_rect(fill = NA, color = "grey60"),
              panel.grid.major = element_line(linetype = 'dashed', 
                                              colour = "grey70",
                                              size = .1),
              legend.position="bottom",
              legend.title.align=0.5,
              legend.background = element_rect(fill = "white")
        )

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


# Scale of values
sca <- scale_fill_gradientn(
        colours = rev(c("#A84C00", "#D97D27", "#F5BD44", "#FFD561", "#FFF291",
                        "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695")),
        limits = c(0,1),
        guide = step.guide,
        labels = c("Low", "", "", "", "High"))



# Generate plots ----
##### Current ----
# Rectangles to highlight
rects <- list(data.frame(x1 = -1600, x2 = -100, y1 = 2300, y2= 3300, label = "a"),
              data.frame(x1 = 2511, x2 = 4000, y1 = -1800, y2= -2808, label = "b"))

(pc <- ggplot()+
         # Base maps
         geom_sf(data = base,
                 color = c("#CFCFCF"),
                 size = 0.6,
                 fill = "white") +
         # Get the raster
         geom_raster(data = curr.v, aes(x = x, y = y, fill = val)) +
         sca + # Add color scale
         # Add occurrence points
         # geom_point(data = data.frame(occ@coords),
         #            aes(x = decimalLongitude, y = decimalLatitude), size = 1,
         #            alpha = .2, shape = 16)+
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
                    size = 10, fontface = "bold", color = "grey30",
                    label.size = 0)+
         # Remove axis labels and add theme
         xlab("Easting") + ylab("Northing") + wlt +
         # Add grid
         scale_x_discrete(position = "top", breaks = c(-2000, 0, 2000, 4000)) 
)

# Prepare insets
(pc.i1 <- ggplot()+
                # Base maps
                geom_sf(data = base,
                        color = c("#CFCFCF"),
                        size = 0.6,
                        fill = "white") +
                # Get the raster
                geom_raster(data = curr.v, aes(x = x, y = y, fill = val)) +
                sca + # Add color scale
                # Establish area
                coord_sf(xlim = c(rects[[1]]$x1, rects[[1]]$x2),
                         ylim = c(rects[[1]]$y1, rects[[1]]$y2),
                         datum = st_crs(base),
                         expand = F) + int)

(pc.i2 <- ggplot()+
                # Base maps
                geom_sf(data = base,
                        color = c("#CFCFCF"),
                        size = 0.6,
                        fill = "white") +
                # Get the raster
                geom_raster(data = curr.v, aes(x = x, y = y, fill = val)) +
                sca + # Add color scale
                # Establish area
                coord_sf(xlim = c(rects[[2]]$x1, rects[[2]]$x2),
                         ylim = c(rects[[2]]$y1, rects[[2]]$y2),
                         datum = st_crs(base),
                         expand = F) + int)

pcf <- pc +
        annotate(
                "segment",
                x = c(rects[[1]]$x2, 500),
                y = c((rects[[1]]$y1+rects[[1]]$y2)/2,
                      (rects[[2]]$y1+rects[[2]]$y2)/2),
                xend = c(1600, rects[[2]]$x1),
                yend = c((rects[[1]]$y1+rects[[1]]$y2)/2,
                         (rects[[2]]$y1+rects[[2]]$y2)/2),
                lineend = "round",
                colour = "grey20",
                size = 0.3
        ) +
        annotation_custom(
                grob = ggplotGrob(pc.i1),
                xmin = 1400,
                xmax = 4100,
                ymin = 3200,
                ymax = 1200) +
        annotation_custom(
                grob = ggplotGrob(pc.i2),
                xmin = -2500,
                xmax = 1000,
                ymin = -3100,
                ymax = -900)


ggsave(paste0("figures/", sp, "_lgcp_current_m", m,".jpg"), pcf,
       width = 16, height = 18, units = "cm")



##### SSP1 ----
(ps1 <- ggplot()+
                # Base maps
                geom_sf(data = base,
                        color = c("#CFCFCF"),
                        size = 0.6,
                        fill = "white") +
                # Get the raster
                geom_raster(data = ssp1.v, aes(x = x, y = y, fill = val)) +
                sca + # Add color scale
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
                           size = 10, fontface = "bold", color = "grey30",
                           label.size = 0)+
                # Remove axis labels and add theme
                xlab(NULL) + ylab(NULL) + wlt +
                # Add grid
                scale_x_discrete(position = "top", breaks = c(-2000, 0, 2000, 4000))
)

# Prepare insets
(ps1.i1 <- ggplot()+
                # Base maps
                geom_sf(data = base,
                        color = c("#CFCFCF"),
                        size = 0.6,
                        fill = "white") +
                # Get the raster
                geom_raster(data = ssp1.v, aes(x = x, y = y, fill = val)) +
                sca + # Add color scale
                # Establish area
                coord_sf(xlim = c(rects[[1]]$x1, rects[[1]]$x2),
                         ylim = c(rects[[1]]$y1, rects[[1]]$y2),
                         datum = st_crs(base),
                         expand = F) + int)

(ps1.i2 <- ggplot()+
                # Base maps
                geom_sf(data = base,
                        color = c("#CFCFCF"),
                        size = 0.6,
                        fill = "white") +
                # Get the raster
                geom_raster(data = ssp1.v, aes(x = x, y = y, fill = val)) +
                sca + # Add color scale
                # Establish area
                coord_sf(xlim = c(rects[[2]]$x1, rects[[2]]$x2),
                         ylim = c(rects[[2]]$y1, rects[[2]]$y2),
                         datum = st_crs(base),
                         expand = F) + int)

ps1f <- ps1 +
        annotate(
                "segment",
                x = c(rects[[1]]$x2, 500),
                y = c((rects[[1]]$y1+rects[[1]]$y2)/2,
                      (rects[[2]]$y1+rects[[2]]$y2)/2),
                xend = c(1600, rects[[2]]$x1),
                yend = c((rects[[1]]$y1+rects[[1]]$y2)/2,
                         (rects[[2]]$y1+rects[[2]]$y2)/2),
                lineend = "round",
                colour = "grey20",
                size = 0.3
        ) +
        annotation_custom(
                grob = ggplotGrob(ps1.i1),
                xmin = 1400,
                xmax = 4100,
                ymin = 3200,
                ymax = 1200) +
        annotation_custom(
                grob = ggplotGrob(ps1.i2),
                xmin = -2500,
                xmax = 1000,
                ymin = -3100,
                ymax = -900)

ps1f.s <- ps1f + theme(legend.position = "none")

ggsave(paste0("figures/", sp, "_lgcp_ssp1_m", m,".jpg"), ps1f.s,
       width = 16, height = 18, units = "cm")




##### SSP2 ----
(ps2 <- ggplot()+
                # Base maps
                geom_sf(data = base,
                        color = c("#CFCFCF"),
                        size = 0.6,
                        fill = "white") +
                # Get the raster
                geom_raster(data = ssp2.v, aes(x = x, y = y, fill = val)) +
                sca + # Add color scale
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
                           size = 10, fontface = "bold", color = "grey30",
                           label.size = 0)+
                # Remove axis labels and add theme
                xlab(NULL) + ylab(NULL) + wlt +
                # Add grid
                scale_x_discrete(position = "top", breaks = c(-2000, 0, 2000, 4000)) 
)

# Prepare insets
(ps2.i1 <- ggplot()+
                # Base maps
                geom_sf(data = base,
                        color = c("#CFCFCF"),
                        size = 0.6,
                        fill = "white") +
                # Get the raster
                geom_raster(data = ssp2.v, aes(x = x, y = y, fill = val)) +
                sca + # Add color scale
                # Establish area
                coord_sf(xlim = c(rects[[1]]$x1, rects[[1]]$x2),
                         ylim = c(rects[[1]]$y1, rects[[1]]$y2),
                         datum = st_crs(base),
                         expand = F) + int)

(ps2.i2 <- ggplot()+
                # Base maps
                geom_sf(data = base,
                        color = c("#CFCFCF"),
                        size = 0.6,
                        fill = "white") +
                # Get the raster
                geom_raster(data = ssp2.v, aes(x = x, y = y, fill = val)) +
                sca + # Add color scale
                # Establish area
                coord_sf(xlim = c(rects[[2]]$x1, rects[[2]]$x2),
                         ylim = c(rects[[2]]$y1, rects[[2]]$y2),
                         datum = st_crs(base),
                         expand = F) + int)

ps2f <- ps2 +
        annotate(
                "segment",
                x = c(rects[[1]]$x2, 500),
                y = c((rects[[1]]$y1+rects[[1]]$y2)/2,
                      (rects[[2]]$y1+rects[[2]]$y2)/2),
                xend = c(1600, rects[[2]]$x1),
                yend = c((rects[[1]]$y1+rects[[1]]$y2)/2,
                         (rects[[2]]$y1+rects[[2]]$y2)/2),
                lineend = "round",
                colour = "grey20",
                size = 0.3
        ) +
        annotation_custom(
                grob = ggplotGrob(ps2.i1),
                xmin = 1400,
                xmax = 4100,
                ymin = 3200,
                ymax = 1200) +
        annotation_custom(
                grob = ggplotGrob(ps2.i2),
                xmin = -2500,
                xmax = 1000,
                ymin = -3100,
                ymax = -900)

ps2f.s <- ps2f + theme(legend.position = "none")

ggsave(paste0("figures/", sp, "_lgcp_ssp2_m", m,".jpg"), ps2f.s,
       width = 16, height = 18, units = "cm")




##### SSP3 ----
(ps3 <- ggplot()+
                # Base maps
                geom_sf(data = base,
                        color = c("#CFCFCF"),
                        size = 0.6,
                        fill = "white") +
                # Get the raster
                geom_raster(data = ssp3.v, aes(x = x, y = y, fill = val)) +
                sca + # Add color scale
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
                geom_label(aes(label = "D", x = 4000, y = 4100),
                           size = 10, fontface = "bold", color = "grey30",
                           label.size = 0)+
                # Remove axis labels and add theme
                xlab(NULL) + ylab(NULL) + wlt +
                # Add grid
                scale_x_discrete(position = "top", breaks = c(-2000, 0, 2000, 4000)) 
)

# Prepare insets
(ps3.i1 <- ggplot()+
                # Base maps
                geom_sf(data = base,
                        color = c("#CFCFCF"),
                        size = 0.6,
                        fill = "white") +
                # Get the raster
                geom_raster(data = ssp3.v, aes(x = x, y = y, fill = val)) +
                sca + # Add color scale
                # Establish area
                coord_sf(xlim = c(rects[[1]]$x1, rects[[1]]$x2),
                         ylim = c(rects[[1]]$y1, rects[[1]]$y2),
                         datum = st_crs(base),
                         expand = F) + int)

(ps3.i2 <- ggplot()+
                # Base maps
                geom_sf(data = base,
                        color = c("#CFCFCF"),
                        size = 0.6,
                        fill = "white") +
                # Get the raster
                geom_raster(data = ssp3.v, aes(x = x, y = y, fill = val)) +
                sca + # Add color scale
                # Establish area
                coord_sf(xlim = c(rects[[2]]$x1, rects[[2]]$x2),
                         ylim = c(rects[[2]]$y1, rects[[2]]$y2),
                         datum = st_crs(base),
                         expand = F) + int)

ps3f <- ps3 +
        annotate(
                "segment",
                x = c(rects[[1]]$x2, 500),
                y = c((rects[[1]]$y1+rects[[1]]$y2)/2,
                      (rects[[2]]$y1+rects[[2]]$y2)/2),
                xend = c(1600, rects[[2]]$x1),
                yend = c((rects[[1]]$y1+rects[[1]]$y2)/2,
                         (rects[[2]]$y1+rects[[2]]$y2)/2),
                lineend = "round",
                colour = "grey20",
                size = 0.3
        ) +
        annotation_custom(
                grob = ggplotGrob(ps3.i1),
                xmin = 1400,
                xmax = 4100,
                ymin = 3200,
                ymax = 1200) +
        annotation_custom(
                grob = ggplotGrob(ps3.i2),
                xmin = -2500,
                xmax = 1000,
                ymin = -3100,
                ymax = -900)

ps3f.s <- ps3f + theme(legend.position = "none")

ggsave(paste0("figures/", sp, "_lgcp_ssp3_m", m,".jpg"), ps3f.s,
       width = 16, height = 18, units = "cm")

# Save composite ----
final <- pcf + ps1f + ps2f + ps3f + plot_layout(nrow = 1, guides = "collect") &
        theme(legend.position='bottom')

ggsave(paste0("figures/", sp, "_lgcp_comp_m", m,".jpg"), final,
       width = 50, height = 18, units = "cm")
