library(ggplot2)
library(sf)
library(raster)
library(patchwork)

curr <- raster(paste0("results/lyva/predictions/lyva_mean_current.tif"))
fut <- raster(paste0("results/lyva/predictions/lyva_mean_ssp1.tif"))

dif.f <- (fut-curr) / curr

curr <- (curr - minValue(curr))/(maxValue(curr) - minValue(curr))

curr <- as(curr, "SpatialPixelsDataFrame")

curr <- as.data.frame(curr)
colnames(curr) <- c("val", "x", "y")

dif.f <- as(dif.f, "SpatialPixelsDataFrame")

dif.f <- as.data.frame(dif.f)
colnames(dif.f) <- c("val", "x", "y")

retas <- data.frame(x = c(1256.5727, 1256.5727, 4890.4723, 4890.4723),
                    xend = c(1256.5727, 4890.4723,4890.4723, 2325.4219),
                    y = c(4470.2438, 1996.0985, -373.9687, -993.9687),
                    yend = c(1996.0985, -373.9687,-993.9687, -4848.4605))

#curr <- st_as_sf(curr)

# curr <- data.frame(rasterToPoints(curr))
# colnames(curr) <- c("x", "y", "val")
#curr$val <- curr$val * 1e6

# America
base1 <- shapefile("gis/basemaps/ne_110m_land_edited.shp")
base1 <- buffer(base1, 0)
base1 <- crop(base1, extent(-120, -10, -55, 65))
proj <- "+proj=laea +lat_0=0 +lon_0=-70 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs"
base1 <- spTransform(base1, CRS(proj))

# Convert to SF
base1 <- st_as_sf(base1)

step.guide <- guide_coloursteps(title = "Intensity",
                                show.limits = TRUE,
                                barheight = unit(0.12, "in"),
                                barwidth = unit(3.5, "in"),
                                ticks = T,
                                ticks.colour = "grey20",
                                frame.colour = "grey20",
                                #title.hjust = unit(5, "in"),
                                title.position = "top")

# Maps scales
main.s.x <-
        scale_x_continuous(breaks = seq(-35,-100,-10),
                           limits = c(-100,-29))
main.s.y <-
        scale_y_continuous(breaks = seq(-40, 40, 10),
                           limits = c(-42.5, 42.5))

rects <- list(data.frame(x1 = -85, x2 = -65, y1 = 20, y2= 30, label = "A"),
              data.frame(x1 = -45, x2 = -33, y1 = -25, y2= 0, label = "B"))

# Generate current map
(pdif <- ggplot()+
                # Base maps
                geom_sf(data = base1,
                        color = NA,
                        size = 0.6,
                        fill = NA) +
                # Get the raster
                geom_tile(data = dif.f, aes(x = x, y = y, fill = val))+
                scale_fill_stepsn(breaks = c(-10, -5, -2.5, -1, 0,
                                             1, 2.5, 5, 10),
                                  limits = c(-20, 20),
                                  values = scales::rescale(c(-20,-10, -5, -2.5, -1, 0,
                                             1, 2.5, 5, 10, 20)),
                                  colors = RColorBrewer::
                                          brewer.pal(11, "BrBG"),
                                  na.value = NA,
                                  guide = step.guide)+
                # Establish area
                coord_sf(xlim = c(-3239.984, 5200.636),
                         ylim = c(-4614.575, 4596.965),
                         datum = st_crs(base1),
                         expand = F,
                         label_axes = list(
                                 top = "E",
                                 left = "N",
                                 top = ""
                         ))+
                geom_segment(data = retas,
                             aes(x = x, xend = xend, y = y, yend = yend),
                             size = 2, color = "grey90",
                             lineend = "round")+
                # coord_sf(
                #         xlim = c(-3239.984, 4510.636),
                #         ylim = c(-4905.32, 4906.17),
                #         expand = FALSE,
                #         crs = CRS(proj),
                #         label_axes = list(
                #                 top = "E",
                #                 left = "N",
                #                 top = ""
                #         )
                # ) +
        # draw three rectangles
        # geom_rect(data = rects[[1]],
        #           aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
        #           fill = NA, color = "grey20")+
        # geom_text(data = rects[[1]],
        #            aes(x = (x2 - 2), y = (y2 - 2), label = label))+
        # draw three rectangles
        # geom_rect(data = rects[[2]],
        #           aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
        #           fill = NA, color = "grey20")+
        # geom_text(data = rects[[2]],
        #           aes(x = (x2 - 2), y = (y2 - 2), label = label))+
        geom_text(aes(label = "SSP1", x = -2700, y = 4400),
                  size = 10, fontface = "bold", color = "grey30")+
        # Add plot details
        xlab(NULL) + ylab(NULL)+
                #main.s.x + main.s.y+
                theme_classic()+
                theme(axis.line = element_blank(),
                      #axis.text = element_blank(),
                      axis.text = element_blank(),
                      axis.ticks = element_blank(),
                      #panel.border = element_rect(colour = "grey70", fill = NA),
                      legend.position="none",
                      legend.title.align=0.5,
                      # panel.grid.major = element_line(linetype = 'dashed', 
                      #                                 colour = "grey90",
                      #                                 size = .05)
                      panel.background = element_blank(),
                      plot.background = element_blank()
                      ))

(pdif2 <- ggplot()+
                # Base maps
                geom_sf(data = base1,
                        color = NA,
                        size = 0.6,
                        fill = NA) +
                # Get the raster
                geom_tile(data = dif.f, aes(x = x, y = y, fill = val))+
                scale_fill_stepsn(breaks = c(-10, -5, -2.5, -1, 0,
                                             1, 2.5, 5, 10),
                                  limits = c(-20, 20),
                                  values = scales::rescale(c(-20,-10, -5, -2.5, -1, 0,
                                                             1, 2.5, 5, 10, 20)),
                                  colors = RColorBrewer::
                                          brewer.pal(11, "BrBG"),
                                  na.value = NA,
                                  guide = step.guide)+
                # Establish area
                coord_sf(xlim = c(-3239.984, 5200.636),
                         ylim = c(-4614.575, 4596.965),
                         datum = st_crs(base1),
                         expand = F,
                         label_axes = list(
                                 top = "E",
                                 left = "N",
                                 top = ""
                         ))+
                geom_segment(data = retas,
                             aes(x = x, xend = xend, y = y, yend = yend),
                             size = 2, color = "grey90",
                             lineend = "round")+
                # coord_sf(
                #         xlim = c(-3239.984, 4510.636),
                #         ylim = c(-4905.32, 4906.17),
                #         expand = FALSE,
                #         crs = CRS(proj),
                #         label_axes = list(
                #                 top = "E",
                #                 left = "N",
                #                 top = ""
                #         )
                # ) +
        # draw three rectangles
        # geom_rect(data = rects[[1]],
        #           aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
        #           fill = NA, color = "grey20")+
        # geom_text(data = rects[[1]],
        #            aes(x = (x2 - 2), y = (y2 - 2), label = label))+
        # draw three rectangles
        # geom_rect(data = rects[[2]],
        #           aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
        #           fill = NA, color = "grey20")+
        # geom_text(data = rects[[2]],
        #           aes(x = (x2 - 2), y = (y2 - 2), label = label))+
        geom_text(aes(label = "SSP1", x = -2700, y = 4400),
                  size = 10, fontface = "bold", color = "grey30")+
        # Add plot details
        xlab(NULL) + ylab(NULL)+
                #main.s.x + main.s.y+
                theme_classic()+
                theme(axis.line = element_blank(),
                      #axis.text = element_blank(),
                      axis.text = element_blank(),
                      axis.ticks = element_blank(),
                      #panel.border = element_rect(colour = "grey70", fill = NA),
                      legend.position="bottom",
                      legend.title.align=0.5,
                      # panel.grid.major = element_line(linetype = 'dashed', 
                      #                                 colour = "grey90",
                      #                                 size = .05)
                      panel.background = element_blank(),
                      plot.background = element_blank()
                ))



(pdif3 <- ggplot()+
                # Base maps
                geom_sf(data = base1,
                        color = NA,
                        size = 0.6,
                        fill = NA) +
                # Get the raster
                geom_tile(data = dif.f, aes(x = x, y = y, fill = val))+
                scale_fill_stepsn(breaks = c(-10, -5, -2.5, -1, 0,
                                             1, 2.5, 5, 10),
                                  limits = c(-20, 20),
                                  values = scales::rescale(c(-20,-10, -5, -2.5, -1, 0,
                                                             1, 2.5, 5, 10, 20)),
                                  colors = RColorBrewer::
                                          brewer.pal(11, "BrBG"),
                                  na.value = NA,
                                  guide = step.guide)+
                # Establish area
                coord_sf(xlim = c(-3239.984, 5200.636),
                         ylim = c(-4614.575, 4596.965),
                         datum = st_crs(base1),
                         expand = F,
                         label_axes = list(
                                 top = "E",
                                 left = "N",
                                 top = ""
                         ))+
                # geom_segment(data = retas,
                #              aes(x = x, xend = xend, y = y, yend = yend),
                #              size = 2, color = "grey90",
                #              lineend = "round")+
                # coord_sf(
                #         xlim = c(-3239.984, 4510.636),
                #         ylim = c(-4905.32, 4906.17),
                #         expand = FALSE,
                #         crs = CRS(proj),
                #         label_axes = list(
                #                 top = "E",
                #                 left = "N",
                #                 top = ""
                #         )
                # ) +
        # draw three rectangles
        # geom_rect(data = rects[[1]],
        #           aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
        #           fill = NA, color = "grey20")+
        # geom_text(data = rects[[1]],
        #            aes(x = (x2 - 2), y = (y2 - 2), label = label))+
        # draw three rectangles
        # geom_rect(data = rects[[2]],
        #           aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
        #           fill = NA, color = "grey20")+
        # geom_text(data = rects[[2]],
        #           aes(x = (x2 - 2), y = (y2 - 2), label = label))+
        
        # Add plot details
        geom_text(aes(label = "SSP1", x = -2700, y = 4400),
                  size = 10, fontface = "bold", color = "grey30")+
        xlab(NULL) + ylab(NULL)+
                #main.s.x + main.s.y+
                theme_classic()+
                theme(axis.line = element_blank(),
                      #axis.text = element_blank(),
                      axis.text = element_blank(),
                      axis.ticks = element_blank(),
                      #panel.border = element_rect(colour = "grey70", fill = NA),
                      legend.position="none",
                      legend.title.align=0.5,
                      # panel.grid.major = element_line(linetype = 'dashed', 
                      #                                 colour = "grey90",
                      #                                 size = .05)
                      panel.background = element_blank(),
                      plot.background = element_blank()
                ))

pdif2



(p2 <- ggplot()+
                # Base maps
                geom_sf(data = base1,
                        color = c("#DEDEDE"),
                        size = 0.6,
                        fill = "white") +
                # Get the raster
                geom_raster(data = curr, aes(x = x, y = y, fill = val))+
                scale_fill_stepsn(breaks = seq(10,90,10),
                                  limits = c(0, 100),
                                  colors = RColorBrewer::
                                          brewer.pal(10, "Spectral"),
                                  na.value = NA,
                                  guide = step.guide)+
                # Establish area
                coord_sf(
                        xlim = c(rects[[1]]$x1, rects[[1]]$x2),
                        ylim = c(rects[[1]]$y1, rects[[1]]$y2),
                        expand = FALSE#,
                        # label_axes = list(
                        #         top = "E",
                        #         left = "N",
                        #         top = ""
                        # )
                ) +
                # add text
                geom_text(data = rects[[1]],
                          aes(x = (x2 - 1), y = (y2 - 1), label = label),
                          size = 15)+
                # Add plot details
                xlab(NULL) + ylab(NULL)+
                main.s.x + main.s.y+
                theme_classic()+
                theme(axis.line = element_blank(),
                      axis.text = element_blank(),
                      #axis.text = element_text(colour = "grey60"),
                      axis.ticks = element_blank(),
                      panel.border = element_rect(colour = "grey70", fill = NA),
                      legend.position="none",
                      #legend.title.align=0.5,
                      # panel.grid.major = element_line(linetype = 'dashed', 
                      #                                 colour = "grey90",
                      #                                 size = .05)
                      )
        )


(p3 <- ggplot()+
                # Base maps
                geom_sf(data = base1,
                        color = c("#DEDEDE"),
                        size = 0.6,
                        fill = "white") +
                # Get the raster
                geom_raster(data = curr, aes(x = x, y = y, fill = val))+
                scale_fill_stepsn(breaks = seq(10,90,10),
                                  limits = c(0, 100),
                                  colors = RColorBrewer::
                                          brewer.pal(10, "Spectral"),
                                  na.value = NA,
                                  guide = step.guide)+
                # Establish area
                coord_sf(
                        xlim = c(rects[[2]]$x1-7, rects[[2]]$x2),
                        ylim = c(rects[[2]]$y1, rects[[2]]$y2),
                        expand = FALSE#,
                        # label_axes = list(
                        #         top = "E",
                        #         left = "N",
                        #         top = ""
                        # )
                ) +
                # add text
                geom_text(data = rects[[2]],
                          aes(x = (x2 - 1), y = (y2 - 1), label = label),
                          size = 15)+
                # Add plot details
                xlab(NULL) + ylab(NULL)+
                main.s.x + main.s.y+
                theme_classic()+
                theme(axis.line = element_blank(),
                      axis.text = element_blank(),
                      #axis.text = element_text(colour = "grey60"),
                      axis.ticks = element_blank(),
                      panel.border = element_rect(colour = "grey70", fill = NA),
                      legend.position="none",
                      #legend.title.align=0.5,
                      # panel.grid.major = element_line(linetype = 'dashed', 
                      #                                 colour = "grey90",
                      #                                 size = .05)
                )
)

library(patchwork)

comp1 <- p2 / p3
comp1

comp2 <- p1 + comp1
comp2





layout <- c(
        patchwork::area(t = 1.1, l = 1  , b = 5.5, r = 4),
        patchwork::area(t = 1  , l = 3  , b = 5  , r = 6),
        patchwork::area(t = 1.1, l = 5  , b = 5.5  , r = 8)
)
pdif + pdif2 + pdif3 +
        plot_layout(design = layout)
