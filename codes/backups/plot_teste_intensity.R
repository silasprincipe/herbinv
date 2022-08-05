library(ggplot2)
library(sf)
library(raster)

curr <- raster(paste0("results/lyva/predictions/lyva_mean_current.tif"))
curr <- (curr - minValue(curr))/(maxValue(curr) - minValue(curr))
curr <- data.frame(rasterToPoints(curr))
colnames(curr) <- c("x", "y", "val")
#curr$val <- curr$val * 1e6

# America
base1 <- shapefile("gis/basemaps/ne_110m_land_edited.shp")
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
(p1 <- ggplot()+
        # Base maps
        geom_sf(data = base1,
                color = c("#DEDEDE"),
                size = 0.6,
                fill = "white") +
        # Get the raster
        geom_raster(data = curr, aes(x = x, y = y, fill = val))+
        scale_fill_stepsn(breaks = seq(10,round(max(curr$val)),length.out = 9),
                          limits = c(0, round(max(curr$val))),
                          colors = RColorBrewer::
                                  brewer.pal(9, "Blues"),
                          na.value = NA,
                          guide = step.guide)+
        # Establish area
        coord_sf(
                xlim = c(-3239.984, 4510.636),
                ylim = c(-4905.32, 4906.17),
                expand = FALSE,
                label_axes = list(
                        top = "E",
                        left = "N",
                        top = ""
                )
        ) +
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
        xlab(NULL) + ylab(NULL)+
        #main.s.x + main.s.y+
        theme_classic()+
        theme(axis.line = element_blank(),
              #axis.text = element_blank(),
              axis.text = element_text(colour = "grey60"),
              axis.ticks = element_blank(),
              #panel.border = element_rect(colour = "grey70", fill = NA),
              legend.position="bottom",
              legend.title.align=0.5,
              panel.grid.major = element_line(linetype = 'dashed', 
                                              colour = "grey90",
                                              size = .05)))


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
