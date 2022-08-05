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
sp <- "trve"

# Load rasters generated before
curr <- raster(paste0("results/", sp, "/predictions/", sp, "_mean_current.tif"))
ssp1 <- raster(paste0("results/", sp, "/predictions/", sp, "_mean_ssp1.tif"))
ssp2 <- raster(paste0("results/", sp, "/predictions/", sp, "_mean_ssp2.tif"))
ssp3 <- raster(paste0("results/", sp, "/predictions/", sp, "_mean_ssp3.tif"))

# Get differences as percentage
ssp1.difs <- ssp1-curr
ssp2.difs <- ssp2-curr
ssp3.difs <- ssp3-curr

# Normalize all for easier interpretation

maxv <- max(maxValue(curr), maxValue(ssp1), maxValue(ssp2), maxValue(ssp3)) 
minv <- min(minValue(curr), minValue(ssp1), minValue(ssp2), minValue(ssp3)) 

curr2 <- round((curr - minValue(curr))/(maxValue(curr)-minValue(curr)),2)
ssp32 <- round((ssp3 - minValue(ssp3))/(maxValue(ssp3)-minValue(ssp3)),2)

th <- min(extract(curr2, pts))

curr <- (curr - minv)/(maxv - minv)
ssp1 <- (ssp1 - minv)/(maxv - minv)
ssp2 <- (ssp2 - minv)/(maxv - minv)
ssp3 <- (ssp3 - minv)/(maxv - minv)

# Convert to data.frame
get.val <- function(x){
        temp <- as(x, "SpatialPixelsDataFrame")
        temp <- as.data.frame(temp)
        colnames(temp) <- c("val", "x", "y")
        temp
}

curr.v <- get.val(curr)
ssp1.v <- get.val(ssp1)
ssp2.v <- get.val(ssp2)
ssp3.v <- get.val(ssp3)

# Get if reduct or increased
get.cg <- function(x){
        temp <- get.val(x)
        cg <- ifelse(temp$val < 0, "L",
                     ifelse(temp$val > 0, "I", "K"))
        cg
}

ssp1.cg <- get.cg(ssp1.difs)
ssp2.cg <- get.cg(ssp2.difs)
ssp3.cg <- get.cg(ssp3.difs)

ssp1.v$change <- ssp1.cg
ssp2.v$change <- ssp2.cg
ssp3.v$change <- ssp3.cg

# Prepare theme/ploting stuff ----
# Guide for legend
step.guide <- guide_coloursteps(title = "Intensity",
                                show.limits = TRUE,
                                barheight = unit(0.12, "in"),
                                barwidth = unit(3.5, "in"),
                                ticks = T,
                                ticks.colour = "grey20",
                                frame.colour = "grey20",
                                title.position = "top")

# step.guide.b <- guide_legend(title = "Change in intensity",
#                                 show.limits = TRUE,
#                                 barheight = unit(0.12, "in"),
#                                 barwidth = unit(3.5, "in"),
#                                 ticks = T,
#                                 ticks.colour = "grey20",
#                                 frame.colour = "grey20",
#                                 title.position = "top",
#                              label.position = "bottom")

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
sca <- scale_fill_stepsn(breaks = seq(0.1,0.9,length.out = 9),
                         limits = c(0, 1),
                         colors = RColorBrewer::
                                 brewer.pal(9, "GnBu"),
                         na.value = NA,
                         guide = step.guide)

scab <- scale_fill_stepsn(breaks = seq(0,1,length.out = 5),
                         limits = c(0, 1),
                         colors = c(RColorBrewer::
                                 brewer.pal(4, "GnBu")),
                         na.value = NA,
                         guide = step.guide)

sca.cg <- scale_fill_stepsn(breaks = seq(0.1,0.9,length.out = 9),
                         limits = c(0, 1),
                         colors = RColorBrewer::
                                 brewer.pal(9, "RdPu"),
                         na.value = NA,
                         guide = step.guide)

# sca.b <- scale_fill_stepsn(breaks = c(seq(-25, 25, by = 5)),
#                          limits = c(-25, 25),
#                          colors = RColorBrewer::
#                                  brewer.pal(10, "BrBG"),
#                          na.value = "white",
#                          guide = step.guide)


rects <- list(data.frame(x1 = -1500, x2 = 0, y1 = 2000, y2= 3000, label = "a"),
              data.frame(x1 = 2511, x2 = 4000, y1 = -1800, y2= -2808, label = "b"))


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
                   fill = NA, color = "grey20")+
         # geom_text(data = rects[[1]],
         #            aes(x = (x2 - 100), y = (y2 - 100), label = label))+
         geom_rect(data = rects[[2]],
                   aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
                   fill = NA, color = "grey20")+
         # geom_text(data = rects[[1]],
         #            aes(x = (x2 - 100), y = (y2 - 100), label = label))+
         # Add title
         geom_label(aes(label = "A", x = 4000, y = 4100),
                    size = 10, fontface = "bold", color = "grey30",
                    label.size = 0)+
         # Remove axis labels and add theme
         xlab("Easting") + ylab("Northing") + wlt +
         # Add grid
         scale_x_discrete(position = "top") 
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
                size = 0.5
        ) +
        annotation_custom(
                grob = ggplotGrob(pc.i1),
                xmin = 1500,
                xmax = 4000,
                ymin = 3000,
                ymax = 1400) +
        annotation_custom(
                grob = ggplotGrob(pc.i2),
                xmin = -2500,
                xmax = 1000,
                ymin = -3100,
                ymax = -900)


ggsave(paste0("figures/", sp, "_lgcp_current.jpg"), pcf,
       width = 16, height = 18, units = "cm")







###################

ssp1.pol <- ssp1.difs
ssp1.pol <- calc(ssp1.pol, function(x){
        if (is.na(x)) {
                NA
        }else{
                x <- ifelse(x < 0, -1,
                            ifelse(x > 0, 1, 0))
                x
        }
})

pols <- function(x){
        st <- list()
        for (i in 1:2) {
                val <- c(-1, 1)[i]
                temp <- calc(x, function(z){
                        if (is.na(z)) {
                                NA
                        }else{
                                ifelse(z == val, 1, NA)
                        }
                })
                st[[i]] <- rasterToPolygons(temp, dissolve = T)
        }
        names(st) <- c("loss", "gain")
        st
}

ssp1.pol.t <- pols(ssp1.pol)
ssp1.pol.t <- lapply(ssp1.pol.t, st_as_sf)


# SSP1
scab <- scale_fill_stepsn(breaks = seq(0,1,length.out = 5),
                          limits = c(0, 1),
                          colors = c("#DAE6E4", "#47909E", "#1C6C94", "#034F82"),
                          na.value = NA,
                          guide = step.guide)

(p1 <- ggplot()+
                # Base maps
                geom_sf(data = base,
                        color = c("#CFCFCF"),
                        size = 0.6,
                        fill = "white") +
                # Get the raster
                geom_raster(data = curr.v,
                            aes(x = x, y = y, fill = val)) +
                sca+
                #geom_sf(data = ssp1.pol.t$loss, color = "red", size = .3, fill = NA)+
                #geom_sf(data = ssp1.pol.t$gain, color = "blue", size = .3, fill = NA)+
                # geom_raster(data = ssp1.v[ssp1.v$change == "K" | ssp1.v$change == "I",],
                #             aes(x = x, y = y, fill = val)) +
                # sca+
                # new_scale_fill()+
                # geom_raster(data = ssp1.v[ssp1.v$change == "K" | ssp1.v$change == "L",],
                #             aes(x = x, y = y, fill = val)) +
                # sca.b+
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
                xlab(NULL) + ylab(NULL) + wlt 
 # +
 #                theme(plot.margin = margin(1, 1, 0, 0, "pt"),
 #                      plot.background = 
 #                              element_rect(color = "grey30", fill = 'white'),
 #                      axis.text = element_blank(),
 #                      panel.grid.major = element_blank()
 #                      )
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
                sca.b + # Add color scale
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
                xlab(NULL) + ylab(NULL) + wlt #+
                # theme(plot.margin = margin(1, 1, 0, 0, "pt"),
                #       plot.background = 
                #               element_rect(color = "grey30", fill = 'white'),
                #       axis.text = element_blank(),
                #       panel.grid.major = element_blank()
                # )
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
                sca.b + # Add color scale
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
                xlab(NULL) + ylab(NULL) + wlt #+
                # theme(plot.margin = margin(1, 1, 0, 0, "pt"),
                #       plot.background = 
                #               element_rect(color = "grey30", fill = 'white'),
                #       axis.text = element_blank(),
                #       panel.grid.major = element_blank()
                # )
)

# Get final plot ----
(fp <- p1 + p2 + p3)

ggsave(paste0("figures/", sp, "_sst_composite.jpg"), fp,
       width = 48, height = 18, units = "cm")

ggsave(paste0("figures/", sp, "_sst_mainmap.jpg"), pc,
       width = 21, height = 18, units = "cm")
