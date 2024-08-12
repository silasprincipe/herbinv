#### Modelling of coral reef herbivorous invertebrates ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2022 ##

# Plot of LGCP results - all species [supplementary material]

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

# Define species
species <- c("lyva", "eclu", "trve")

for (sp in species) {
 
  # Load results ----
  
  # Load species occurrence data
  occ <- SpatialPoints(read.csv(paste0("data/", sp, "/", sp, "_filt.csv"))[,1:2],
                       proj4string = CRS(proj))
  
  # Define model
  allf <- list.files(paste0("results/", sp, "/predictions/"), full.names = T)
  
  m <- stringr::str_extract(allf[1], "(?<=m)\\d+")
  
  # Load rasters generated before
  curr <- stack(allf[grepl("current", allf)])
  ssp1 <- stack(allf[grepl("ssp1", allf)])
  ssp2 <- stack(allf[grepl("ssp2", allf)])
  ssp3 <- stack(allf[grepl("ssp3", allf)])
  
  # Convert to data.frame
  get.val <- function(x, add.col = NULL){
    temp <- as(x, "SpatialPixelsDataFrame")
    temp <- as.data.frame(temp)
    colnames(temp) <- c("val", "x", "y")
    if (!is.null(add.col)) {
      temp$group <- add.col 
    }
    return(temp)
  }
  
  # Prepare theme/ploting stuff ----
  # Guide for legend
  step.guide <- guide_colorbar(title = "ROR\n(contrast)",
                               show.limits = TRUE,
                               barwidth = unit(0.16, "in"),
                               barheight = unit(1.5, "in"),
                               ticks = F,
                               ticks.colour = "grey20",
                               frame.colour = "grey20",
                               title.position = "top")
  
  step.guide.int <- guide_colorbar(title = "ROR\n(integral)",
                               show.limits = TRUE,
                               barwidth = unit(0.16, "in"),
                               barheight = unit(1.5, "in"),
                               ticks = F,
                               ticks.colour = "grey20",
                               frame.colour = "grey20",
                               title.position = "top")
  
  step.guide.sd <- guide_colorbar(title = "ROR\n(contrast SD)",
                               show.limits = TRUE,
                               barwidth = unit(0.16, "in"),
                               barheight = unit(1.5, "in"),
                               ticks = F,
                               ticks.colour = "grey20",
                               frame.colour = "grey20",
                               title.position = "top")
  
  
  # Themes
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
          axis.text = element_text(colour = "grey60", size = 12),
          axis.title.x = element_text(vjust = 2, size = 12, colour = "grey60"),
          axis.title.y.left = element_text(size = 12, colour = "grey60"),
          plot.title = element_text(size = 16),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "grey60"),
          panel.grid.major = element_line(linetype = 'dashed', 
                                          colour = "grey70",
                                          size = .1),
          legend.position="right",
          legend.title.align=0.5,
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.background = element_rect(fill = "white"),
          strip.background = element_blank(),
          strip.text = element_text(size = 14)
          
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
  sca <- function(lims, guide){
    scale_fill_gradientn(
      colours = rev(c("#A84C00", "#D97D27", "#F5BD44", "#FFD561", "#FFF291",
                      "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695")),
      limits = lims,
      guide = guide,
      #labels = c("Low", "", "", "", "High"),
      na.value = "#381900")
  }
  
  sca.sd <- function(lims, guide){
    scale_fill_distiller(palette = "GnBu", direction = 1,
                         limits = lims, guide = guide)
  }
  
  
  
  # Generate plots ----
  non.cont <- as.data.frame(rbind(
    get.val(curr[[names(curr)[grepl(paste0("mean_m", m, "_int"), names(curr))]]], "Current"),
    get.val(ssp1[[names(ssp1)[grepl(paste0("mean_m", m, "_int"), names(ssp1))]]], "SSP1"),
    get.val(ssp2[[names(ssp2)[grepl(paste0("mean_m", m, "_int"), names(ssp2))]]], "SSP2"),
    get.val(ssp3[[names(ssp3)[grepl(paste0("mean_m", m, "_int"), names(ssp3))]]], "SSP3")
  ))
  non.cont$variant <- "Mean"
  
  (pint <- ggplot()+
      # Base maps
      geom_sf(data = base,
              color = c("#CFCFCF"),
              size = 0.6,
              fill = "white") +
      # Get the raster
      geom_raster(data = non.cont, aes(x = x, y = y, fill = val)) +
      sca(c(min(non.cont[,1]), max(non.cont[,1])), step.guide.int) + # Add color scale
      # Establish area
      coord_sf(xlim = c(-3239.984, 4510.636),
               ylim = c(-4614.575, 4596.965),
               datum = st_crs(base),
               expand = F,
               label_axes = list(
                 bottom = "E",
                 left = "N"
               )) +
      # Remove axis labels and add theme
      xlab("Easting") + ylab("Northing") + wlt + theme(strip.text.x = element_blank()) +
      ggtitle("Integral") +
      # Add grid
      scale_x_discrete(position = "bottom", breaks = c(-2000, 0, 2000, 4000)) +
      facet_grid(cols = vars(group), rows = vars(variant))
  )
  
  cont <- as.data.frame(rbind(
    get.val(curr[[names(curr)[grepl(paste0("mean_m", m, "_cont"), names(curr))]]], "Current"),
    get.val(ssp1[[names(ssp1)[grepl(paste0("mean_m", m, "_cont"), names(ssp1))]]], "SSP1"),
    get.val(ssp2[[names(ssp2)[grepl(paste0("mean_m", m, "_cont"), names(ssp2))]]], "SSP2"),
    get.val(ssp3[[names(ssp3)[grepl(paste0("mean_m", m, "_cont"), names(ssp3))]]], "SSP3")
  ))
  
  sdr <- as.data.frame(rbind(
    get.val(curr[[names(curr)[grepl(paste0("sd_m", m, "_cont"), names(curr))]]], "Current"),
    get.val(ssp1[[names(ssp1)[grepl(paste0("sd_m", m, "_cont"), names(ssp1))]]], "SSP1"),
    get.val(ssp2[[names(ssp2)[grepl(paste0("sd_m", m, "_cont"), names(ssp2))]]], "SSP2"),
    get.val(ssp3[[names(ssp3)[grepl(paste0("sd_m", m, "_cont"), names(ssp3))]]], "SSP3")
  ))
  
  qhigh <- as.data.frame(rbind(
    get.val(curr[[names(curr)[grepl(paste0("q0.975_m", m, "_cont"), names(curr))]]], "Current"),
    get.val(ssp1[[names(ssp1)[grepl(paste0("q0.975_m", m, "_cont"), names(ssp1))]]], "SSP1"),
    get.val(ssp2[[names(ssp2)[grepl(paste0("q0.975_m", m, "_cont"), names(ssp2))]]], "SSP2"),
    get.val(ssp3[[names(ssp3)[grepl(paste0("q0.975_m", m, "_cont"), names(ssp3))]]], "SSP3")
  ))
  
  qlow <- as.data.frame(rbind(
    get.val(curr[[names(curr)[grepl(paste0("q0.025_m", m, "_cont"), names(curr))]]], "Current"),
    get.val(ssp1[[names(ssp1)[grepl(paste0("q0.025_m", m, "_cont"), names(ssp1))]]], "SSP1"),
    get.val(ssp2[[names(ssp2)[grepl(paste0("q0.025_m", m, "_cont"), names(ssp2))]]], "SSP2"),
    get.val(ssp3[[names(ssp3)[grepl(paste0("q0.025_m", m, "_cont"), names(ssp3))]]], "SSP3")
  ))
  
  cont$variant <- "Mean"
  sdr$variant <- "SD"
  qhigh$variant <- "Q0.975"
  qlow$variant <- "Q0.025"
  
  fulldata <- data.frame(rbind(cont, qhigh, qlow))
  
  fulldata$variant <- factor(fulldata$variant, levels = c("Mean", "Q0.025", "Q0.975"))
  
  # minmax <- range(cont$val[cont$group == "Current"])
  # 
  # fulldata$val <- (fulldata$val - min(minmax))/(max(minmax) - min(minmax))
  options(scipen = 999)
  (pcont <- ggplot()+
      # Base maps
      geom_sf(data = base,
              color = c("#CFCFCF"),
              size = 0.6,
              fill = "white") +
      # Get the raster
      geom_raster(data = fulldata[fulldata$variant == "Mean",], aes(x = x, y = y, fill = val)) +
      sca(c(min(fulldata[,1]), max(fulldata[,1])), step.guide) + # Add color scale
      # Establish area
      coord_sf(xlim = c(-3239.984, 4510.636),
               ylim = c(-4614.575, 4596.965),
               datum = st_crs(base),
               expand = F,
               label_axes = list(
                 bottom = "",
                 left = "N"
               )) +
      # Remove axis labels and add theme
      xlab(NULL) + ylab(NULL) + wlt + theme(legend.position = "none") +
      ggtitle("Contrast") +
      # Add grid
      scale_x_discrete(position = "bottom", breaks = c(-2000, 0, 2000, 4000)) +
      facet_grid(cols = vars(group), rows = vars(variant))
  )
  
  (pcontqlow <- ggplot()+
      # Base maps
      geom_sf(data = base,
              color = c("#CFCFCF"),
              size = 0.6,
              fill = "white") +
      # Get the raster
      geom_raster(data = fulldata[fulldata$variant == "Q0.025",], aes(x = x, y = y, fill = val)) +
      sca(c(min(fulldata[,1]), max(fulldata[,1])), step.guide) + # Add color scale
      # Establish area
      coord_sf(xlim = c(-3239.984, 4510.636),
               ylim = c(-4614.575, 4596.965),
               datum = st_crs(base),
               expand = F,
               label_axes = list(
                 bottom = "",
                 left = "N"
               )) +
      # Remove axis labels and add theme
      xlab(NULL) + ylab(NULL) + wlt + theme(strip.text.x = element_blank()) +
      #ggtitle("Contrast") +
      # Add grid
      scale_x_discrete(position = "bottom", breaks = c(-2000, 0, 2000, 4000)) +
      facet_grid(cols = vars(group), rows = vars(variant))
  )
  
  
  (pcontqhigh <- ggplot()+
      # Base maps
      geom_sf(data = base,
              color = c("#CFCFCF"),
              size = 0.6,
              fill = "white") +
      # Get the raster
      geom_raster(data = fulldata[fulldata$variant == "Q0.975",], aes(x = x, y = y, fill = val)) +
      sca(c(min(fulldata[,1]), max(fulldata[,1])), step.guide) + # Add color scale
      # Establish area
      coord_sf(xlim = c(-3239.984, 4510.636),
               ylim = c(-4614.575, 4596.965),
               datum = st_crs(base),
               expand = F,
               label_axes = list(
                 bottom = "",
                 left = "N"
               )) +
      # Remove axis labels and add theme
      xlab(NULL) + ylab(NULL) + wlt + theme(strip.text.x = element_blank()) + theme(legend.position = "none") +
      #ggtitle("Contrast") +
      # Add grid
      scale_x_discrete(position = "bottom", breaks = c(-2000, 0, 2000, 4000)) +
      facet_grid(cols = vars(group), rows = vars(variant))
  )
  
  (pcontsd <- ggplot()+
      # Base maps
      geom_sf(data = base,
              color = c("#CFCFCF"),
              size = 0.6,
              fill = "white") +
      # Get the raster
      geom_raster(data = sdr, aes(x = x, y = y, fill = val)) +
      sca.sd(c(min(sdr[,1]), max(sdr[,1])), step.guide.sd) + # Add color scale
      # Establish area
      coord_sf(xlim = c(-3239.984, 4510.636),
               ylim = c(-4614.575, 4596.965),
               datum = st_crs(base),
               expand = F,
               label_axes = list(
                 bottom = "",
                 left = "N"
               )) +
      # Remove axis labels and add theme
      xlab(NULL) + ylab(NULL) + wlt +
      theme(strip.text.x = element_blank(), 
            legend.text = element_text(size = 12), 
            legend.title = element_text(size = 12)) +
      #ggtitle("Contrast") +
      # Add grid
      scale_x_discrete(position = "bottom", breaks = c(-2000, 0, 2000, 4000)) +
      facet_grid(cols = vars(group), rows = vars(variant))
  )
  
  
  # Save composite ----
  final <- pcont + pcontqlow + pcontqhigh + pcontsd + pint + plot_layout(ncol = 1)# &
   # theme(legend.position='bottom')
  
  ggsave(paste0("figures/", sp, "_lgcp_supp_comp_m", m, ".jpg"), final,
         width = 25, height = 30, units = "cm", quality = 100)
  
}



# Plot spatial component -----
step.guide <- guide_colorbar(title = "Spatial component effect",
                             show.limits = TRUE,
                             barheight = unit(0.12, "in"),
                             barwidth = unit(3.5, "in"),
                             ticks = F,
                             ticks.colour = "grey20",
                             frame.colour = "grey20",
                             title.position = "top")

sca <- function(lims, guide){
  scale_fill_gradientn(
    colours = rev(c("#A84C00", "#D97D27", "#F5BD44", "#FFD561", "#FFF291",
                    "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695")),
    limits = lims,
    breaks = seq(round(min(lims)), round(max(lims)), by = 0.5),
    guide = guide,
    #labels = c("Low", "", "", "", "High"),
    na.value = "#381900")
}


# Load spatial component
allf <- list.files("results", full.names = T, recursive = T)
allf <- allf[grepl("spatial_mean", allf)]

spatcomp <- stack(allf)

spatcomp.data <- data.frame(rbind(
  get.val(spatcomp$lyva_m4_spatial_mean_effect, "Lytechinus variegatus"),
  get.val(spatcomp$eclu_m4_spatial_mean_effect, "Echinometra lucunter"),
  get.val(spatcomp$trve_m6_spatial_mean_effect, "Tripneustes ventricosus")
))

spatcomp.data$group <- factor(spatcomp.data$group,
                              levels = c("Lytechinus variegatus", "Echinometra lucunter", "Tripneustes ventricosus"))

species.pts <- lapply(species, function(sp){
  occ <- SpatialPoints(read.csv(paste0("data/", sp, "/", sp, "_filt.csv"))[,1:2],
                       proj4string = CRS(proj))
  occ <- sf::st_as_sf(occ)
  occ$group <- ifelse(sp == "lyva", "Lytechinus variegatus",
                      ifelse(sp == "eclu", "Echinometra lucunter", "Tripneustes ventricosus"))
  occ
})
species.pts <- do.call("rbind", species.pts)

(pspat <- ggplot()+
    # Base maps
    geom_sf(data = base,
            color = c("#CFCFCF"),
            size = 0.6,
            fill = "white") +
    # Get the raster
    geom_raster(data = spatcomp.data, aes(x = x, y = y, fill = val)) +
    sca(c(min(spatcomp.data[,1]), max(spatcomp.data[,1])), step.guide) + # Add color scale
    geom_sf(data = species.pts, size = .4, alpha = .2) +
    # Establish area
    coord_sf(xlim = c(-3239.984, 4510.636),
             ylim = c(-4614.575, 4596.965),
             datum = st_crs(base),
             expand = F,
             label_axes = list(
               bottom = "",
               left = "N"
             )) +
    # Remove axis labels and add theme
    xlab(NULL) + ylab(NULL) + wlt + theme(legend.position = "bottom",
                                          legend.text = element_text(size = 9),
                                          legend.title = element_text(size = 11),
                                          strip.text = element_text(size = 14, face = "italic")) +
    #ggtitle("Contrast") +
    # Add grid
    scale_x_discrete(position = "bottom", breaks = c(-2000, 0, 2000, 4000)) +
    facet_grid(cols = vars(group))
)

ggsave(paste0("figures/lgcp_supp_spatial.jpg"), pspat,
       width = 27, height = 14, units = "cm", quality = 100)
