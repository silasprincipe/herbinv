#### Modelling of coral reef herbivorous invertebrates ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2022 ##

# Plot of the MESS maps and difference in temperature

# Load needed packages ----
library(ecospat)
library(inlabru)
library(ggplot2)
library(sf)
library(raster)
library(patchwork)

# Set species
species <- c("lyva", "eclu", "trve")

# Load quadrature points
list2env(readRDS('data/lgcp_data.rds'),globalenv())

ips <- ipoints(starea, mesh)

# Load environmental data
env <- stack(list.files("data/env/ready_layers", pattern = "_cur", full.names = T))
r12 <- stack(list.files("data/env/ready_layers", pattern = "_r12", full.names = T))
r24 <- stack(list.files("data/env/ready_layers", pattern = "_r24", full.names = T))
r37 <- stack(list.files("data/env/ready_layers", pattern = "_r37", full.names = T))

# Include the distance to coast layer
env <- stack(env, raster("data/env/ready_layers/distcoast.tif"))
r12[[8]] <- r24[[8]] <- r37[[8]] <- env[[8]]

# Change names
names(env) <- c("chl", "coldm", "ph", "pho", "sal", "sst", "warmm", "dist")
names(r12) <- names(r24) <- names(r37) <- names(env)

# Load non-scaled SST
env.temp <- stack("data/env/crop_layers/BO21_tempmean_ss.tif",
                  "data/env/bioclim_layers/mtemp_coldm_current_hr.tif")
r12.temp <- stack("data/env/proj_layers/ssp126/BO21_tempmean_ss.tif",
                  "data/env/bioclim_layers/mtemp_coldm_ssp126_hr.tif")
r24.temp <- stack("data/env/proj_layers/ssp245/BO21_tempmean_ss.tif",
                  "data/env/bioclim_layers/mtemp_coldm_ssp245_hr.tif")
r37.temp <- stack("data/env/proj_layers/ssp370/BO21_tempmean_ss.tif",
                  "data/env/bioclim_layers/mtemp_coldm_ssp370_hr.tif")

names(env.temp) <- names(r12.temp) <- names(r24.temp) <- names(r37.temp) <- c("sst", "coldm")

env.temp <- mask(env.temp, env.temp[[1]])
r12.temp <- mask(r12.temp, r12.temp[[1]])
r24.temp <- mask(r24.temp, r24.temp[[1]])
r37.temp <- mask(r37.temp, r37.temp[[1]])

env.temp <- projectRaster(env.temp, crs = crs(env))
r12.temp <- projectRaster(r12.temp, crs = crs(env))
r24.temp <- projectRaster(r24.temp, crs = crs(env))
r37.temp <- projectRaster(r37.temp, crs = crs(env))

r12.temp.delta <- r12.temp - env.temp
r24.temp.delta <- r24.temp - env.temp
r37.temp.delta <- r37.temp - env.temp

r12.temp.delta <- rasterToPoints(r12.temp.delta)
r24.temp.delta <- rasterToPoints(r24.temp.delta)
r37.temp.delta <- rasterToPoints(r37.temp.delta)

env.temp <- rasterToPoints(env.temp)

# Set layers used in each model
layers <- list(
  c("sal", "dist", "sst", "ph"),
  c("sal", "dist", "coldm", "ph"),
  c("sal", "dist", "coldm", "ph", "chl", "pho")
)

# Set function to get data for integration points using IDW
getd <- function(rast, ip){
  library(gstat)
  rast <- extend(rast, (extent(min(mesh$loc[,1]),
                               max(mesh$loc[,1]),
                               min(mesh$loc[,2]),
                               max(mesh$loc[,2])))+
                   c(-200, 200, -200, 200))
  
  epts <- raster::extract(rast, ip)
  epts <- data.frame(epts, coordinates(ip))
  
  tofill <- epts[is.na(epts[,1]),]
  
  epts <- data.frame(rasterToPoints(rast))
  
  coordinates(epts) <- ~x+y
  coordinates(tofill) <- ~x+y
  
  for (i in 1:nlayers(rast)) {
    
    mod <- gstat(formula = as.formula(paste(names(epts)[i],"~ 1")),
                 data = epts, nmax = 12)
    
    pred <- predict(mod, tofill)
    
    rast[[i]][cellFromXY(rast[[i]], coordinates(tofill))] <- pred$var1.pred
  }
  
  rast
}

env <- getd(env, ip = ips)
r12 <- getd(r12, ip = ips)
r24 <- getd(r24, ip = ips)
r37 <- getd(r37, ip = ips)

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

base.env <- stack(list.files("data/env/ready_layers", pattern = "_cur", full.names = T)[1])


# Get mess for each species ----
for (i in 1:length(species)) {
  
  sp <- species[i]
  
  # Load occurrence data
  proj <- "+proj=laea +lat_0=0 +lon_0=-70 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs"
  occ <- SpatialPoints(read.csv(paste0("data/", sp, "/", sp, "_filt.csv"))[,1:2],
                       proj4string = CRS(proj))
  
  # Select used environmental layers
  env.used <- env[[layers[[i]]]]
  r12.used <- r12[[layers[[i]]]]
  r24.used <- r24[[layers[[i]]]]
  r37.used <- r37[[layers[[i]]]]
  
  # Get environmental information on points
  env.used.pts <- raster::extract(env.used, rbind(occ, ips))
  
  env.mess <- ecospat.mess(rasterToPoints(env.used), cbind(coordinates(rbind(occ, ips)), env.used.pts))
  r12.mess <- ecospat.mess(rasterToPoints(r12.used), cbind(coordinates(rbind(occ, ips)), env.used.pts))
  r24.mess <- ecospat.mess(rasterToPoints(r24.used), cbind(coordinates(rbind(occ, ips)), env.used.pts))
  r37.mess <- ecospat.mess(rasterToPoints(r37.used), cbind(coordinates(rbind(occ, ips)), env.used.pts))
  
  env.mess.t <- base.env
  r12.mess.t <- base.env
  r24.mess.t <- base.env
  r37.mess.t <- base.env
  
  env.mess.t[cellFromXY(env.mess.t, env.mess[,1:2])] <- env.mess[,5]
  r12.mess.t[cellFromXY(r12.mess.t, r12.mess[,1:2])] <- r12.mess[,5]
  r24.mess.t[cellFromXY(r24.mess.t, r24.mess[,1:2])] <- r24.mess[,5]
  r37.mess.t[cellFromXY(r37.mess.t, r37.mess[,1:2])] <- r37.mess[,5]
  
  env.mess <- mask(env.mess.t, base.env)
  r12.mess <- mask(r12.mess.t, base.env)
  r24.mess <- mask(r24.mess.t, base.env)
  r37.mess <- mask(r37.mess.t, base.env)
  
  env.mess <- as.data.frame(rasterToPoints(env.mess))
  r12.mess <- as.data.frame(rasterToPoints(r12.mess))
  r24.mess <- as.data.frame(rasterToPoints(r24.mess))
  r37.mess <- as.data.frame(rasterToPoints(r37.mess))
  
  env.mess[,3] <- ifelse(env.mess[,3] > 0, "Extrapolation", "No extrapolation")
  r12.mess[,3] <- ifelse(r12.mess[,3] > 0, "Extrapolation", "No extrapolation")
  r24.mess[,3] <- ifelse(r24.mess[,3] > 0, "Extrapolation", "No extrapolation")
  r37.mess[,3] <- ifelse(r37.mess[,3] > 0, "Extrapolation", "No extrapolation")
  
  colnames(env.mess)[3] <- "MESSneg"
  colnames(r12.mess)[3] <- "MESSneg"
  colnames(r24.mess)[3] <- "MESSneg"
  colnames(r37.mess)[3] <- "MESSneg"
  
  # Load predictions species
  preds <- list.files(paste0("results/", sp, "/predictions"), full.names = T)
  preds <- preds[grepl("cont", preds)]
  preds <- preds[grepl("mean", preds)]
  
  env.pred <- stack(preds[grepl("current", preds)])
  r12.pred <- stack(preds[grepl("ssp1", preds)])
  r24.pred <- stack(preds[grepl("ssp2", preds)])
  r37.pred <- stack(preds[grepl("ssp3", preds)])
  
  maxval <- maxValue(env.pred)
  
  r12.pred[r12.pred <= maxval] <- NA
  r24.pred[r24.pred <= maxval] <- NA
  r37.pred[r37.pred <= maxval] <- NA
  
  r12.pred[!is.na(r12.pred)] <- 1
  r24.pred[!is.na(r24.pred)] <- 1
  r37.pred[!is.na(r37.pred)] <- 1
  
  r12.pred <- aggregate(buffer(rasterToPolygons(r12.pred, dissolve = T), 0.0001))
  r24.pred <- aggregate(buffer(rasterToPolygons(r24.pred, dissolve = T), 0.0001))
  r37.pred <- aggregate(buffer(rasterToPolygons(r37.pred, dissolve = T), 0.0001))
  
  r12.pred <- sf::st_as_sf(r12.pred)
  r24.pred <- sf::st_as_sf(r24.pred)
  r37.pred <- sf::st_as_sf(r37.pred)
  
  # Prepare theme/ploting stuff ----
  # Guide for legend
  step.guide <- guide_legend(title = NULL,
                             # show.limits = TRUE,
                             keyheight = unit(0.2, "in"),
                             keywidth = unit(0.7, "in"),
                             ticks = T,
                             ticks.colour = "grey20",
                             frame.colour = "grey20",
                             title.position = "left",
                             title.hjust = 0.5,
                             label.position = "bottom",
                             label.hjust = 0.5)
  
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
          axis.text = element_text(colour = "grey60", size = 14),
          axis.title.x.top = element_text(vjust = 2, size = 16),
          axis.title.y.left = element_text(size = 16),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "grey60"),
          panel.grid.major = element_line(linetype = 'dashed', 
                                          colour = "grey70",
                                          size = .1),
          legend.position="bottom",
          legend.title.align=0.5,
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          legend.background = element_rect(fill = "white"),
          legend.spacing.x = unit(0.1, 'cm')
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
  
  sca <- scale_fill_manual(values = c("#DC267F", "grey80"),
                           name = NULL,
                           #guide = step.guide, 
                           na.value = "grey70")
  
  # sca <- scale_fill_stepsn(breaks = c(),
  #                          limits = c(8, 30),
  #                          colors = rev(RColorBrewer::brewer.pal(10, "RdYlBu")),
  #                          na.value = NA,
  #                          guide = step.guide)
  # 
  
  # Generate plots ----
  ##### Current ----
  (pc <- ggplot()+
      # Base maps
      geom_sf(data = base,
              color = c("#CFCFCF"),
              size = 0.6,
              fill = "white") +
      # Get the raster
      geom_raster(data = env.mess, aes(x = x, y = y, fill = as.factor(MESSneg))) +
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
      geom_label(aes(label = "A", x = 4000, y = 4100),
                 size = 10, fontface = "bold", color = "grey30",
                 label.size = 0)+
      # Remove axis labels and add theme
      xlab(NULL) + ylab("Northing") + wlt +
      # Add grid
      scale_x_discrete(position = "top", breaks = c(-2000, 0, 2000, 4000))
  )
  
  
  ##### SSP1 ----
  (ps1 <- ggplot()+
     # Base maps
     geom_sf(data = base,
             color = c("#CFCFCF"),
             size = 0.6,
             fill = "white") +
     # Get the raster
     geom_raster(data = r12.mess, aes(x = x, y = y, fill = as.factor(MESSneg))) +
     sca + # Add color scale
     geom_sf(data = r12.pred, color = "grey20", fill = NA) +
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
     # Add title
     geom_label(aes(label = "B", x = 4000, y = 4100),
                size = 10, fontface = "bold", color = "grey30",
                label.size = 0)+
     # Remove axis labels and add theme
     xlab(NULL) + ylab(NULL) + wlt +
     # Add grid
     scale_x_discrete(position = "top", breaks = c(-2000, 0, 2000, 4000))
  )
  
  
  
  ##### SSP2 ----
  (ps2 <- ggplot()+
     # Base maps
     geom_sf(data = base,
             color = c("#CFCFCF"),
             size = 0.6,
             fill = "white") +
     # Get the raster
     geom_raster(data = r24.mess, aes(x = x, y = y, fill = as.factor(MESSneg))) +
     sca + # Add color scale
     geom_sf(data = r24.pred, color = "grey20", fill = NA) +
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
     # Add title
     geom_label(aes(label = "C", x = 4000, y = 4100),
                size = 10, fontface = "bold", color = "grey30",
                label.size = 0)+
     # Remove axis labels and add theme
     xlab(NULL) + ylab(NULL) + wlt +
     # Add grid
     scale_x_discrete(position = "top", breaks = c(-2000, 0, 2000, 4000)) 
  )
  
  ##### SSP3 ----
  (ps3 <- ggplot()+
     # Base maps
     geom_sf(data = base,
             color = c("#CFCFCF"),
             size = 0.6,
             fill = "white") +
     # Get the raster
     geom_raster(data = r37.mess, aes(x = x, y = y, fill = as.factor(MESSneg))) +
     sca + # Add color scale
     geom_sf(data = r37.pred, color = "grey20", fill = NA) +
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
     # Add title
     geom_label(aes(label = "D", x = 4000, y = 4100),
                size = 10, fontface = "bold", color = "grey30",
                label.size = 0)+
     # Remove axis labels and add theme
     xlab("Easting") + ylab(NULL) + wlt +
     # Add grid
     scale_x_discrete(position = "top", breaks = c(-2000, 0, 2000, 4000)) 
  )
  
  
  # Save composite ----
  final <- pc + ggtitle("MESS") + theme(plot.title = element_text(size = 20)) +
    ps1 + ps2 + ps3 + plot_layout(nrow = 1, guides = "collect") &
    theme(legend.position='bottom') 
  
  # ggsave(paste0(paste0("figures/mess_", sp, ".jpg")), final,
  #        width = 50, height = 18, units = "cm", quality = 100)
  # 
  
  
  
  ### Temperature plots ----
  step.guide.a <- guide_colorbar(title = "Temperature (°C)",
                               show.limits = TRUE,
                               barheight = unit(0.12, "in"),
                               barwidth = unit(3.5, "in"),
                               ticks = F,
                               ticks.colour = "grey20",
                               frame.colour = "grey20",
                               title.position = "top",
                               label.theme = element_text(size = 10))
  
  sca.a <- scale_fill_stepsn(breaks = if ("sst" %in% names(env.used)) {seq(8,30,2)} else {seq(1,28.5,2.5)},
                           limits = if ("sst" %in% names(env.used)) {c(8,30)} else {c(1,28.5)},
                           colors = rev(RColorBrewer::brewer.pal(10, "RdYlBu")),
                           na.value = NA,
                           guide = step.guide.a)
  
  step.guide.b <- guide_colorbar(title = "Difference in temperature (°C)",
                                 show.limits = TRUE,
                                 barheight = unit(0.12, "in"),
                                 barwidth = unit(3.5, "in"),
                                 ticks = F,
                                 ticks.colour = "grey20",
                                 frame.colour = "grey20",
                                 title.position = "top")
  
  sca.b <- scale_fill_distiller(
    type = "seq",
    palette = "RdPu",
    direction = 1,
    limits = if ("sst" %in% names(env.used)) {c(0.3,4.2)} else {c(0.6,6.7)},
    na.value = NA,
    guide = step.guide.b)
  
  ##### Current ----
  var <- if ("sst" %in% names(env.used)) {"sst"} else {"coldm"}
  temp.comp <- as.data.frame(env.temp[,c("x", "y", var)])
  colnames(temp.comp)[3] <- "values"
  
  (pc <- ggplot()+
     # Base maps
     geom_sf(data = base,
             color = c("#CFCFCF"),
             size = 0.6,
             fill = "white") +
     # Get the raster
     geom_raster(data = temp.comp,
                 aes(x = x, y = y, fill = values)) +
     sca.a + # Add color scale
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
     geom_label(aes(label = "E", x = 4000, y = 4100),
                size = 10, fontface = "bold", color = "grey30",
                label.size = 0)+
     # Remove axis labels and add theme
     xlab(NULL) + ylab(NULL) + wlt +
     # Add grid
     scale_x_discrete(position = "top", breaks = c(-2000, 0, 2000, 4000))
  )
  
  
  ##### SSP1 ----
  temp.comp <- as.data.frame(r12.temp.delta[,c("x", "y", var)])
  colnames(temp.comp)[3] <- "values"
  
  (ps1 <- ggplot()+
     # Base maps
     geom_sf(data = base,
             color = c("#CFCFCF"),
             size = 0.6,
             fill = "white") +
     # Get the raster
     geom_raster(data = temp.comp,
                  aes(x = x, y = y, fill = values)) +
     sca.b + # Add color scale
     geom_sf(data = r12.pred, color = "grey20", fill = NA) +
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
     # Add title
     geom_label(aes(label = "F", x = 4000, y = 4100),
                size = 10, fontface = "bold", color = "grey30",
                label.size = 0)+
     # Remove axis labels and add theme
     xlab(NULL) + ylab(NULL) + wlt +
     # Add grid
     scale_x_discrete(position = "top", breaks = c(-2000, 0, 2000, 4000))
  )
  
  
  
  ##### SSP2 ----
  temp.comp <- as.data.frame(r24.temp.delta[,c("x", "y", var)])
  colnames(temp.comp)[3] <- "values"
  
  (ps2 <- ggplot()+
     # Base maps
     geom_sf(data = base,
             color = c("#CFCFCF"),
             size = 0.6,
             fill = "white") +
     # Get the raster
     geom_raster(data = temp.comp,
                 aes(x = x, y = y, fill = values)) +
     sca.b + # Add color scale
     geom_sf(data = r24.pred, color = "grey20", fill = NA) +
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
     # Add title
     geom_label(aes(label = "G", x = 4000, y = 4100),
                size = 10, fontface = "bold", color = "grey30",
                label.size = 0)+
     # Remove axis labels and add theme
     xlab(NULL) + ylab(NULL) + wlt +
     # Add grid
     scale_x_discrete(position = "top", breaks = c(-2000, 0, 2000, 4000)) 
  )
  
  ##### SSP3 ----
  temp.comp <- as.data.frame(r37.temp.delta[,c("x", "y", var)])
  colnames(temp.comp)[3] <- "values"
  
  (ps3 <- ggplot()+
     # Base maps
     geom_sf(data = base,
             color = c("#CFCFCF"),
             size = 0.6,
             fill = "white") +
     # Get the raster
     geom_raster(data = temp.comp,
                 aes(x = x, y = y, fill = values)) +
     sca.b + # Add color scale
     geom_sf(data = r37.pred, color = "grey20", fill = NA) +
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
     # Add title
     geom_label(aes(label = "H", x = 4000, y = 4100),
                size = 10, fontface = "bold", color = "grey30",
                label.size = 0)+
     # Remove axis labels and add theme
     xlab(NULL) + ylab(NULL) + wlt +
     # Add grid
     scale_x_discrete(position = "top", breaks = c(-2000, 0, 2000, 4000)) 
  )
  
  
  # Save composite ----
  final.temp <- pc + ggtitle("Temperature") + theme(plot.title = element_text(size = 20)) +
    ps1 + ps2 + ps3 + plot_layout(nrow = 1, guides = "collect") &
    theme(legend.position='bottom') 
  
  final.both <- final / final.temp
  
  ggsave(paste0(paste0("figures/mess_temperature_", sp, ".jpg")), final.both,
         width = 50, height = 36, units = "cm", quality = 100)
  
}
