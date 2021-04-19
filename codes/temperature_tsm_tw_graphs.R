species <- c("lyva", "eclu", "trve")
round.code <- "ih_c2"


temps <- data.frame("species" = c("lyva", "eclu", "trve"),
                    "tmax" = c(34.5, 36, 34),
                    "topt" = c(27.2, 29.4, 30.7))

all.data <- list()

for (i in 1:length(species)) {
        
        curr <- raster(paste0(
                species[i],
                "/proj_current_", species[i],"_",
                round.code,"/individual_projections/",
                species[i], "_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img"
        ))
        
        thr.app <- function(x){
                x[x < 667] <- 0
                x[x >= 667] <- 1
                x
        }
        
        curr <- calc(curr, thr.app)
        
        pts <- data.frame(rasterToPoints(curr))
        
        pts <- pts[pts[,3] == 1,]
        
        colnames(pts) <- c("x", "y", "val")
        
        samp.pts <- sample_n(pts, 1000)
        
        sst.r <- lapply(c("current", "ssp126", "ssp245", "ssp585"), function(x){
                
                if (x == "current") {
                        r <- raster("data/env/crop_layers/BO21_tempmean_ss.tif")
                }
                
                if (x != "current") {
                        r <- raster(paste0("data/env/proj_layers/",
                                           x,
                                           "/BO21_tempmean_ss.tif"))
                }
                
                r
        })
        
        
        tws <- lapply(sst.r, function(rast){
                
                tdata <- raster::extract(rast,
                                         samp.pts[,1:2])
                
                tdata <- (tdata - temps[i, 2]) * -1
                
                tdata <- data.frame(cbind(samp.pts[,1:2],
                                          tdata))
                
                tdata[,3] <- cut(tdata[,3], breaks = c(-10, 0, 2, 5, 30),
                                 include.lowest = T)
                
                levels(tdata[,3]) <- c("< 0",
                                       "0 - 2",
                                       "2 - 5",
                                       "> 5")
                tdata
        })
        
        
        
        tsms <- lapply(sst.r, function(rast){
                
                tdata <- raster::extract(rast,
                                         samp.pts[,1:2])
                
                tdata <- (tdata - temps[i, 3]) * -1
                
                tdata <- data.frame(cbind(samp.pts[,1:2],
                                          tdata))
                
                # tdata[,3] <- cut(tdata[,3], breaks = c(-10, -1, 0, 2, 5, 10),
                #                  include.lowest = T)
                # 
                # levels(tdata[,3]) <- c("< -1",
                #                        "-1 - 0",
                #                        "0 - 2",
                #                        "2 - 5",
                #                        "> 5")
                
                tdata$class <- ifelse(tdata[,3] < 0,
                                      "R",
                                      "N")
                
                tdata$group <- names(rast)
                
                tdata
        })
        
        all.data[[i]] <- list(tws, tsms)
        
        names(all.data[[i]]) <- c("tw", "tsm")
        
}

names(all.data) <- species


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
                legend.key=element_blank(),
                # legend.title = element_text(size = 10),
                # legend.text = element_text(size = 10),
                #legend.key.height = unit(7, 'mm'),
                #legend.key.width = unit(7, 'mm'),
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

step.guide <- guide_coloursteps(title = "TSM",
                                show.limits = TRUE,
                                barheight = unit(2.4, "in"),
                                barwidth = unit(0.18, "in"),
                                ticks = T,
                                ticks.colour = "grey20",
                                frame.colour = "grey20",
                                title.vjust = unit(0.1, "in"))


teste <- all.data$lyva$tsm[[1]] %>% 
        mutate(binned = cut(tdata, breaks = c(-10, -1, 0, 2, 5, 10),
                            include.lowest = T))

levels(teste$binned) <- c("< -1",
                          "-1 - 0",
                          "0 - 2",
                          "2 - 5",
                          "> 5")

colscale <- scale_color_manual(values = c("#C10000",
                                          "#FF7F11",
                                          "#6CA4FF",
                                          "#0053D7",
                                          "#002E77"
                                          ), name = "TSM",
                               breaks = c("< -1",
                                          "-1 - 0",
                                          "0 - 2",
                                          "2 - 5",
                                          "> 5"),
                               labels = c("< -1",
                                          "-1 - 0",
                                          "0 - 2",
                                          "2 - 5",
                                          "> 5"),
                               drop = FALSE)




la <- list(list(bottom = "", left = "N", top = "E"),
           list(bottom = "E", right = "", top = ""),
           list(bottom = "", left = "", top = "E"),
           list(bottom = "E", right = "N", top = "")
           )

for (z in 1:length(species)) {
        
        tsm.plotl <- list()
        
        for (i in 1:4) {
                
                p <- ggplot()+
                        # Base maps
                        geom_sf(data = base2,
                                color = c.sea,
                                fill = c.sea) +
                        geom_sf(data = base1,
                                color = c.land,
                                fill = c.land) +
                        # Add points
                        geom_point(data = all.data[[species[z]]]$tsm[[i]],
                                   aes(x = x, y = y, color = tdata),
                                   size = 1)+
                        colscale+
                        # Establish area
                        coord_sf(
                                xlim = c(-99, -29),
                                ylim = c(-42.5, 42.5),
                                expand = FALSE,
                                label_axes = la[[i]]
                        ) +
                        # Add plot details
                        xlab(NULL) + ylab(NULL)+
                        main.theme + main.s.x + main.s.y +
                        guides(colour = guide_legend(override.aes = list(size = 3)))
                
                if (i != 1) {
                        p <- p+theme(legend.position = 'none')
                }
                
                tsm.plotl[[i]] <- p
                
        }
        
        pf <- tsm.plotl[[1]] | tsm.plotl[[2]] | tsm.plotl[[3]] | tsm.plotl[[4]]
        
        
        ggsave(paste0("figures/", species[z], "_tsm.tiff"),
               pf,
               width = 50, height = 15, dpi = 300, units = 'cm')
}



### TW

colscale <- scale_color_manual(values = c("#C10000",
                                          "#9DC2FF",
                                          "#3B86FF",
                                          "#0053D7",
                                          "#002E77"), 
                                name = "TW",
                                breaks = c("< 0",
                                           "0 - 2",
                                           "2 - 5",
                                           "> 5"),
                                labels = c("< 0",
                                           "0 - 2",
                                           "2 - 5",
                                           "> 5"),
                                drop = FALSE)


for (z in 1:length(species)) {
        
        tw.plotl <- list()
        
        for (i in 1:4) {
                
                p <- ggplot()+
                        # Base maps
                        geom_sf(data = base2,
                                color = c.sea,
                                fill = c.sea) +
                        geom_sf(data = base1,
                                color = c.land,
                                fill = c.land) +
                        # Add points
                        geom_point(data = all.data[[species[z]]]$tw[[i]],
                                   aes(x = x, y = y, color = tdata),
                                   size = 1)+
                        colscale+
                        # Establish area
                        coord_sf(
                                xlim = c(-99, -29),
                                ylim = c(-42.5, 42.5),
                                expand = FALSE,
                                label_axes = la[[i]]
                        ) +
                        # Add plot details
                        xlab(NULL) + ylab(NULL)+
                        main.theme + main.s.x + main.s.y +
                        guides(colour = guide_legend(override.aes = list(size = 3)))
                
                if (i != 1) {
                        p <- p+theme(legend.position = 'none')
                }
                
                tw.plotl[[i]] <- p
                
        }
        
        pf <- tw.plotl[[1]] | tw.plotl[[2]] | tw.plotl[[3]] | tw.plotl[[4]]
        
        
        ggsave(paste0("figures/", species[z], "_tw.tiff"),
               pf,
               width = 50, height = 15, dpi = 300, units = 'cm')
}





main.theme <-
        theme(
                panel.background = element_rect(fill = "white", 
                                                colour = "gray50"),
                # panel.grid.major = element_line(linetype = 'solid', 
                #                                 colour = "lightgrey"),
                panel.grid.major = element_blank(), #If wanted to remove line grid
                panel.grid.minor = element_blank(),
                axis.ticks.length = unit(0.8, "mm"),
                axis.text = element_text(size = 10),
                #legend.justification = c(1, 0),
                #legend.position = c(0.18, 0.02),
                #legend.background = element_blank(),
                # legend.title = element_text(size = 10),
                # legend.text = element_text(size = 10),
                # legend.key.height = unit(7, 'mm'),
                # legend.key.width = unit(7, 'mm'),
                # legend.key = element_rect(colour = "#696969", size =
                #                                   2.5),
                plot.margin = unit(c(0.3, 0.35, 0.3, 0.35), 'mm'),
                strip.background = element_blank(),
                axis.ticks.y.right = element_blank(),
                axis.title.y.left = element_blank(),
                axis.text.y.right = element_blank(),
                axis.title.y.right = element_text(face = "italic",
                                                  margin = margin(l = 10))
                )


# Maps scales
main.s.x <-
        scale_x_continuous(breaks = seq(-5,7.5,2.5),
                           limits = c(-5.5,8))
main.s.y <-
        scale_y_continuous(breaks = seq(-40, 40, 10),
                           labels = c(paste0(seq(-40, -10, 10), "°S"),
                                      "0°",
                                      paste0(seq(10, 40, 10), "°N")),
                           limits = c(-42.5, 42.5),
                           sec.axis = dup_axis())

all.plots <- list()


for (i in 1:3) {
        
        all.data[[species[i]]]$tsm[[1]]$group <- "Current" 
        all.data[[species[i]]]$tsm[[2]]$group <- "SSP 1"
        all.data[[species[i]]]$tsm[[3]]$group <- "SSP 2"
        all.data[[species[i]]]$tsm[[4]]$group <- "SSP 5"
        
        data.full <- rbind(all.data[[species[i]]]$tsm[[1]],
                           all.data[[species[i]]]$tsm[[2]],
                           all.data[[species[i]]]$tsm[[3]],
                           all.data[[species[i]]]$tsm[[4]])
        
        means <- data.full %>% group_by(group) %>%
                summarise("mean" = mean(tdata))
        
        
        p <-    ggplot(data.full)+
                geom_point(aes(x = tdata, y = y, color = class),
                           size = 1, alpha = 0.5)+
                scale_color_manual(values = c("#1874CD", "#EE7600"))+
                geom_vline(xintercept = 0, color = "gray30")+
                geom_vline(data = means, aes(xintercept = mean), color = "#1874CD",
                           linetype = 2)+
                main.theme +
                theme(legend.position = "none")+
                main.s.y+
                ylab(c("Lytechinus variegatus",
                       "Echinometra lucunter",
                       "Tripneustes ventricosus")[i])+
                main.s.x+
                xlab(NULL)+
                facet_wrap(~ group, nrow = 1)
        
        if (i == 1 | i == 2) {
                
                p <- p + theme(axis.text.x=element_blank(),
                               axis.ticks.x=element_blank())
        }
        
        if (i != 1) {
                
                p <- p + theme(strip.text.x = element_blank())
        }
        
        all.plots[[i]] <- p
        
}

fp <- all.plots[[1]] / all.plots[[2]] / all.plots[[3]]

ggsave(paste0("figures/tsm_points.tiff"),
       fp,
       width = 30, height = 25, dpi = 300, units = 'cm')

