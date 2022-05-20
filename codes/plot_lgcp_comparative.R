ssp2.v$y <- round(ssp2.v$y)

# Convert to data.frame
get.val <- function(x){
        temp <- as(x, "SpatialPixelsDataFrame")
        temp <- as.data.frame(temp)
        colnames(temp) <- c("val", "x", "y")
        temp$val <- (temp$val - min(temp$val))/(max(temp$val) - min(temp$val))
        temp
}

curr.v <- get.val(curr)
ssp1.v <- get.val(ssp1)
ssp2.v <- get.val(ssp2)
ssp3.v <- get.val(ssp3)


curr.v <- curr.v %>% group_by(y) %>%
        summarise(mvl = mean(val),
                  upq = quantile(val, .95),
                  lpq = quantile(val, .05))
ssp1.v <- ssp1.v %>% group_by(y) %>%
        summarise(mvl = mean(val),
                  upq = quantile(val, .95),
                  lpq = quantile(val, .05))
ssp2.v <- ssp2.v %>% group_by(y) %>%
        summarise(mvl = mean(val),
                  upq = quantile(val, .95),
                  lpq = quantile(val, .05))
ssp3.v <- ssp3.v %>% group_by(y) %>%
        summarise(mvl = mean(val),
                  upq = quantile(val, .95),
                  lpq = quantile(val, .05))

total <- rbind(
        cbind(curr.v, scen = "curr"),
        cbind(ssp1.v, scen = "ssp1"),
        cbind(ssp2.v, scen = "ssp2"),
        cbind(ssp3.v, scen = "ssp3")
)

total <- as.data.frame(total)

ggplot(total)+
        geom_line(aes(x = y, y = mvl))+
        geom_ribbon(aes(x = y, ymin = lpq, ymax = upq), alpha = .5)+
        coord_flip()+
        facet_wrap(~scen, nrow = 1, ncol =4)

tit <- data.frame(y = rep(1, 4),
                  x = rep(4900, 4),
                  t = c("Current", "SSP1", "SSP2", "SSP3"),
                  scen = c("curr", "ssp1", "ssp2", "ssp3"))

stguide <- guide_coloursteps(title = "Mean intensity",
                             show.limits = TRUE,
                             barheight = unit(0.1, "in"),
                             barwidth = unit(1.1, "in"),
                             ticks = T,
                             ticks.colour = "grey20",
                             frame.colour = "grey20",
                             title.position = "top")

ggplot(total)+
        geom_tile(aes(x = y, y = 1, fill = mvl))+
        scale_fill_stepsn(breaks = seq(0.1,0.9,length.out = 9),
                          limits = c(0, 1),
                          colors = c(
                                  "#F7FCFD", "#E0ECF4", "#BFD3E6", "#9EBCDA",
                                  "#8C96C6", "#8C6BB1",
                                  "#88419D", "#810F7C", "#4D004B", "#2F002E"
                          ),
                          na.value = NA,
                          guide = stguide,
                          labels = function(x){
                                  ifelse(x != 0 & x != 0.5 & x != 1, "", x)
                                  })+
        geom_text(data = tit, aes(x = x, y = y, label = t), angle = 90,
                  hjust = 0)+
        coord_flip()+
        xlab("Northing")+ ylab(NULL)+
        theme_classic()+
        scale_x_continuous(limits = c(-4700, 6500),
                           expand = c(0,0),
                           breaks = seq(-4000, 4000, by = 1000))+
        theme(
                panel.grid = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.line = element_blank(),
                strip.background = element_blank(),
                strip.text.x = element_blank(),
                legend.position="bottom"
        )+
        facet_wrap(~scen, nrow = 1, ncol = 4)

ggsave(paste0("figures/", sp, "_lgcp_comparative.jpg"),
       width = 5, height = 18, units = "cm")
