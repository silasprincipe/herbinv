sst <- raster("data/env/ready_layers/sst_cur.tif")
sst3 <- raster("data/env/ready_layers/sst_r37.tif")

ps1 <- extract(sst, pts)
ps3 <- extract(sst3, pts)


data <- data.frame(rbind(
        data.frame(sst = ps1, scen = "current"),
        data.frame(sst = ps3, scen = "ssp3")
))

data$sit <- ifelse(data$sst < 0, "L", "H")

# VER SE DÁ PRA NÃO DIVIDIR O HALFEYE OU SEJA, APROVEITAR O MESMO FORMATO
# E SO TROCAR A COR!!!!
ggplot(data, aes(x = scen, y = sst
                 #, color = scen
                 )) + 
        ggdist::stat_halfeye(
                adjust = .5,
                width = .3,
                .width = 0,
                justification = -.2,
                point_colour = NA,
                aes(fill = sit),
                fill_type = "segments"
        ) +
        geom_boxplot(
                width = .15, 
                outlier.shape = NA#,
                #aes(fill = scen),
                
        ) +
        geom_point(
                ## draw horizontal lines instead of points
                shape = "|",
                size = 5,
                alpha = .2,
                aes(color = sit)
        ) +
        scale_color_manual(values = c("#0B59BF", "#D65600"))+
        scale_fill_manual(values = c("#0B59BF", "#D65600"))+
        #scale_color_viridis_d()+
        #scale_fill_viridis_d()+
        #coord_cartesian(xlim = c(1.2, NA), clip = "off") +
        coord_flip()
