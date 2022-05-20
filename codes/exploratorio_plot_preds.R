test <- st_as_sf(pred.comp$spatial)
test <- st_transform(test, 4326)
pal <- colorNumeric("viridis", c(-2.8, 5.4), na.color = "transparent")
leaflet(test) %>%
        addTiles()%>%
        addCircles(color = pal(test$mean))%>%
        addCircles(lng = pts2$decimalLongitude, lat = pts2$decimalLatitude)
# Predict model without spatial component
pred.comp <- predict(m[[7]], data = df,
                     formula = ~ exp(Intercept + sal + sst + spatial))

pred.comp2 <- predict(m[[8]], data = df,
                      formula = ~ exp(Intercept + sal + sst + spatial))


ggplot()+
        gg(pred.comp2)+ # change here according to layer
        coord_equal()


for (i in 8:10) {
        
        pred.f <- as.formula(paste0(
                "~ exp(",
                paste(if (is.null(names(
                        m[[i]]$summary.random))) {
                        m[[i]]$names.fixed
                } else{
                        c(m[[i]]$names.fixed,
                          names(m[[i]]$summary.random))
                }
                , collapse = "+"), ")"
        ))
        
        pred.t <- predict(m[[i]], data = df,
                          formula = pred.f)
        
        print(ggplot()+
                      gg(pred.t)+
                      gg(pts, color = "orange", size = .2)+
                      scale_fill_viridis()+# change here according to layer
                      coord_equal()+ggtitle(paste0("M ", i)))
        
        
}



occ <- lapply(c("lyva", "eclu", "trve"), function(x){
        read.csv(paste0("data/", x, "/", x, "_cell.csv"))[,1:2]
})

plot(r, ylim = c(-20, -10), xlim = c(-45, -30))
abline(h = -16.3, lwd = 2)
abline(h = -17.8, lwd = 2, col = "red")
abline(h = -12.9, lwd = 2, col = "blue")
points(occ[[1]])
