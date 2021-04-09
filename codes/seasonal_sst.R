files <- list.files("data/oisst/", pattern = "oisst_noaa", full.names = T)

oisst <- stack(files)

cropbox <-c(-60,-40,-3,20)   #colocar os limites da area de estudo
adg443 <- crop(adg443, cropbox)
adg443[adg443 < 0]= NA

# atribuir uma data para cada raster stacked
dates <- seq(as.Date('2010-01-01'), as.Date('2020-12-31'), 'day') #colocar a serie temporal
names(oisst) <- dates

dim(oisst)

#obter a data dos nomes das camadas e extrair o numero do mes
indices <- format(as.Date(names(oisst), format = "X%Y.%m.%d"), format = "%m")
indices <- as.numeric(indices)

#sum layers
monthly.mean <- stackApply(oisst, indices, fun = mean) #Essa Função vai somar todos os meses de índice igual
names(monthly.mean) <- format(seq(as.Date('2010-01-01'), as.Date('2010-12-31'), 'month'), "%m")


mmean.data <- data.frame(rasterToPoints(monthly.mean))

mmean.data <- pivot_longer(mmean.data,
                           3:14,
                           names_to = "month",
                           values_to = "sst_mean")

getMonth <- function(month, sst_mean){
        x <- month[which(sst_mean == max(sst_mean))]
        x
}

max.data <- mmean.data %>%
        group_by(x, y) %>%
        summarise("max" = max(sst_mean),
                  "month" = getMonth(month, sst_mean)) %>%
        mutate(month = str_remove(month, "X"))

season.sst <- rasterFromXYZ(max.data[,1:3])
maxmonth.sst <- rasterFromXYZ(max.data[,c(1:2, 4)])

plot(season.sst)

# We calculate the max just to compare

max.sst <- calc(oisst, max)

crs(season.sst) <- crs("+proj=longlat +datum=WGS84 +no_defs")
crs(maxmonth.sst) <- crs("+proj=longlat +datum=WGS84 +no_defs")


#Save

writeRaster(season.sst, "data/oisst/sst_maxhotmonth.tif")
writeRaster(maxmonth.sst, "data/oisst/sst_maxhotmonth_hm.tif")
writeRaster(max.sst, "data/oisst/sst_max.tif")
