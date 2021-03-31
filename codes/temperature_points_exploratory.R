#### Thermal limits ####

# Load libraries ----

library(raster)
library(tidyverse)
library(boot)
library(rasterVis)
library(RColorBrewer)

# rasterVis plot parameters
mapTheme <- rasterTheme(region = rev(brewer.pal(10, "RdBu")))

# Load species data ----

sp.data <- read.csv("data/lyva/lyva_cell.csv")

# Load SST data ----

oisst <- raster(list.files("data/temperature/noaa_oisst", full.names = T), 
                varname = 'sst')

oisst <- rotate(oisst)

# Extract temperature data ----

sp.temp <- raster::extract(oisst, sp.data[,1:2])

hist(sp.temp)


# bootstrapping with 1000 replications
results <- boot(
        data = sp.temp,
        statistic = function(data, indices) {
                m <- mean(data[indices], na.rm = T)
                return(m)
        },
        R = 1000
)

boot.ci(results, type="bca")

plot(results)


# Get min and max

sp.min <- min(sp.temp, na.rm = T)
sp.max <- max(sp.temp, na.rm = T)

# Restrict the map to the temperature

sp.map <- calc(oisst, fun = function(x){
        x[x >= sp.max] <- 0
        x[x <= sp.min] <- 0
        x[x != 0] <- 1
        x
})

cutpts <- c(0,0.5,1)
(plt <- levelplot(subset(sp.map, 1), margin = F, at=cutpts, cuts=11,
                  pretty=TRUE, par.settings = mapTheme,
                  main="OI SST limits"))


# Load GEE data

sp.gee <- read.csv("data/temperature/gee_data/lyva_points_test.csv")

sp.gee$system.index <- seq(1:nrow(sp.gee))

#Remove last 2 columns
sp.gee <- sp.gee[,-((ncol(sp.gee)-1):ncol(sp.gee))]

sp.gee <- sp.gee %>% 
        pivot_longer(cols = 2:ncol(sp.gee), 
                     names_to = "date",
                     values_to = "sst")

sp.gee$date <- gsub(pattern = "[.]", replacement = "_", x = sp.gee$date)

sp.gee <- sp.gee %>% separate(col = date, into = c("year", "month", "day"),
                               sep = "_")

sp.gee$sst[sp.gee$sst == "No data"] <- NA

sp.gee <- sp.gee[!is.na(sp.gee$sst),]

sp.gee$sst <- as.numeric(as.character(sp.gee$sst))

sp.gee.m <- sp.gee %>% group_by(year, month) %>% summarise(sst_mean = mean(sst))

ggplot(sp.gee, aes(x = month, y = sst))+
        stat_summary(
                mapping = aes(x = month, y = sst),
                fun.ymin = min,
                fun.ymax = max,
                fun.y = median
        )



### Simple sdm based solely on temperature

bath <- raster("data/env/bath_layers/bath_2_100.tif")

sst.layer <- crop(oisst, bath)
bath <- aggregate(bath, 3)
sst.layer <- mask(sst.layer, bath)

library(sdm)

sp.pts <- sp.data
colnames(sp.pts) <- c("x", "y", "lyva")

sst.layer <- stack(sst.layer)

coordinates(sp.pts) <- ~x+y

sdm.data <- sdmData(lyva~Daily.sea.surface.temperature, 
                    train = sp.pts, predictors = sst.layer,
                    bg=list(n=1000,method='gRandom',remove=TRUE))


sdm.model <- sdm(lyva~.,data=sdm.data,methods=c('brt', 'rf'))

sdm.pred <- predict(sdm.model, sst.layer)


#Extract Mahalonobis values from the occ. points
occPredVals <- raster::extract(sdm.pred$id_1.sp_1.m_brt, sp.data[,1:2])

occPredVals <- occPredVals[!is.na(occPredVals)]

#Get the 10th percentile value
if(length(occPredVals) < 10){
        p10 <- floor(length(occPredVals) * 0.9)
} else {
        p10 <- ceiling(length(occPredVals) * 0.9)
}

#Get the threshold
thresh <- rev(sort(occPredVals))[p10]

binarized.brt <- sdm.pred$id_1.sp_1.m_brt
binarized.brt <- calc(binarized.brt, function(x){
        x[x < thresh] <- 0
        x[x >= thresh] <- 1
        x
})




#linear models

full.data <- cbind(sp.data[,1:2], extract(sst.layer, sp.data[,1:2]))
full.data <- full.data[!is.na(full.data$Daily.sea.surface.temperature),]

xy.total <- as.data.frame(xyFromCell(sst.layer, 1:ncell(sst.layer)))

colnames(xy.total) <- c('decimalLongitude', 'decimalLatitude')

dups <- duplicated(rbind(xy.total, sp.data[,1:2]))

xy.total <- xy.total[!dups,]

sst.data <- cbind(xy.total, extract(sst.layer, xy.total))

sst.data <- sst.data[!is.na(sst.data$Daily.sea.surface.temperature),]

rp <- sample_n(sst.data, 1000)

full.data$species <- 1

rp$species <- 0

full.data <- rbind(full.data, rp)

lm1 <- lm(species~Daily.sea.surface.temperature, data = full.data)                   

plot(lm1)
