library(raster)
library(fs)

layers <- data.frame(
        bior = c("ph", "chlomax", "chlomean", "dissoxmean", 
                 "salinitymean", "silicatemax", "silicatemean",
                 "tempmax", "tempmean"),
        cmip = c("phos", "chlos", "chlos", "o2os", "sos",
                 "sios", "sios", "tos", "tos"),
        var = c("mean", "max", "mean", "mean",  "mean", "max",
                "mean", "max", "mean")
)

dir_create("data/env/proj_layers/sd")
dir_create(paste0("data/env/proj_layers/",
                  c("ssp126", "ssp245", "ssp370", "ssp585")))

for (k in 1:4) {
        
        ssp <- c("ssp126", "ssp245", "ssp370", "ssp585")[k]
        
        for (i in 1:nrow(layers)) {
                bio.f <- list.files("data/env/crop_layers/",
                                    pattern = layers[i, 1],
                                    full.names = T)
                
                bio.f <- raster(bio.f)
                
                cmip.f <- list.files(paste0("data/cmip/", layers[i, 2], "/"),
                                     pattern = layers[i, 2])
                
                cmip.f <- cmip.f[grep(layers[i, 3], cmip.f)]
                
                cmip.f <- cmip.f[grep(ssp, cmip.f)]
                
                cmip.f <- paste0("data/cmip/", layers[i, 2], "/", cmip.f)
                
                cmip.f <- stack(cmip.f)
                
                cmip.f <- extend(cmip.f, extent(bio.f))
                
                rel.change <- bio.f * cmip.f
                
                change <- bio.f + rel.change
                
                change.m <- calc(change, mean)
                
                change.sd <- calc(change, sd)
                
                delta <- change.m - bio.f
                
                par(mfrow = c(1,2))
                plot(change.m, main = layers[i, 2])
                plot(delta)
                
                print(summary(change.m))
                print(summary(bio.f))
                print(summary(delta))
                
                cat("Max/min change mean")
                print(cellStats(change.m, max))
                print(cellStats(change.m, min))
                cat("Max/min original")
                print(cellStats(bio.f, max))
                print(cellStats(bio.f, min))
                
                
                writeRaster(change.m,
                            paste0("data/env/proj_layers/",
                                   ssp, "/",
                                   names(bio.f), ".tif"),
                            overwrite = T)
                
                writeRaster(change.sd,
                            paste0("data/env/proj_layers/sd/",
                                   names(bio.f), "_", ssp,
                                   "_sd.tif"),
                            overwrite = T)
        }
}
