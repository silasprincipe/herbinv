#### Modelling of Herbivorous invertebrates of coral reefs ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2021

### Modeling - Sea-urchins ###

# Load packages ----
library(raster)
library(snowfall)

# Set species
splist <- c("lyva", "eclu", "trve")


# Extract from current period

noaa <- list.files("data/sst_limits/noaa/", pattern = "\\.nc",
                   full.names = T)

noaa <- noaa[grep(paste0(2000:2014, collapse = "|"), noaa)]

noaa <- lapply(noaa, brick)

noaa <- stack(noaa)

noaa <- rotate(noaa)

ext.cur.temp <- function(species, sst){
        
        sp.data <- read.csv(paste0("data/", species,
                                   "/", species, "_cell_noaa.csv"))
        
        sp.data <- sp.data[sp.data[,3] == 1,]
        
        temps <- extract(sst, sp.data[,1:2])
        
        temps <- cbind(sp.data[, 1:2], temps)
        
        write.csv(temps, paste0("data/sst_limits/", species,
                                "_sst_noaa_current.csv"),
                  row.names = F)
        
        paste(species, "done.")
}

lapply(splist, ext.cur.temp, noaa)

# # Creates a function so it can be done in parallel
# extract.temp <- function(species, ssp){
#         
#         files <- list.files(paste0("data/sst_limits/cmip6/",
#                                    ssp, "/mean"))
#         
#         files <- files[grep(paste0(2086:2100, collapse = "|"), files)]
#         
#         sp.data <- read.csv(paste0("data/", species,
#                                    "/", species, "_cell_noaa.csv"))
#         
#         sp.data <- sp.data[sp.data[,3] == 1,]
#         
#         temps <- data.frame(matrix(ncol = length(files),
#                                    nrow = nrow(sp.data)))
#         
#         nam <- gsub("\\.", "_",
#                     gsub(".tif", "",
#                 gsub(paste0("tos_", ssp, "_X"), "", files)))
#         
#         colnames(temps) <- nam
#         
#         files <- paste0("data/sst_limits/cmip6/", ssp, "/mean/", files)
#         
#         for (i in 1:length(files)) {
#                 
#                 r <- raster(files[i])
#                 
#                 ext <- extract(r, sp.data[,1:2])
#                 
#                 temps[,i] <- ext
#         }
#         
#         temps <- cbind(sp.data[,1:2], temps)
#         
#         write.csv(temps, paste0("data/sst_limits/", species,
#                                 "_sst_noaa_", ssp, ".csv"),
#                   row.names = F)
#         
#         paste(ssp, species, "done.")
#         
# }
# 
# # Extract for each ssp
# 
# sfInit(parallel = T, cpus = 3)
# 
# sfLibrary(raster)
# 
# sfLapply(splist, extract.temp, "ssp126")
# 
# sfLapply(splist, extract.temp, "ssp245")
# 
# sfLapply(splist, extract.temp, "ssp370")
# 
# sfLapply(splist, extract.temp, "ssp585")
# 
# sfStop()
# 
