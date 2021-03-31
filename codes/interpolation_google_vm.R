library(raster)
library(automap)
library(gstat)
library(dplyr)
library(snowfall)


#Creates a fuction for interpolation
interpol <- function(x, r.lay, outf, fnam) {
        
        #Get the name of the layer
        lname <- names(r.lay[[x]])
        
        #Get the layer from the list
        r <- r.lay[[x]]
        
        #Create a dataframe with data from raster
        tmp <-
                data.frame(cbind(xyFromCell(r, 1:length(r[[1]])), values(r[[1]])))
        
        #change colnames
        colnames(tmp) <- c('lon', 'lat', 'var')
        
        #Remove NAs
        tmp <- tmp[!is.na(tmp$var),]
        
        #Crop to the area of study
        tmp <-
                tmp %>% filter(lon >= -90 &
                                       lon <= 22) %>% filter(lat >= -42.5 &
                                                                     lat <= 42.5)
        
        #Convert to SpatialGrid
        coordinates(tmp) <- ~ lon + lat
        
        #Correct CRS
        crs(tmp) <- crs(base)
        crs(r) <- crs(base)
        
        #Create variograms
        #v <- variogram(var ~ 1, tmp)
        #Use automap to fit the variogram automatically
        vario <-
                automap::autofitVariogram(formula = var ~ 1, input_data = tmp)
        m <- vario$var_model
        
        #If wanted, variogram can be fitted mannualy
        #m <- fit.variogram(v, vgm(1, "Sph", 500, 1))
        
        #Create the kriging object, using 12 nearest cells
        gform <- gstat(NULL, "var", var ~ 1, tmp, model = m, nmax = 12)
        
        #Interpolate using the basemap as model
        imap <- interpolate(base, gform)
        
        #Correct raster name
        names(imap) <- lname
        
        #Mask with the base map
        imap <- mask(imap, base)
        
        #Return the object
        writeRaster(imap, paste0("data/",
                                 outf, "/",
                                 fnam, names(imap), ".tif"))
        
}


#Define study extent
ext <- extent(-99, -29, -42.5, 42.5)

base <- raster("data/temperature/noaa_oisst/oisst-avhrr-v02r01.20201101.nc")

base <- rotate(base)

base <- crop(base, ext)


for (i in 1:4) {
        
        file <- paste0("data/cmip/",
                       c("tos_hist_miroc6_20000101-20091231.nc",
                         "tos_hist_miroc6_20100101-20141231.nc",
                         "tos_ssp585_miroc6_20850101-20941231.nc",
                         "tos_ssp585_miroc6_20950101-21001231.nc")[i])
        
        tos <- brick(file)    
        
        tos <- crop(tos, ext)
        
        nt <- nlayers(tos)

        outfolder <- c("hist_1", "hist_2", "ssp585_1", "ssp585_2")[i]
        
        fname <- c("hist1_", "hist2_", "ssp1_", "ssp2_")[i]
        
        #initiate parallel computing
        sfInit(parallel = TRUE, cpus = 8) ##I have 4 cores, using 3 here
        
        ## Export packages
        sfLibrary('raster', character.only = TRUE)
        sfLibrary('automap', character.only = TRUE)
        sfLibrary('gstat', character.only = TRUE)
        sfLibrary('dplyr', character.only = TRUE)
        
        ## Export variables
        sfExport('tos')
        sfExport('base')
        sfExport('outfolder')
        sfExport('fname')
        
        cat("Starting interpolation. Running in parallel using 3 cores.",
            "\n")
        
        ## Interpolate
        sfLapply(1:nt, interpol, r.lay = tos, outf = outfolder, fnam = fname)
        
        ## stop snowfall
        sfStop(nostop = FALSE)
        
}
