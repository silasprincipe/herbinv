#### Modelling of Herbivorous invertebrates of coral reefs ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2020

### Kernel density maps generation and occurrence thinning - Sea-urchins ###


# Kernel density maps generation ----

#### Creates the kernel density maps for the thinning procedure #####
###Note: thinning was done using the "OccurrenceThinner v1.04" from
###Dr. Heroen Verbruggen /// https://github.com/hverbruggen/OccurrenceThinner

# NOTE: this code may take a while to run, depending on your PC

# Load libraries
library(KernSmooth)
library(raster)
library(tidyverse)

# Define species acronyms
grp <- c("lyva", "eclu", "trve")

# Run in loop
for (i in 1:length(grp)) {
        
        # Load species occurrence data
        pts <- read.csv(paste0("data/", grp[i], "/", grp[i], "_cell.csv"))
        pts <- pts[pts[,3]==1,1:2]
        pts <- cbind(data.frame(species = rep(grp[i], nrow(pts))),
                     pts)
        
        # Generate the density map
        dens <- bkde2D(pts[,2:3], 
                       bandwidth=c(3,3), 
                       gridsize=c(4320,2160),
                       range.x=list(c(-180,180),c(-90,90)))
        
        # Create a raster to hold results
        dens.r <- raster(list(x = dens$x1, y = dens$x2, z = dens$fhat))
        projection(dens.r) <- CRS("+init=epsg:4326")
        xmin(dens.r) <- -180
        xmax(dens.r) <- 180
        ymin(dens.r) <- -90
        ymax(dens.r) <- 90
        
        # Write the final raster
        writeRaster(dens.r,
                    paste0("data/thinning/", grp[i], "_density.asc"),
                    "ascii",
                    overwrite = T)
        
        # Create a csv file with occurrences on the format used by Occ.Thinner
        colnames(pts) <- c("species", "longitude", "latitude")
        write.csv(pts, paste0("data/thinning/", grp[i], "_occ.csv"), 
                  row.names = F, quote = F)
        
        cat(grp[i], "density map done. \n")
        
        ###END
}

# Thinning ----

###Note: this is just a wrapper to run "OccurrenceThinner v1.04" directly
###from R. This code also adjust the files to the format used in the analysis.
###OccurrenceThinner is a java based program which was created by
###Dr. Heroen Verbruggen and is freely available at 
### https://github.com/hverbruggen/OccurrenceThinner

#Note: Thinning by OccurrenceThinner is random. Thus, running this code
#will certainly produce final files different from ours. Backuping to compare
#is advised. Anyway, final results should not be considerably different.

for(i in 1:length(grp)){
        
        #Define the working directory
        setwd("data/thinning/")
        
        #Create the running code for the OccurrenceThinner java prog.
        java.run <- paste0("java -jar OccurrenceThinner104.jar -i ", 
                           grp[i], "_occ.csv -r ", grp[i],
                           "_density.asc -nr 1 -t1 0.5 -t2 1")
        #These are the thresholds/options
        
        #Run the OccurrenceThinner java program
        java.out <- system(java.run, intern = TRUE)
        
        #Save output
        write.table(java.out, paste0(grp[i], "_javaoutput.txt"))
        
        #Rename file produced by OccurrenceThinner
        file.rename(paste0(grp[i], "_occ.csv.thinned1"), 
                    paste0(grp[i], "_occ_thinned.csv"))
        
        #Open the file produced
        pts.t <- read.csv(paste0(grp[i], "_occ_thinned.csv"))
        
        #Create a 'presence' column
        pts.t$pdata <- 1
        
        pts.t <- pts.t[,-1]
        
        #Change column names
        colnames(pts.t) <- c("decimalLongitude", "decimalLatitude", grp[i])
        
        #Return to the original directory
        setwd("../../")
        
        #Open the original points
        pts.o <- read.csv(paste0("data/", grp[i], "/", grp[i], "_cell.csv"))
        
        #Save an object with the number of points (P/A) just to register
        pts.dif <- nrow(pts.o)
        
        #Select only absences (if any)
        pts.o <- pts.o[pts.o[,3] == 0,]
        
        #Bind with the new presences selected by the thinning procedure
        pts.final <- rbind(pts.o, pts.t)
        
        #Get the difference in number of points
        pts.dif <- pts.dif - nrow(pts.final)
        
        #Save the final file
        write.csv(pts.final,
                  paste0("data/", grp[i], "/", grp[i], "_cell_thinned.csv"),
                  row.names = F)
        
        #Print results
        cat(grp[i], "thinning done.", pts.dif, "points were excluded. \n")
        
        #END
}
