###### Reef-builders modelling #######
#Silas C. Principe - 2020
#silasprincipe@yahoo.com.br

#### Applying change factor to future layers #####

####Read first:
#To create the future layers we used the change-factor approach,
#the same used by the Bio-ORACLE team. In this approach, coarse layers from
#CMIP are used to obtain min,max or mean from a reference condition
#(in our case 2000-2014) and a future condition (2091-2100)
#for each scenario of RCP. Then we get the absolute difference between the reference
#condition and the new one. All this procedure was done using CDO (Climate Data Operators).
#More details of this software (that is opensource) can be found on
#https://code.mpimet.mpg.de/projects/cdo/
#Codes for this step are supplied (embedded on R using "system", but CDO is needed)
#Then the change factor layer was interpolated to the same resolution of our climate layers
#using kriging interpolation in R (code is on the file interpolation.r)
#Now we will apply this change factor to our climate layers from Bio-ORACLE.
####

# Load libraries ----
library(raster)
library(stringr)

# Apply change factor ----

# Establish the scenarios
scenarios <- c("ssp126","ssp245","ssp585")

# Establish which variables will be used
var.bior <- read.table("data/env/env_layers.txt")[,2]
var.bior <- var.bior[-3] # We remove lightbottom, no future projection

# Names are different for CMIP files. Have to be in the same order!
var.cmip <- c("tos",
              "sos",
              "phos",
              "chlos",
              "sios",
              "o2os")

# Establish the type of cf to apply (relative or absolute)
cft <- c("ac",
         "rc",
         "ac",
         "rc",
         "ac",
         "ac")

# Load bath layer
bath <- raster("data/env/bath_layers/bath_2_100.tif")

# Apply CF in loop
for (i in 1:length(scenarios)) {

  # Establish the scenario
  set <- scenarios[i]
  
  # Create dir if still not created
  if (!dir.exists(paste0("data/env/proj_layers/", set))) {
          dir.create(paste0("data/env/proj_layers/", set)) 
  }
  
  cat("Files for scenario", set, "\n")
  
  # Apply for each variable
  for (z in 1:length(var.cmip)) {
          
          # Open the CMIP file
          cfact <- stack(paste0("data/cmip/", var.cmip[z],"/",
                                paste(var.cmip[z], "gfdl",
                                      set,
                                      # Here we use conditionals to see if mean
                                      # max, or min
                                      # If the variable is ph, it's mean
                                      ifelse(var.cmip[z] == "phos",
                                             "mean",
                                             ifelse(str_detect(var.bior[z],
                                                               "mean"),
                                                    "mean",
                                                    ifelse(str_detect(var.bior[z],
                                                                      "max"),
                                                           "max",
                                                           "min"))),
                                      cft[z],
                                      sep = "_"),
                                ".tif"))
          
          # Mask with the bath layer
          cfact <- mask(cfact, bath)
          
          # Open the Bio-ORACLE file
          orig <- raster(paste0("data/env/crop_layers/",
                                var.bior[z],
                                ".tif"))
          
          # If absolute change
          if (cft[z] == "ac") {
                  
                  cat("Combining", var.bior[z],
                      var.cmip[z], "using absolute change", "\n")
                  
                  # Sum layers
                  final <- orig + cfact
                  
                  # Write final raster
                  writeRaster(final,
                              paste0("data/env/proj_layers/",set,"/",
                                     names(orig),".tif"),
                              overwrite = T)
                  
          } 
          
          if (cft[z] == "rc") {
                  
                  cat("Combining", var.bior[z],
                      var.cmip[z], "using relative change", "\n")
                  
                  # Get relative change factor then sum layers
                  relc <- orig * cfact
                  
                  final <- orig + relc
                  
                  # Write final raster
                  writeRaster(final,
                              paste0("data/env/proj_layers/",set,"/",
                                     names(orig),".tif"),
                              overwrite = T)
                  
          } 
  }
}

##END##
