#### Modelling of coral reef herbivorous invertebrates ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2021

### Sea-urchins data download and cleaning ###
# We also aggregate the data from bibliographic sources

# Data downloaded on 2021-10-19

# Load packages ----
library(tidyverse)
library(robis)
library(rgbif)
library(obistools)
library(fuzzySim)



# Download data | OBIS ----

# List species that will be downloaded
sp.list <- c("Echinometra lucunter",
             "Lytechinus variegatus",
             "Tripneustes ventricosus")

# Create a unique code for each species
sp.codes <- tolower(spCodes(sp.list, nchar.gen = 2, nchar.sp = 2))

# Download from OBIS
obis.data <- lapply(sp.list, occurrence)

lapply(obis.data, nrow) # Verify number of records

# Add codes of species to the list to easy handling
names(obis.data) <- sp.codes




# Download data | GBIF ----

# Download from GBIF

gbif.data <- list()

for (i in 1:3) {
        gbif <- occ_data(scientificName = sp.list[i], hasCoordinate = T, limit = 30000)
        gbif.data[[i]] <- as.data.frame(gbif$data)
        rm(gbif)
}

lapply(gbif.data, nrow) # Verify number of records

# Add codes of species to the list to easy handling
names(gbif.data) <- sp.codes




# Data cleaning | OBIS ----

# Define year, record type and area of work
years <- as.character(c(1950:2021))
record <- c("HumanObservation", "Occurrence", "PreservedSpecimen")
area <- c(-42.5, 42.5, -99, 22) #Lat 1, Lat 2, Long 1, Long 2

# Remove points that don't meet requirements
obis.data <- lapply(obis.data, function(x){
        x <- x %>%
                dplyr::select(-flags) %>%
                filter(year %in% years) %>%
                filter(basisOfRecord %in% record) %>%
                filter(decimalLatitude >= area[1] & 
                               decimalLatitude <= area[2]) %>%
                filter(decimalLongitude >= area[3] & 
                               decimalLongitude <= area[4])
                
})

lapply(obis.data, nrow)

# Remove points on land
obis.land <- lapply(obis.data, check_onland)

for (i in 1:3) {
        obis.data[[i]] <- obis.data[[i]] %>% 
                anti_join(obis.land[[i]])
}

# Check each dataset for source (this one is done individually)

# E. lucunter
unique(obis.data$eclu$institutionCode)
#We use this piece of code to get the institutions codes

obis.data$eclu <- obis.data$eclu %>% 
        filter(!is.na(institutionCode)) #Here we remove just those that don't 
                                        #have institution info. But we could
                                        #remove others just adding other
                                        #filters
                                

# L. variegatus
unique(obis.data$lyva$institutionCode)

obis.data$lyva <- obis.data$lyva %>% 
        filter(!is.na(institutionCode))


# T. ventricosus
unique(obis.data$trve$institutionCode)

obis.data$trve <- obis.data$trve %>% 
        filter(!is.na(institutionCode)) %>%
        filter(institutionCode != "Diveboard") #Remove Diveboard data

# Optional: plot each data to verify
plot_map_leaflet(obis.data$trve)




# Data cleaning | GBIF ----

# Define year, record type and area of work
years <- as.character(c(1950:2021))
record <- c("HUMAN_OBSERVATION", 
            "PRESERVED_SPECIMEN", 
            "MACHINE_OBSERVATION",
            "LIVING_SPECIMEN")
area <- c(-42.5, 42.5, -99, 22) #Lat 1, Lat 2, Long 1, Long 2

# Remove points that don't meet requirements
gbif.data <- lapply(gbif.data, function(x){
        x <- x %>%
                filter(year %in% years) %>%
                filter(basisOfRecord %in% record) %>%
                filter(decimalLatitude >= area[1] & 
                               decimalLatitude <= area[2]) %>%
                filter(decimalLongitude >= area[3] & 
                               decimalLongitude <= area[4]) %>%
                select(-networkKeys)
})

lapply(gbif.data, nrow)

# Remove points on land
gbif.land <- lapply(gbif.data, check_onland)

for (i in 1:3) {
        gbif.data[[i]] <- gbif.data[[i]] %>% 
                anti_join(gbif.land[[i]])
}

# Check each dataset for source (this one is done individually)

# E. lucunter
unique(gbif.data$eclu$institutionCode)

gbif.data$eclu <- gbif.data$eclu %>% 
        filter(!is.na(institutionCode)) %>%
        filter(!institutionCode %in% c("iNaturalist", "BioObs", "Diveboard",
                                       "naturgucker"))


# L. variegatus
unique(gbif.data$lyva$institutionCode)

gbif.data$lyva <- gbif.data$lyva %>% 
        filter(!is.na(institutionCode)) %>%
        filter(!institutionCode %in% c("iNaturalist", "BioObs"))


# T. ventricosus
unique(gbif.data$trve$institutionCode)

gbif.data$trve <- gbif.data$trve %>% 
        filter(!is.na(institutionCode)) %>%
        filter(!institutionCode %in% c("iNaturalist", "BioObs", "Diveboard",
                                       "naturgucker", "SEAK", "NO DISPONIBLE",
                                       "INCONNU"))

# Optional: plot each data to verify
plot_map_leaflet(gbif.data$trve)

## End of data cleaning 




# Save full datasets ----
for (i in 1:3) {
        
        if (dir.exists(paste0("data/", sp.codes[i])) == F) {
                dir.create(paste0("data/", sp.codes[i]), recursive = T)
        }
        
        
        write.csv(obis.data[[i]], paste0("data/",
                                         sp.codes[i],
                                         "/", sp.codes[i],
                                         "_obis_fulldata.csv"),
                  row.names = F)
        
        write.csv(gbif.data[[i]], paste0("data/",
                                         sp.codes[i],
                                         "/", sp.codes[i],
                                         "_gbif_fulldata.csv"),
                  row.names = F)
}





# Bind data from bibliographic sources and databases ----

# Open data from bibliographic sources
biblio.data <- lapply(paste0("data/", sp.codes, "/", sp.codes, "_biblio.csv"), read.csv)

# Add a column with the species name (acronym)
for (i in 1:3) {
        biblio.data[[i]]$scientificName <- sp.codes[i]
}

# Load Scielo data (bibliographic source too)
# NOTE: the file is only on the ECLU folder, but contains data for all species
scielo.data <- read.csv("data/eclu/scielo_locations.csv")

scielo.data <- lapply(sp.codes, function(x){filter(scielo.data, Species == x) %>%
        mutate(scientificName = Species)})

# Select relevant collumns
biblio.data <- lapply(biblio.data, function(x){
        x <- x %>% select(scientificName, decimalLongitude, decimalLatitude)
})

scielo.data <- lapply(scielo.data, function(x){
        x <- x %>% select(scientificName, decimalLongitude, decimalLatitude)
})

# Select relevant collumns of OBIS and GBIF data
obis.data <- lapply(obis.data, function(x){
        x <- x %>% select(scientificName, decimalLongitude, decimalLatitude)
})

gbif.data <- lapply(gbif.data, function(x){
        x <- x %>% select(scientificName, decimalLongitude, decimalLatitude)
})

# Bind datasets
final.data <- list() # Create an empty list to bind datasets

#Change the names of OBIS and GBIF data to the acronyms and bind all data
for (i in 1:3) {
        obis.data[[i]]$scientificName <- sp.codes[i]
        gbif.data[[i]]$scientificName <- sp.codes[i]
        final.data[[i]] <- rbind(obis.data[[i]],
                                 gbif.data[[i]],
                                 biblio.data[[i]],
                                 scielo.data[[i]])
}

# Remove points on land (we do this again because we added the bibliographic
# data).
final.land <- lapply(final.data, check_onland)

for (i in 1:3) {
        final.data[[i]] <- final.data[[i]] %>% 
                anti_join(final.land[[i]])
}

sapply(final.land, nrow);sapply(final.data, nrow)


# Save final datasets ----
for (i in 1:3) {
        
        write.csv(final.data[[i]], paste0("data/",
                                         sp.codes[i],
                                         "/", sp.codes[i],
                                         "_final.csv"),
                  row.names = F)
}

#END of code