# Based on heatwaveR vignette 
# (https://cran.r-project.org/web/packages/heatwaveR/vignettes/OISST_preparation.html)


# The packages we will use
library(dplyr) # A staple for modern data management in R
library(lubridate) # Useful functions for dealing with dates
library(ggplot2) # The preferred library for data visualisation
library(tidync) # For easily dealing with NetCDF data
library(rerddap) # For easily downloading subsets of data
library(doParallel) # For parallel processing

# The information for the NOAA OISST data
rerddap::info(datasetid = "ncdcOisst21Agg_LonPM180", url = "https://coastwatch.pfeg.noaa.gov/erddap/")

# This function downloads and prepares data based on user provided start and end dates
OISST_sub_dl <- function(time_df){
        OISST_dat <- griddap(x = "ncdcOisst21Agg_LonPM180", 
                             url = "https://coastwatch.pfeg.noaa.gov/erddap/", 
                             time = c(time_df$start, time_df$end), 
                             zlev = c(0, 0),
                             latitude = c(-42.5, 42.5),
                             longitude = c(-99, -29),
                             fields = "sst")$data %>% 
                mutate(time = as.Date(stringr::str_remove(time, "T00:00:00Z"))) %>% 
                dplyr::rename(t = time, temp = sst) %>% 
                select(lon, lat, t, temp) %>% 
                na.omit()
}

dl_years <- data.frame(date_index = 1,
                       start = as.Date(c("2018-01-01")),
                       end = as.Date(c("2018-12-31")))

# Download all of the data with one nested request
# The time this takes will vary greatly based on connection speed
system.time(
        OISST_data <- dl_years %>% 
                group_by(date_index) %>% 
                group_modify(~OISST_sub_dl(.x)) %>% 
                ungroup() %>% 
                select(lon, lat, t, temp)
) # 636 seconds, ~127 seconds per batch

saveRDS(data2019, "oisst_2019.Rds")

OISST_data %>% 
        filter(t == "2019-12-01") %>% 
        ggplot(aes(x = lon, y = lat)) +
        geom_tile(aes(fill = temp)) +
        #borders() + # Activate this line to see the global map
        scale_fill_viridis_c() +
        coord_quickmap(expand = T) +
        labs(x = NULL, y = NULL, fill = "SST (°C)") +
        theme(legend.position = "bottom")
