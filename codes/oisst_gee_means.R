#### Modelling of coral reef herbivorous invertebrates ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2021

### Extract max and min temperature experience by sea urchins ###
# Here we use daily data obtained from GEE (extract for each presence point)
# to establish three main data:
# 1) mean +- SD temperature experienced in the coldest presence point
# 2) mean +- SD temperature experienced in the hottest presence point
# 3) time that those temperatures are experienced (hot+SD; cold-SD) during the
# time span being used here.

# Whit this data we can project areas that are suitable (according to this time
# threshold) now and in the future.

# Load libraries ----
library(tidyverse)

# Create a function that can be applied in all species ----
getThresholds <- function(species){
        
        # Load data from GEE. For each presence point we got daily SST data
        # from the Optimum Interpolation, from NOAA.
        tdata <- read_csv(paste0("data/oisst/",
                                 species,
                                 "_daily_oisst.csv"))
        
        # GEE exports data in an awkward format. We need to do some 
        # data-wrangling here...
        tdata <- tdata %>% separate(`system:index`, into = c("date", "code"),
                                    sep = "_", extra = "drop") %>%
                # This is to extract the locations, not used here, but can be needed
                extract(.geo, into = c("locations"),
                        regex = "(\\[.[[:digit:]]*\\.[[:digit:]]*,[[:digit:]]*\\.[[:digit:]]*\\])") %>%
                separate(locations, into = c(NA, "locations"), sep = "\\[") %>%
                separate(locations, into = c("locations", NA), sep = "\\]") %>%
                separate(locations, 
                         into = c("decimalLongitude", "decimalLatitude"),
                         sep = ",") %>%
                # This is to separate dates in columns
                separate("date", into = c("year", "month_day"), sep = 4) %>%
                separate("month_day", into = c("month", "day"), sep = 2)
        
        # We remove years 2021 and 1981 for which data is incomplete
        #tdata <- filter(tdata, year < 2021 & year > 1981)
        
        # Correct the code column
        tdata$code <- as.factor(tdata$code)
        levels(tdata$code) <- seq(1:length(unique(tdata$code)))
        
        # Correct the temperature using the offset from the original file
        tdata$sst <- tdata$sst * 0.01
        
        # Get mean and SD for each point
        tdata.mean <- tdata %>% group_by(code) %>%
                summarise(mean_sst = mean(sst), sd = sd(sst))
        
        # Get the max and min values for each point and then the mean
        tdata.mm <- tdata %>% group_by(year, month, code) %>%
                summarise(max_sst = max(sst), min_sst = min(sst)) %>%
                group_by(code) %>%
                summarise(mean_max = mean(max_sst), mean_min = mean(min_sst))
        
        # Extract data for just the points where there is the maximum and
        # minimum mean temperatures (hottest and coolest points)
        # when more points are returned we just get the first
        max.sp <- tdata.mm[which(tdata.mm$mean_max == max(tdata.mm$mean_max)),][1,]
        min.sp <- tdata.mm[which(tdata.mm$mean_min == min(tdata.mm$mean_min)),][1,]
        
        # Then get the mean temperature + SD for the hottest and
        # mean - SD for the coolest (these are our thresholds)
        max.limit <- tdata.mean$mean_sst[max.sp$code]+tdata.mean$sd[max.sp$code]
        min.limit <- tdata.mean$mean_sst[min.sp$code]-tdata.mean$sd[min.sp$code]
        
        # Classify daily temperatures according to threshold
        # 1 - inside the range, 0 - outside the range
        tdata.range <- tdata %>%
                mutate(sst,
                       ran = ifelse(sst > max.limit | sst < min.limit,
                                    0,
                                    1))
        
        # Get the total number of days
        ntotal <- tdata %>% group_by(code) %>% summarise(total = n())
        unique(ntotal$total)
        ntotal <- unique(ntotal$total)
        
        # Sum the thresholded daily temperature and divide by total
        # number of days (thus getting the percentage of time the 
        # point was inside the temperature)
        tdata.range <- tdata.range %>% group_by(code) %>%
                summarise(summed = sum(ran)) %>%
                mutate(percent_time = summed/ntotal)
        
        # Get the percent of time that temperature was inside the range for the
        # hottest and coolest points
        max.time <- tdata.range$percent_time[max.sp$code]
        min.time <- tdata.range$percent_time[min.sp$code]
        
        # Export all data as a list
        all.data <- list("maximum_species" = max.sp,
                         "minimum_species" = min.sp,
                         "maximum_limit" = max.limit,
                         "minimum_limit" = min.limit,
                         "total_cases" = ntotal,
                         "time_per_point" = tdata.range,
                         "time_inrange_hottest_point" = max.time,
                         "time_inrange_coolest_point" = min.time)
        all.data
}

# Get thresholds for all species -----
# This may take sometime
thresholds <- lapply(c("eclu", "lyva", "trve"), getThresholds)

# Change names to easy handling
names(thresholds) <- c("eclu", "lyva", "trve")

# Access individual data
thresholds$eclu

# Save as RDS file with all data
saveRDS(thresholds, "data/oisst/allspecies_oisst_thvalues.rds")

### END