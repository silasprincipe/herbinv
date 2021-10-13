#### Modelling of coral reef herbivorous invertebrates ####
## Silas C. Principe - silasprincipe@yahoo.com.br - 2021

### Extract max and min temperature experienced by sea urchins ###
library(tidyverse)
library(lubridate)

thresholds <- list()

for (i in 1:3) {
        
        species <- c("lyva", "eclu", "trve")[i]
        
        tdata <- read.csv(paste0("data/sst_limits/", species, "_sst_noaa_current.csv"))
        
        tdata <- tdata[,-c(1:2)]
        
        tdata$code <- 1:nrow(tdata)
        
        
        tdata <- pivot_longer(tdata, cols = 1:(length(tdata)-1), names_to = "date",
                              values_to = "sst")
        
        
        tdata$date <- as_date(tdata$date, format = "X%Y.%m.%d")
        
        
        tdata.mean <- tdata %>% group_by(code) %>%
                summarise(mean_sst = mean(sst), sd = sd(sst))
        
        tdata.mm <- tdata %>%
                group_by(code, year(date), month(date)) %>%
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
        
        thresholds[[i]] <- all.data
}


# Change names to easy handling
names(thresholds) <- c("lyva", "eclu", "trve")

# Access individual data
thresholds$eclu

# Save as RDS file with all data
save(thresholds, file = "data/sst_limits/allspecies_oisst_thvalues.RData")
