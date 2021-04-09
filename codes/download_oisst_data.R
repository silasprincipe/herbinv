# Download OISST data

library(rerddap)
library(fs)

dates.df <- data.frame("init" = paste0(seq(2010, 2020), "-01-01"),
                       "final" = paste0(seq(2010, 2020), "-12-31"))

rlist <- list()

for (i in 1:11) {
        
        rlist[[i]] <- griddap(x = "ncdcOisst21Agg_LonPM180", 
                              url = "https://coastwatch.pfeg.noaa.gov/erddap/", 
                              time = c(dates.df$init[i], dates.df$final[i]), 
                              zlev = c(0, 0),
                              latitude = c(-40, 40),
                              longitude = c(-101, -27),
                              fields = "sst",
                              read = F)
}


for (i in 1:11) {
        
        file_move(rlist[[i]]$summary$filename,
                  paste0("data/oisst/daily_oisst_noaa_",
                         seq(2010, 2020)[i],
                         ".nc"))
        
}
