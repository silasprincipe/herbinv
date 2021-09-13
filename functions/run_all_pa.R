run.all.pa <- function(species, model.list, nb.rep, env){
        
        results <- list()
        
        for (i in 1:length(model.list)) {
                results[[i]] <- list(
                        data.frame(matrix(ncol = 5, nrow = nb.rep)),
                        data.frame(matrix(ncol = 5, nrow = nb.rep)))
                
                names(results[[i]]) <- c("block_cv", "lat_cv")
                
                colnames(results[[i]][[1]]) <- c("AUC", "Spec",
                                                 "Sens", "TSS", "Boyce")
                
                colnames(results[[i]][[2]]) <- c("AUC", "Spec",
                                                 "Sens", "TSS", "Boyce")
                
                results[[i]][[3]] <- list()
        }
        
        names(results) <- model.list
        
        for (i in 1:nb.rep) {
                # Load species data
                pts <- read.csv(paste0("data/", species, "/", species, "_pts_", i,".csv"))
                
                dsp <- read.csv(paste0("data/", species, "/", species, "_dsp_", i,".csv"))
                dsp[is.na(dsp$dsp),] <- 0 # this data is in the same format to use with
                # BIOMOD, where PA is treated as NA
                
                # Load split table (block cross validation)
                blocks <- read.csv(paste0("data/", species, "/",
                                          species, "_foldindex_", i,".csv"))
                blocks <- blocks[,1]
                
                # The same for latitudinal blocks
                blocks.lat <- read.csv(paste0("data/", species, "/",
                                              species, "_foldindex_lat_",i,".csv"))
                blocks.lat <- blocks.lat[,1]
                
                # Extract the data for the points
                data <- raster::extract(env, pts)
                
                # Bind with prsence/absence data
                data <- cbind(dsp, data)
                
                # Creates a weigthed response
                wt <- function(resp, prev = 0.5){
                        
                        pres <- sum(resp)
                        
                        absc <- length(resp) - pres
                        
                        wgth <- rep(1, length(resp))
                        
                        wgth[which(resp > 0)] <- (prev * absc)/(pres * (1 - prev))
                        
                        wgth <- round(wgth)
                        
                        wgth
                }
                
                wt.dsp <- wt(dsp$dsp)
                
                # Generate a list with all data
                full <- list(
                        "points" = pts,
                        "pa" = dsp,
                        "blocks" = blocks,
                        "lat_blocks" = blocks.lat,
                        "data" = data,
                        "weigths" = wt.dsp
                )
                
                sp.data <- full
                
                
                #####
                if ("gbm" %in% model.list) {
                        
                        gbm.m <- gbm.modeling(sp.data$data,
                                              sp.data$blocks,
                                              sp.data$lat_blocks,
                                              sp.data$weigths)
                        
                        results[["gbm"]][[1]][i,] <- apply(gbm.m$block_cv,
                                                       2, mean, na.rm = T)
                        
                        results[["gbm"]][[2]][i,] <- apply(gbm.m$latitudinal_cv,
                                                       2, mean, na.rm = T)
                        
                        #cat("GBM finished.")
                        
                }
                
                if ("rf" %in% model.list) {
                        
                        rf.m <- rf.modeling(data = sp.data$data,
                                            sp.data$blocks,
                                            sp.data$lat_blocks)
                        
                        results[["rf"]][[1]][i,] <- apply(rf.m$block_cv,
                                                       2, mean, na.rm = T)
                        
                        results[["rf"]][[2]][i,] <- apply(rf.m$latitudinal_cv,
                                                       2, mean, na.rm = T)
                        
                }
                
                if ("glm" %in% model.list) {
                        
                        glm.m <- glm.modeling(" dsp ~ 1 + BO21_tempmean_ss +
                I(BO21_tempmean_ss^2) + BO21_salinitymean_ss + 
                BO21_lightbotmax_bdmean + I(BO21_lightbotmax_bdmean^2) + 
                BO_ph + I(BO_ph^2) + BO21_chlomean_ss + I(BO21_chlomean_ss^2) + 
                BO21_silicatemean_ss + I(BO21_silicatemean_ss^2) +
                BO21_dissoxmean_ss + I(BO21_dissoxmean_ss^2)",
                                              sp.data$data,
                                              sp.data$blocks,
                                              sp.data$lat_blocks,
                                              sp.data$weigths)
                        
                        results[["glm"]][[1]][i,] <- apply(glm.m$block_cv,
                                                          2, mean, na.rm = T)
                        
                        results[["glm"]][[2]][i,] <- apply(glm.m$latitudinal_cv,
                                                          2, mean, na.rm = T)
                        
                        pred <- predict(env, glm.m$model,
                                                              type = "response")
                        
                        pred[pred < glm.m$tss_thresh] <- 0
                        pred[pred >= glm.m$tss_thresh] <- 1
                        
                        results[["glm"]][[3]][[i]] <- pred
                        
                }
                
                
                if ("gam" %in% model.list) {
                        
                        gam.m <- gam.modeling("dsp ~ 1 + s(BO21_tempmean_ss, k = 4) + s(BO21_salinitymean_ss, 
    k = 4) + s(BO21_lightbotmax_bdmean, k = 4) + s(BO_ph, k = 4) + 
    s(BO21_chlomean_ss, k = 4) + s(BO21_silicatemean_ss, k = 4) + 
    s(BO21_dissoxmean_ss, k = 4)",
                                              sp.data$data,
                                              sp.data$blocks,
                                              sp.data$lat_blocks,
                                              sp.data$weigths)
                        
                        results[["gam"]][[1]][i,] <- apply(gam.m$block_cv,
                                                           2, mean, na.rm = T)
                        
                        results[["gam"]][[2]][i,] <- apply(gam.m$latitudinal_cv,
                                                           2, mean, na.rm = T)
                        
                }
                
                cat("Run", i, "done. \n")
                
        }
        
        results
        
}
