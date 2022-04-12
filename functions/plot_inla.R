# Plot random effects
plot.random <- function(model, CI=T, adpt = T, disable.par = F) {
        lr <- length(model$summary.random)
        if (lr-1 == 0) {
                lr <- 2
        }
        if (!disable.par) {
                if (lr-1 <= 3) {
                        par(mfrow = c(1,(lr-1)))
                } else {par(mfrow = c(2,3))}   
        }
        for (i in 1:(lr-1)) {
                
                nam <- gsub("\\, n = [[:alnum:]][[:alnum:]]", "",
                            gsub("[\\(\\)]", "",
                                 gsub("inla.group", "",
                                      names(model$summary.random))))
                
                if(CI){
                        ylim <- c(min(c(model$summary.random[[i]]$`0.025quant`)),
                                  max(c(model$summary.random[[i]]$`0.975quant`)))
                } else {ylim <- NULL}
                
                plot(c(model$summary.random[[i]]$mean) ~
                             model$summary.random[[i]]$ID, type = "l",
                     main = nam[i], ylim = ylim,
                     xlab = "", ylab = "Effect")
                abline(h = 0, lty = "dashed", lwd = 1, col = "gray70")
                
                if(CI){
                        lines(x=model$summary.random[[i]]$ID,
                              y=c(model$summary.random[[i]]$`0.025quant`),
                              lwd=2, lty="dashed", col="blue")
                        lines(x=model$summary.random[[i]]$ID,
                              y=c(model$summary.random[[i]]$`0.975quant`),
                              lwd=2, lty="dashed", col="blue")
                }
                if(adpt){
                        points(model$summary.random[[i]][,1:2], pch = 20)
                }
                
                
        }
}


# plot fixed effects
plot.fixed <- function(model){
        
        par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
        
        fef <- data.frame(
                eff_name = model$names.fixed,
                mean = c(model$summary.fixed[,"mean"]),
                upper = c(model$summary.fixed[,"0.975quant"]),
                lower = c(model$summary.fixed[,"0.025quant"])
        )
        
        plot(y = 1:nrow(fef), x = fef$mean, yaxt = "n", ylab = "", xlab = "Effect",
             xlim=c(min(fef$lower),max(fef$upper)), ylim=c(0.5,(nrow(fef)+0.5)), pch = 19, cex = 1)
        lines(rbind(fef$lower,fef$upper,NA), rbind(1:nrow(fef),1:nrow(fef),NA))
        abline(v = 0, lty = 2 )
        axis(side=2,at=1:nrow(fef),labels=fef$eff_name)
}

# plot posteriors
plot.post <- function(model){
        
        if (length(model$marginals.fixed) <= 3) {
                par(mfrow = c(1, length(model$marginals.fixed)))
        } else {
                par(mfrow = c(2,3))
        }
        
        for (i in 1:length(model$marginals.fixed)) {
                plot(model$marginals.fixed[[i]], type="l", main = model$names.fixed[i])
                abline(v=0, col="red", lwd=2)
        }
}


# Plot spatial effect
plot.spatial <- function(model, mesh, pts = NULL, addb = T){
        library(fields)
        # Get projection mesh
        projmesh <- inla.mesh.projector(mesh, xlim = starea@bbox[1, ], 
                                        ylim = starea@bbox[2, ], dims=c(300, 300))
        
        par(mfrow = c(1,2), mar = c(5, 5, 4, 5) + 0.1)
        mean.g <- inla.mesh.project(projmesh, model$summary.random$spatial$mean)
        sd.g <- inla.mesh.project(projmesh, model$summary.random$spatial$sd)
        
        # Plot mean spatial effect
        image.plot(projmesh$x, 
                   projmesh$y,
                   mean.g, axes=TRUE, xlab="Longitude", ylab="Latitude", main = "Spatial effect")
        if (addb) {
                plot(poly.barrier,add=T, axes=TRUE, col = "grey")
        }
        
        if (!is.null(pts)) {
                points(pts, pch = 16, cex = .5, col = "white")
        }
        
        # Plot SD of spatial effect
        image.plot(projmesh$x, 
                   projmesh$y,
                   sd.g, axes=TRUE, xlab="Longitude", ylab="Latitude", main = "SD")
        if (addb) {
                plot(poly.barrier,add=T, axes=TRUE, col = "grey")
        }
}


# Plot cpo or pit
plot.cpo <- function(modellist, mode = "cpo", n = 1000, pts = T, ptsn = 300){
        
        if (pts) {
                idx <- length(modellist[[1]]$cpo$cpo):(length(modellist[[1]]$cpo$cpo) - ptsn)
        }else{
                idx <- sample(1:length(modellist[[1]]$cpo$cpo), n)
        }
        
        ccol <- list()
        
        if (length(modellist) <= 4) {
                ccol[[1]] <- rgb(red = 0, green = 0, blue = 1, alpha = 0.2)
                ccol[[2]] <- rgb(red = 1, green = 0, blue = 0, alpha = 0.2)
                ccol[[3]] <- rgb(red = 0, green = 1, blue = 0.1, alpha = 0.2)
                ccol[[4]] <- rgb(red = 0.5, green = 0, blue = 1, alpha = 0.2)
        }else{
                for (i in 1:length(modellist)) {
                        ccol[[i]] <- rgb(red = sample(seq(0,1, by = 0.1), 1),
                                         green = sample(seq(0,1, by = 0.1), 1),
                                         blue = sample(seq(0,1, by = 0.1), 1),
                                         alpha = 0.2)
                }
        }
        
        for (i in 1:length(modellist)) {
                
                if (mode == "cpo") {
                        tab <- data.frame(point = 1:length(modellist[[i]]$cpo$cpo), 
                                          cpo = modellist[[i]]$cpo$cpo)
                }else{
                        tab <- data.frame(point = 1:length(modellist[[i]]$cpo$cpo), 
                                          cpo = modellist[[i]]$cpo$pit)
                }
                if (i == 1) {
                        plot(y = tab$cpo[idx], x = tab$point[idx], pch = 16,
                             col = ccol[[i]], xlab = "Point", ylab = toupper(mode))
                }else{
                        points(y = tab$cpo[idx], x = tab$point[idx], pch = 16,
                               col = ccol[[i]])
                }

        }
        
        legend("topleft", legend = paste("M", 1:length(modellist)), col = unlist(ccol),
               pch = 16)
}


# Plot ROC curves ----
plot.roc <- function(real, pred){
        plot(1 - real$specificities, real$sensitivities,
             xlim = c(0, 1), ylim = c(0, 1), col = "cyan3", type = "l",
             xlab = "1-Specificity", ylab = "Sensitivity", main = "ROC curve")
        lines(1 - pred$specificities, pred$sensitivities)
        legend("bottomright", legend = c("Fit", "CV"), col = c("cyan3", "black"),
               lty = c(1, 1))
}