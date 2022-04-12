##### Get integration points weights for LGCP
# In the case of inlabru, this is just a wrapper for the ipoints function,
# but that saves the ipoints, so the initialization is faster when running
# again.

# MODES AVAILABLE:
# ips - integration points returned as SpatialPoints (for inlabru use)
# simple - integration points calculated using inlabru ipoints, but returned
# in the format for lgcp through plan INLA code.
# other - integration points weights calculated using dual mesh from Elias K.

get.weights <- function(mesh, starea, mode = "ips", save = T, folder = "data/"){
        
        lf <- list.files(folder)
        
        if ("w.txt" %in% lf) {
                cat("Saved weights vector found. Loading it. \n")
                w <- read.table(paste0(folder, "w.txt"), header = T)
                
                if (mode == "ips") {
                        if (length(w) <= 1) {
                                stop("Previous weight file was of different type. Delete it and redo.\n")
                        }else{
                                crs.w <- crs(starea)
                                w <- w[,1:3]
                                w <- SpatialPointsDataFrame(w[,1:2], data.frame(weight = w[,3]),
                                                            proj4string = crs.w)
                        }
                        
                }else{
                        w <- w[,1]
                }
                
        }else{
                if (mode == "simple") {
                        cat("Generating weights through the inlabru mode ---- \n")
                        
                        vlocs <- SpatialPoints(mesh$loc[,1:2])
                        crs(vlocs) <- crs(starea)
                        
                        ips <- ipoints(starea, mesh)
                        
                        ipb <- buffer(ips, 1)
                        
                        ov <- over(vlocs, ipb)
                        ov2 <- ifelse(is.na(ov), FALSE, TRUE)
                        
                        tloc <- data.frame(vlocs@coords)
                        tloc$w <- 0
                        tloc$w[ov2] <- ips$weight
                        
                        w <- tloc$w
                        
                        cat("\tNOTE: sum of weights approximately the total area.
                    To get equal area, use mode 'dmesh'. \n")
                }
                
                if (mode == "ips"){
                        cat("Generating weights for inlabru use ---- \n")
                        
                        w <- ipoints(starea, mesh)
                        
                } else {
                        cat("Generating weights through the dual mesh mode ---- \n")
                        cat("Function created by Elias Krainski \nhttps://becarioprecario.bitbucket.io/spde-gitbook/\n")
                        
                        book.mesh.dual <- function(mesh) {
                                if (mesh$manifold=='R2') {
                                        ce <- t(sapply(1:nrow(mesh$graph$tv), function(i)
                                                colMeans(mesh$loc[mesh$graph$tv[i, ], 1:2])))
                                        library(parallel)
                                        pls <- mclapply(1:mesh$n, function(i) {
                                                p <- unique(Reduce('rbind', lapply(1:3, function(k) {
                                                        j <- which(mesh$graph$tv[,k]==i)
                                                        if (length(j)>0) 
                                                                return(rbind(ce[j, , drop=FALSE],
                                                                             cbind(mesh$loc[mesh$graph$tv[j, k], 1] +
                                                                                           mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 1], 
                                                                                   mesh$loc[mesh$graph$tv[j, k], 2] +
                                                                                           mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 2])/2))
                                                        else return(ce[j, , drop=FALSE])
                                                })))
                                                j1 <- which(mesh$segm$bnd$idx[,1]==i)
                                                j2 <- which(mesh$segm$bnd$idx[,2]==i)
                                                if ((length(j1)>0) | (length(j2)>0)) {
                                                        p <- unique(rbind(mesh$loc[i, 1:2], p,
                                                                          mesh$loc[mesh$segm$bnd$idx[j1, 1], 1:2]/2 +
                                                                                  mesh$loc[mesh$segm$bnd$idx[j1, 2], 1:2]/2, 
                                                                          mesh$loc[mesh$segm$bnd$idx[j2, 1], 1:2]/2 +
                                                                                  mesh$loc[mesh$segm$bnd$idx[j2, 2], 1:2]/2))
                                                        yy <- p[,2]-mean(p[,2])/2-mesh$loc[i, 2]/2
                                                        xx <- p[,1]-mean(p[,1])/2-mesh$loc[i, 1]/2
                                                }
                                                else {
                                                        yy <- p[,2]-mesh$loc[i, 2]
                                                        xx <- p[,1]-mesh$loc[i, 1]
                                                }
                                                Polygon(p[order(atan2(yy,xx)), ])
                                        })
                                        return(SpatialPolygons(lapply(1:mesh$n, function(i)
                                                Polygons(list(pls[[i]]), i))))
                                }
                                else stop("It only works for R2!")
                        }
                        
                        dmesh <- book.mesh.dual(mesh)
                        
                        # Run in parallel (still taking a lot of time, but much better!)
                        # library(parallel)
                        # cl <- makeCluster(3)
                        # clusterEvalQ(cl, library(rgeos))
                        # clusterExport(cl, "starea");clusterExport(cl, "dmesh")
                        # w <- parSapply(cl, 1:length(dmesh), function(i) {
                        #         if (gIntersects(dmesh[i, ], starea))
                        #                 return(gArea(gIntersection(dmesh[i, ], starea)))
                        #         else return(0)
                        # })
                        # stopCluster(cl)
                        w <- sapply(1:length(dmesh), function(i) {
                                if (gIntersects(dmesh[i, ], starea))
                                        return(gArea(gIntersection(dmesh[i, ], starea)))
                                else return(0)
                        })
                }
                
                if (save) {
                        cat("Weights vector saved in folder '", folder, "'\n")
                        
                        if (mode == "ips") {
                                w2 <- cbind(w@coords, w@data)
                                write.table(w2, paste0(folder, "w.txt"), row.names = F)
                        }else{
                                write.table(w, paste0(folder, "w.txt"), row.names = F)
                        }
                }
        }
        
        return(w)
}