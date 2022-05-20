data("gorillas")
crs <- sp::CRS("+proj=utm +zone=32N +datum=WGS84") 

nests <- gorillas$nests

boundary <- gorillas$boundary

bnd <- INLA::inla.sp2segment(boundary)
mesh <- INLA::inla.mesh.2d(
        interior = bnd, max.edge = 222,
        #cutoff = (0.01*111),
        crs = crs
)

ggplot()+gg(mesh)+coord_equal();mesh$n

mesh <- gorillas$mesh
mesh$meta


# Matérn models ----
# Get range and sigma
size <- diff(range(mesh$loc[,1]))
range0 <- round(size / 2)
sigma <- 3


matern <- inla.spde2.pcmatern(mesh,
                              prior.sigma = c(10, 0.01),
                              prior.range = c(5, 0.01)
)

cmp <- coordinates ~ mySmooth(coordinates,
                              model = matern
) +
        Intercept(1)

#ips <- ipoints(boundary, mesh)

fit <- lgcp(cmp, nests, ips = ips, domain = list(coordinates = mesh),
            options = list(control.inla = list(int.strategy = "eb")))

pa <- c(rep(0, length(ips)), rep(1, length(nests)))
sd(pa)

Lambda <- predict(
        fit,
        ipoints(boundary, mesh),
        ~ sum(weight * exp(mySmooth + Intercept))
)
Lambda



# Specify spatial range and variance
sim_range <- 0.5
sim_sd <- 0.1

GRF_sim()


# Function for simulating and plotting a GRF on a unit square
GRF_sim <- function()
{
        # define inla mesh
        mesh <- inla.mesh.2d(loc.domain = cbind(c(0,1,1,0),c(0,0,1,1)),
                             min.angle = 21, max.edge = 0.04, cutoff = 0.02)
        # plot(mesh)
        # mesh$n
        
        # define spde
        spde_obj <- inla.spde2.pcmatern(mesh, constr = T,
                                        prior.range = c(sim_range, 0.5),
                                        prior.sigma = c(sim_sd, 0.5))
        
        # define indices
        index <- inla.spde.make.index('smooth',n.spde = spde_obj$n.spde)
        
        # define prediction matrix
        A_proj_pred <- inla.spde.make.A(mesh, loc = mesh$loc[,c(1,2)])
        
        # simulate GRF
        Q = inla.spde2.precision(spde_obj, theta = c(log(sim_range), log(sim_sd) ))
        field <- inla.qsample(Q=Q)
        field <- field - mean(field)
        
        plot_pixels <- pixels(mesh, mask=as(owin(),'SpatialPolygons'), nx=150, ny=150)
        
        # predict the values onto a high res grid
        A_proj_grid <- inla.spde.make.A(mesh, loc=plot_pixels@coords)
        
        plot_pixels$field <- as.numeric(A_proj_grid %*% field)
        
        colsc <- function(...) {
                scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"RdYlBu")),
                                     limits = range(..., na.rm=TRUE))
        }
        
        ggplot() + gg(plot_pixels) + colsc(plot_pixels$field) +
                ggtitle(paste0('SD = ',sim_sd,'  Range = ',sim_range))  
}






# Function for simulating and plotting a GRF on a unit square
GRF_sim <- function(plot.mesh = F, plot.ips = NULL)
{
        # define inla mesh
        spde <- inla.spde2.pcmatern(mesh, constr = T,
                                    prior.range = c(sim_range, 0.01),
                                    prior.sigma = c(sim_sd, 0.01))
        
        # define indices
        index <- inla.spde.make.index('smooth',n.spde = spde$n.spde)
        
        # define prediction matrix
        A_proj_pred <- inla.spde.make.A(mesh, loc = mesh$loc[,c(1,2)])
        
        # simulate GRF
        Q = inla.spde2.precision(spde, theta = c(log(sim_range), log(sim_sd) ))
        field <- inla.qsample(Q=Q)
        field <- field - mean(field)
        
        plot_pixels <- pixels(mesh,
                              mask=starea,
                              nx=300, ny=300)
        
        # predict the values onto a high res grid
        A_proj_grid <- inla.spde.make.A(mesh, loc=plot_pixels@coords)
        
        plot_pixels$field <- as.numeric(A_proj_grid %*% field)
        
        colsc <- function(...) {
                scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"RdYlBu")),
                                     limits = range(..., na.rm=TRUE))
        }
        
        if (plot.mesh) {
                if (is.null(plot.ips)) {
                        ggplot() + gg(plot_pixels) + colsc(plot_pixels$field) + gg(mesh)+
                                ggtitle(paste0('SD = ',sim_sd,'  Range = ',sim_range))
                }else{
                        ggplot() + gg(plot_pixels) + colsc(plot_pixels$field) + gg(mesh)+ gg(plot.ips, alpha = .1)+
                                ggtitle(paste0('SD = ',sim_sd,'  Range = ',sim_range))
                }
        }else{
                  
                if (is.null(plot.ips)) {
                        ggplot() + gg(plot_pixels) + colsc(plot_pixels$field) + 
                                ggtitle(paste0('SD = ',sim_sd,'  Range = ',sim_range))
                }else{
                        ggplot() + gg(plot_pixels) + colsc(plot_pixels$field) + gg(plot.ips, alpha = .1)+
                                ggtitle(paste0('SD = ',sim_sd,'  Range = ',sim_range))
                }
        }
}


sim_range <- size/2
sim_sd <- 0.05

GRF_sim(plot.mesh=F)

