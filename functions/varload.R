var.load <- function(folder = "crop_layers/", layers = NULL, bath = NULL,
                     bioclim = NULL, bioclimfut = NULL) {
  
  if (length(layers) == 1) {
    if (grepl(".txt", layers)) {
      l.names <- read.table(paste("data/env/", layers, sep = ""), header = F)[,2]
    }else{
      l.names <- list.files(paste0("data/env/",folder))
      l.names <- l.names[grep(layers, l.names)]
    }
  } else{
    l.names <- list.files(paste0("data/env/",folder))
    l.names <- l.names[grep(paste(layers, collapse = "|"), l.names)]
  }
  
  r.names <- paste("data/env/",folder, l.names,
                   ifelse(grepl(".tif", l.names[1]),
                          "", ".tif"), sep = "")
  
  env <- raster::stack(as.character(r.names))
  
  #names(env)
  
  if (!is.null(bath)) {
    
    bath.l <- raster(paste("data/env/bath_layers/bath_",bath,".tif", sep=""))
    
    env <- raster::stack(env, bath.l)
  }
  
  if (!is.null(bioclim)) {
    
    bioc.f <- list.files("data/env/bioclim_layers/", full.names = T)
    
    if (length(bioclim) > 1) {
      bioclim <- paste(bioclim, collapse = "|")
    }
    
    bioc.f <- bioc.f[grep(bioclim, bioc.f)]
    
    if (!is.null(bioclimfut)) {
      bioc.f <- bioc.f[grep(bioclimfut, bioc.f)]
    } else{
      bioc.f <- bioc.f[grep("current", bioc.f)]
    }
    
    bioc.l <- raster::stack(bioc.f)
    
    if (!compareRaster(bioc.l, env, stopiffalse = F)) {
      bioc.l <- raster::extend(bioc.l, env)
    }
    
    env <- raster::stack(env, bioc.l)
    
    env <- raster::mask(env, env[[1]])
  }
  
  return(env)

}