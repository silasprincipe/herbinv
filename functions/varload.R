var.load <- function(folder = NULL, layers = NULL, bath = NULL) {
  
  l.names <- read.table(paste("data/env/", layers, sep = ""), header = F)
  
  r.names <- paste("data/env/",folder, l.names$V2, ".tif", sep = "")
  
  env <- raster::stack(as.character(r.names))
  
  names(env)
  
  rm(l.names, r.names)
  
  if (!is.null(bath)) {
    
    bath.l <- raster(paste("data/env/bath_layers/bath_",bath,".tif", sep=""))
    
    env <- stack(env, bath.l)
  }
  
  return(env)

}