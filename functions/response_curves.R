### Get response curves from inlabru predictions
get.resp.curves <- function(fit, vars, mode = NULL, samp = 100,
                            int = FALSE){
  
  if (class(vars) == "formula") {
    vars <- all.vars(vars)
  }
  
  bvals <- rbind(
    over(ips, env.e),
    over(pts, env.e)
  )
  bvals <- bvals[,vars]
  min.vals <- apply(bvals, 2, min)
  max.vals <- apply(bvals, 2, max)
  mea.vals <- apply(bvals, 2, mean)
  
  pred.list <- list()
  
  resp.data <- matrix(nrow = 100, ncol = length(vars))
  
  for (z in 1:length(vars)) {
    
    resp.data[,] <- rep(mea.vals, each = 100)
    
    cat("Predicting response curve for", vars[z], "\n")
    
    resp.data[,z] <- seq(min.vals[z], max.vals[z], length.out = 100)
    
    if (!is.null(mode)) {
      if (mode == "cloglog") {
        stop("Mode cloglog not available for LGCP only models.")
        pred.resp <- predict(fit,
                             NULL,
                             ~ eval(parse(text =
                                            paste0("1-exp(-exp(",
                                                   paste0(vars,
                                                          "_eval(resp.data[,",
                                                          1:length(vars),
                                                          "])", collapse = "+"),
                                                   ifelse(int, "+intercept_pa_latent", ""),
                                                   "))")
                             )), n.samples = samp
        )
      } else{
        pred.resp <- predict(fit,
                             NULL,
                             ~ eval(parse(text =
                                            paste0(mode, "(",
                                                   paste0(vars,
                                                          "_eval(resp.data[,",
                                                          1:length(vars),
                                                          "])", collapse = "+"),
                                                   ifelse(int, "+Intercept_latent", ""),
                                                   ")")
                             )), n.samples = samp
        )
      }
    } else{
      pred.resp <- predict(fit,
                           NULL,
                           ~ eval(parse(text =
                                          paste0(mode, "(",
                                                 paste0(vars,
                                                        "_eval(resp.data[,",
                                                        1:length(vars),
                                                        "])", collapse = "+"),
                                                 ifelse(int, "+Intercept_latent", ""),
                                                 ")")
                           )), n.samples = samp
      )
    }
    
    pred.resp$base <- resp.data[,z]
    
    pred.list[[z]] <- pred.resp
    
  }
  
  if (!is.null(mode)) {
    tp <- mode
  } else{
    tp <- "lin"
  }
  
  names(pred.list) <- vars
  
  return(structure(pred.list, class = c("respCur", tp, class(pred.list))))
}

# Old version, multiple plots
# plot.respCur <- function(x) {
#   
#   plot.list <- list()
#   
#   sc <- ifelse("lin" %in% class(x), "Linear predictor",
#                ifelse("exp" %in% class(x), "Relative Occurrence Rate", "Probability"))
#   
#   for (z in 1:length(x)) {
#     
#     df <- x[[z]]
#     
#     plot.list[[z]] <- ggplot(df) +
#       geom_line(aes(x = base, y = mean)) +
#       geom_ribbon(aes(x = base, ymin = q0.025, ymax = q0.975), alpha = .4, fill = "#52638E") +
#       scale_x_continuous(expand = c(0,0))+
#       theme_classic() + ylab(sc) + xlab(names(x)[z])
#     
#   }
#   
#   print(eval(parse(text = paste("plot.list[[", 1:length(plot.list), "]]", collapse = "+"))))
#   
#   return(invisible(NULL))
#   
# }

plot.respCur <- function(x, free = NULL, mode = "both") {
  
  sc <- ifelse("lin" %in% class(x), "Linear predictor",
               ifelse("exp" %in% class(x), "Relative Occurrence Rate", "Probability"))
  
  if (is.null(free)) {
    if (sc == "Probability") {
      free = F
    } else {
      free = T
    }
  }
  
  x <-  lapply(names(x), function(z){cbind(x[[z]], var = z)})
  
  dat <- do.call("rbind", x)
  
  if (mode != "both") {
    if (mode == "median") {
      p <- ggplot(dat) +
        geom_line(aes(x = base, y = q0.5))
    } else {
      p <- ggplot(dat) +
        geom_line(aes(x = base, y = mean))
    }
  } else {
    p <- ggplot(dat) +
      geom_line(aes(x = base, y = mean), color = "#233461", linetype = "dashed") +
      geom_line(aes(x = base, y = q0.5))
  }
  
  p <- p + 
    geom_ribbon(aes(x = base, ymin = q0.025, ymax = q0.975), alpha = .4, fill = "#52638E") +
    scale_x_continuous(expand = c(0,0))+
    theme_classic() + ylab(sc) + xlab(NULL) +
    theme(strip.background = element_blank()) +
    facet_wrap(~var, scales = ifelse(free, "free", "free_x"))
  
  
  print(p)
  
  return(invisible(NULL))
  
}

