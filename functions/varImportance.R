var.importance <- function(data, cols = 2:length(data), model, n = 5,
                           get.mean = T){
        
        results <- matrix(nrow = n, ncol = length(cols))
        
        colnames(results) <- colnames(data[,cols])
        
        for (i in 1:n) {
                
                for (k in cols) {
                        
                        x <- data
                        
                        x[,k] <- x[sample.int(nrow(data)), k]
                        
                        pred <- predict(model, data[, cols], type = "response")
                        
                        pred.new <- predict(model, x[, cols], type = "response")
                        
                        if (min(cols) == 1) {
                                results[i, k] <- (1 - cor(pred, pred.new))
                        } else{
                                results[i, (k-1)] <- (1 - cor(pred, pred.new))
                        }
                        
                        
                }
        }
        
        if (isTRUE(get.mean)) {
                round(apply(results, 2, mean),3)
        }else{
                results
        }
}
