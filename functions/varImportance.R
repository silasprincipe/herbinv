var.importance <- function(data, cols = 2:length(data), model, n = 5){
        
        results <- matrix(nrow = n, ncol = length(cols))
        
        colnames(results) <- colnames(data[,cols])
        
        for (i in 1:n) {
                
                for (k in cols) {
                        
                        x <- data
                        
                        x[,k] <- x[sample.int(nrow(data)), k]
                        
                        pred <- predict(model, data[, cols])
                        
                        pred.new <- predict(model, x[, cols])
                        
                        results[i, (k-1)] <- (1 - cor(pred, pred.new))
                }
        }
        
        results
}
