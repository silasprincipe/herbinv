



resp.curves <- function (object, predictors, select.columns = NULL, label = NULL, 
          len = 50, lhsample = 100, lwd = 1, ylab = "Occurrence probabilities", 
          method = "mean", disp = "all", overlay.mean = T, col.curves = "grey", 
          col.novel = "grey", col.mean = "black", lwd.known = 2, lwd.mean = 2,
          expand = FALSE,
          ...) 
{
        if (is.null(select.columns)) 
                select.columns = seq_len(ncol(predictors))
        for (i in select.columns) {
                summaries = data.frame(matrix(0, 6, ncol(predictors)))
                for (iz in 1:ncol(predictors)) {
                        summaries[, iz] = summary(predictors[, iz])
                }
                if (method == "stat3") {
                        summaries.j = as.matrix(summaries[c(1, 4, 6), -i], 
                                                ncol = (ncol(predictors) - 1))
                        comb = min(lhsample, 3^(ncol(predictors) - 1))
                        nc = 3
                }
                else if (method == "stat6") {
                        summaries.j = as.matrix(summaries[, -i], ncol = (ncol(predictors) - 
                                                                                 1))
                        comb = min(lhsample, 6^(ncol(predictors) - 1))
                        nc = 6
                }
                else if (method == "mean") {
                        summaries.j = as.matrix(summaries[4, -i], ncol = (ncol(predictors) - 
                                                                                  1))
                        comb = 1
                        nc = 1
                        overlay.mean = F
                }
                dummy.j = as.matrix(predictors[1:len, -i], ncol = (ncol(predictors) - 
                                                                           1))
                if (comb < lhsample) {
                        mat = vector("list", ncol(dummy.j))
                        for (m in 1:ncol(dummy.j)) mat[[m]] = 1:nc
                        mat = expand.grid(mat)
                }
                else {
                        mat = round(qunif(lhs::randomLHS(lhsample, ncol(dummy.j)), 
                                          1, nrow(summaries.j)), 0)
                }
                if (is.null(label)) {
                        label = names(predictors)
                }
                for (r in 1:nrow(mat)) {
                        for (j in 1:ncol(dummy.j)) {
                                dummy.j[, j] = as.vector(rep(summaries.j[mat[r, 
                                                                             j], j], len))
                        }
                        dummy = data.frame(seq(min(predictors[, i]), max(predictors[, 
                                                                                    i]), length = len), dummy.j)
                        
                        if (isTRUE(expand)) {
                                dummy = data.frame(seq((min(predictors[, i])-(min(predictors[, i])/10)),
                                                       (max(predictors[, i])+(max(predictors[, i])/10)),
                                                       length = len), dummy.j)
                        }
                        
                        names(dummy)[-1] = names(predictors)[-i]
                        names(dummy)[1] = names(predictors)[i]
                        curves <- predict(object, dummy, type = "response")
                        if (disp == "all") {
                                if (r == 1) {
                                        if (i == 1) {
                                                plot(dummy[, names(predictors)[i]], curves, 
                                                     type = "l", ylim = c(0, 1), xlab = label[i], 
                                                     ylab = ylab, lwd = lwd, col = col.curves, 
                                                     ...)
                                                abline(v = min(predictors[, i]), col = 'green')
                                                abline(v = max(predictors[, i]), col = 'green')}
                                                
                                        else {plot(dummy[, names(predictors)[i]], curves, 
                                                  type = "l", ylim = c(0, 1), xlab = label[i], 
                                                  ylab = "", lwd = lwd, col = col.curves, 
                                                  ...)
                                                abline(v = min(predictors[, i]), col = 'green')
                                                abline(v = max(predictors[, i]), col = 'green')
                                                }
                                }
                                else lines(dummy[, names(predictors)[i]], curves, 
                                           lwd = lwd, col = col.curves, ...)
                        }
                        if (disp == "eo.mask") {
                                novel = eo.mask(predictors, dummy)
                                curves.known = curves
                                curves.known[novel == 1] = NA
                                curves.novel = curves
                                curves.novel[novel == 0] = NA
                                if (r == 1) {
                                        if (i == 1) {
                                                plot(dummy[, names(predictors)[i]], curves.known, 
                                                     type = "l", ylim = c(0, 1), xlab = label[i], 
                                                     ylab = ylab, lwd = lwd.known, col = col.curves, 
                                                     ...)
                                                lines(dummy[, names(predictors)[i]], curves.novel, 
                                                      lwd = lwd, col = col.novel, lty = "dotted", 
                                                      ...)
                                        }
                                        else {
                                                plot(dummy[, names(predictors)[i]], curves.known, 
                                                     type = "l", ylim = c(0, 1), xlab = label[i], 
                                                     ylab = "", lwd = lwd.known, col = col.curves, 
                                                     ...)
                                                lines(dummy[, names(predictors)[i]], curves.novel, 
                                                      lwd = lwd, col = col.novel, lty = "dotted", 
                                                      ...)
                                        }
                                }
                                else {
                                        lines(dummy[, names(predictors)[i]], curves.known, 
                                              lwd = lwd.known, col = col.curves, ...)
                                        lines(dummy[, names(predictors)[i]], curves.novel, 
                                              lwd = lwd, col = col.novel, lty = "dotted", 
                                              ...)
                                }
                        }
                }
                if (overlay.mean == T) {
                        dummy = predictors[1:len, ]
                        dummy[, i] = seq(min(predictors[, i]), max(predictors[, 
                                                                              i]), length = len)
                        
                        if (isTRUE(expand)) {
                                dummy[, i] = seq((min(predictors[, i])-(min(predictors[, i])/10)),
                                                       (max(predictors[, i])+(max(predictors[, i])/10)),
                                                       length = len)
                        }
                        
                        for (j in 1:ncol(predictors)) {
                                if (j != i) {
                                        dummy[, j] = rep(mean(predictors[, j]), len)
                                }
                        }
                        curves <- predict(object, dummy, type = "response")
                        lines(dummy[, names(predictors)[i]], curves, lwd = lwd.mean, 
                              col = col.mean, ...)
                }
        }
}
