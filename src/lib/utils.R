#' Compute Pearson's r for a specific cell line experiment
#' and a specific method
#' 
compute_corr <- function(out, iter, cell_line, method){
  out_cor <- vector(mode = "numeric", length = iter)
  for (i in 1:iter){
    out_cor[i] <- stats::cor(out[[i]]$test_pred, out[[i]]$test$y)
  }
  
  corr <- data.frame(as.matrix(out_cor), 
                     as.matrix(rep(cell_line, iter)), 
                     as.matrix(rep(method, iter)),
                     stringsAsFactors = TRUE)
  return(corr)
}