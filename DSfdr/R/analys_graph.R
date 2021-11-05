#' Select features based on mirror stats.
#'
#'
#' @param mm mirror statistics
#' @param ww absolute value of mirror statistics
#' @param q FDR control level
#' @return selected index

analys_graph <- function(mm, ww, q){
  ### mm: mirror statistics
  ### ww: absolute value of mirror statistics
  ### q:  FDR control level
  cutoff_set <- max(ww)
  for(t in ww){
    ps <- length(mm[mm > t])
    ng <- length(na.omit(mm[mm < -t]))
    rto <- (ng)/max(ps, 1)
    if(rto <= q){
      cutoff_set <- c(cutoff_set, t)
    }
  }
  cutoff <- min(cutoff_set)
  selected_index <- which(mm > cutoff)

  return(selected_index)
}
