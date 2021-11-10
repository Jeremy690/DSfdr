#' FDP and power calculation for gaussian graphical model.
#'
#'
#' @param selected_edge Binary matrix of which each element indicates whether the edge is
#' selected by DS or MDS.
#' @param true_edge Binary matrix of which each element indicates whether the edge exists.
#' @export
#' @return A list containing FDP and power.

fdp_power_graph <- function(selected_edge, true_edge){
  num_false_discoveries <- 0
  num_selected_edge <- 0
  p = dim(true_edge)[1]
  for(i in 1:(p - 1)){
    for(j in (i + 1):p){
      if(selected_edge[i, j] == 1 | selected_edge[j, i] == 1){
        num_selected_edge <- num_selected_edge + 1
        if(true_edge[i, j] == 0){
          num_false_discoveries <- num_false_discoveries + 1
        }
      }
    }
  }
  fdp <- num_false_discoveries/num_selected_edge
  power <- (num_selected_edge - num_false_discoveries)/sum(true_edge)*2

  return(list(fdp = fdp, power = power))
}
