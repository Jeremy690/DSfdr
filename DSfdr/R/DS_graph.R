#' FDR control via DS and MDS in graphical model
#'
#' @description The estimation for gaussian graphical model can be recast as a nodewise
#' regression problem. Thus we apply our DS and MDS approach nodewisely and combine all
#' the nodewise selection results via \emph{or} rule. The nominal level for each nodewise
#' selection is set to be \code{q/2} in order to control the FDR for the whole graph properly.
#' @param data n by p matrix
#' @param num_split Repeated number of DS procedure for MDS
#' @param q FDR control level
#' @return A list containing the selected edges using DS or MDS. Each component in the list is
#' a binary p by p matrix A, of which \code{A_{ij} = 1} indicates the edge between node \code{i}
#' and \code{j} is selected.
#' @export
#' @examples
#' require('glmnet')
#' require('MASS')
#' n = 300; p = 50; rho = 8; a = -0.6; c = 1.5
#' q = 0.2
#' num_split = 50
#' precision = matrix(0, nrow = p, ncol = p)
#' edges_set = matrix(0, nrow = p, ncol = p)
#' ###Construct banded graph
#' for(i in 1:p){
#'  for(j in 1:p){
#'    if(i == j){
#'      precision[i, j] <- 1
#'    }
#'    if(i != j & abs(i - j) <= rho){
#'      precision[i, j] <- sign(a)*abs(a)^(abs(i - j)/c)
#'      edges_set[i, j] <- 1
#'    }
#'   }
#' }
#' min_eigen = min(eigen(precision)$values)
#' if(min_eigen < 0){diag(precision) <- diag(precision) + abs(min_eigen) + 0.005}
#' ### generate samples
#' data <- mvrnorm(n, mu = rep(0,p), Sigma = solve(precision))
#' selected = DS_graph(data, q, num_split)
#' DS_selected = selected$DS_selected_edge
#' MDS_selected = selected$MDS_selected_edge
#' DS_result = fdp_power_graph(DS_selected, edges_set)
#' MDS_result = fdp_power_graph(MDS_selected, edges_set)
#' c(DS_result$fdp, MDS_result$fdr)
#' c(DS_result$power, MDS_result$power)


### nodewise data-splitting procedure
DS_graph <- function(data, q, num_split){
  DS_selected_edge  <- matrix(0, nrow = p, ncol = p)
  MDS_selected_edge <- matrix(0, nrow = p, ncol = p)

  for(j in 1:p){
    ### response variable and design matrix
    y <- data[, j]
    X <- data[, -j]

    inclusion_rate <- matrix(0, nrow = num_split, ncol = p - 1)
    num_select <- rep(0, num_split)

    ### multiple data splits
    for(iter in 1:num_split){
      ### randomly split the data
      sample_index1 <- sample(x = c(1:n), size = 0.5 * n, replace = F)
      sample_index2 <- setdiff(c(1:n), sample_index1)

      ### get the penalty lambda for Lasso
      cvfit  <- cv.glmnet(X[sample_index1, ], y[sample_index1], type.measure = "mse", nfolds = 10)
      lambda <- cvfit$lambda.min

      ### run Lasso on the first half of the data
      beta1 <- as.vector(glmnet(X[sample_index1, ], y[sample_index1], family = "gaussian", alpha = 1, lambda = lambda)$beta)
      nonzero_index = which(beta1 != 0)

      ### run OLS on the second half of the data, restricted on the selected features
      if(length(nonzero_index) != 0){
        beta2 <- rep(0, p - 1)
        fit <- lm(y[sample_index2] ~ X[sample_index2, nonzero_index] - 1)
        beta2[nonzero_index] <- as.vector(fit$coeff)
      }

      ### calculate the mirror statistics
      M <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
      selected_index <- analys_graph(M, abs(M), q/2)

      ### the size of the selected neighborhood
      num_select[iter] <- length(selected_index)
      inclusion_rate[iter, selected_index] <- 1/num_select[iter]
    }

    ### single data-splitting result
    DS_selected_edge[j, -j] <- ifelse(inclusion_rate[1, ] > 0, 1, 0)

    ### multiple data-splitting result
    inclusion_rate <- apply(inclusion_rate, 2, mean)
    feature_rank <- order(inclusion_rate)
    feature_rank <- setdiff(feature_rank, which(inclusion_rate == 0))
    null_feature <- numeric()
    for(feature_index in 1:length(feature_rank)){
      if(sum(inclusion_rate[feature_rank[1:feature_index]]) > q/2){
        break
      }else{
        null_feature <- c(null_feature, feature_rank[feature_index])
      }
    }
    selected_index <- rep(0, p - 1)
    selected_index[setdiff(feature_rank, null_feature)] <- 1
    MDS_selected_edge[j, -j] <- selected_index
  }

  return(list(DS_selected_edge = DS_selected_edge, MDS_selected_edge = MDS_selected_edge))
}
