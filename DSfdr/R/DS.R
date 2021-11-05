#' FDR control via DS and MDS in linear regression
#'
#' @description DS procedure first splits the data into two halves. In the first
#' half of the data, we apply LASSO for which the penalization parameter
#' \code{lambda} is chosen via 10-fold cross validation. In the second half
#' of the data, we apply OLS, and get the least squares estimator. We combine
#' the LASSO estimator and OLS estimatr via mirror statistics proposed in
#' \url{https://arxiv.org/abs/2002.08542}.
#' @param X Design matrix
#' @param y Response vector
#' @param num_split Repeated number of DS procedure for MDS
#' @param q FDR control level
#' @return A list containing the selected variables using DS or MDS.
#' @export
#' @examples
#' require('glmnet')
#' require('MASS')
#' n = 500
#' p = 1000
#' p0 = 50
#' delta = 3
#' Sigma = diag(1, p)
#' data = generate_data(n, p, p0, Sigma, delta)
#' X = data$X
#' Y = data$y
#' signal_index = data$signal_index
#' num_split = 50
#' q = 0.1
#' selection = DS(X, Y, num_split, q)
#' DS_selected = selection$DS_feature
#' MDS_selected = selection$MDS_feature
#' DS_result = fdp_power(DS_selected, signal_index)
#' MDS_result = fdp_power(MDS_selected, signal_index)
#' c(DS_result$fdp, MDS_result$fdp)
#' c(DS_result$power, MDS_result$power)

DS <- function(X, y, num_split, q){
  n <- dim(X)[1]; p <- dim(X)[2]
  inclusion_rate <- matrix(0, nrow = num_split, ncol = p)
  fdp <- rep(0, num_split)
  power <- rep(0, num_split)
  num_select <- rep(0, num_split)

  for(iter in 1:num_split){
    ### randomly split the data
    sample_index1 <- sample(x = c(1:n), size = 0.5 * n, replace = F)
    sample_index2 <- setdiff(c(1:n), sample_index1)

    ### get the penalty lambda for Lasso
    cvfit <- cv.glmnet(X[sample_index1, ], y[sample_index1], type.measure = "mse", nfolds = 10)
    lambda <- cvfit$lambda.min

    ### run Lasso on the first half of the data
    beta1 <- as.vector(glmnet(X[sample_index1, ], y[sample_index1], family = "gaussian", alpha = 1, lambda = lambda)$beta)
    nonzero_index <- which(beta1 != 0)
    if(length(nonzero_index)!=0){
      ### run OLS on the second half of the data, restricted on the selected features
      beta2 <- rep(0, p)
      beta2[nonzero_index] <- as.vector(lm(y[sample_index2] ~ X[sample_index2, nonzero_index] - 1)$coeff)

      ### calculate the mirror statistics
      M <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
      # M <- abs(beta1 + beta2) - abs(beta1 - beta2)

      selected_index = analys(M, abs(M), q)
      DS_selected_index = selected_index

      ### number of selected variables
      if(length(selected_index)!=0){
        num_select[iter] <- length(selected_index)
        inclusion_rate[iter, selected_index] <- 1/num_select[iter]
      }
    }else{
      DS_selected_index = NULL
    }
  }


  ### multiple data-splitting (MDS) result
  inclusion_rate <- apply(inclusion_rate, 2, mean)

  ### rank the features by the empirical inclusion rate
  feature_rank <- order(inclusion_rate)
  feature_rank <- setdiff(feature_rank, which(inclusion_rate == 0))
  if(length(feature_rank)!=0){
    null_feature <- numeric()

    ### backtracking
    for(feature_index in 1:length(feature_rank)){
      if(sum(inclusion_rate[feature_rank[1:feature_index]]) > q){
        break
      }else{
        null_feature <- c(null_feature, feature_rank[feature_index])
      }
    }
    MDS_selected_index <- setdiff(feature_rank, null_feature)
  }else{
    MDS_selected_index = NULL
  }
  return(list(DS_feature = DS_selected_index, MDS_feature = MDS_selected_index))
}

