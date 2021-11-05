#' Data generating process
#'
#' @description Generating data following linear regression
#'
#' @param n Sample size
#' @param p Number of feature
#' @param p0 Number of relevent features
#' @param covariance Covaraince matrix for the design matrix
#' @param delta Signal strength, the true signal strength is \deqn{delta*(log p/n)^{1/2}}
#'
#' @return A list containing design matrix, response variable and signal index
#' @export
#'
#' @examples
#' require('MASS')
#' n = 500; p = 1000; p0 = 50; delta = 3
#' Sigma =  diag(1, p)
#' data = generate_data(n, p, p0, Sigma, delta)
#' X = data$X
#' Y = data$y
#' signal_index = data$signal_index
generate_data <- function(n, p, p0, covariance, delta){
  ### randomly generate the true beta
  beta_star <- rep(0, p)
  signal_index = sample(1:p, p0, replace = F)
  beta_star[signal_index] <- rnorm(p0, mean = 0, sd = delta*sqrt(log(p)/n))

  ### generate the design matrix X and response y
  X <- mvrnorm(n, mu = rep(0, p), Sigma = covariance)
  y <- X%*%beta_star + rnorm(n, mean = 0, sd = 1)

  return(list(X = X, y = y, signal_index))
}
