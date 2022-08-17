#' Determine the cointegration rank by Johansen procedure
#'
#' This function loads data matrix X and determines the cointegration rank by
#' likelihood ratio tests proposed by Johansen (trace test or maximum eigenvalue test)
#'
#' @param X Matrix of lagged values of the process
#' @param conf.level confidence level of the test
#' @param type which test statistic should be used, the options are
#' \code{"trace"} (trace test) and \code{"max"} (maximum eigenvalue test).
#' @return A list containing the following components:
#' @slot r cointegration rank
#' @slot lambdas eigenvalues
#' @slot lrts data frame with test statics for ranks 0 to p-1 (rows) and both types of the test (columns "trace" and "max")
#' @export
rank.johansen <- function(X, conf.level = 0.05, type = c("trace", "max")){

  alpha.levels <- c(0.1, 0.05, 0.01)
  P <- dim(X)[2]  # no. of columns

  if (P > 10) print("Dimension exceed the allowed maximum of 10.") else {
    N <- dim(X)[1]  # no. of rows

    cvals <- array(c(6.5, 12.91, 18.9, 24.78, 30.84, 36.25,
                     42.06, 48.43, 54.01, 59, 65.07, 8.18, 14.9, 21.07, 27.14,
                     33.32, 39.43, 44.91, 51.07, 57, 62.42, 68.27, 11.65,
                     19.19, 25.75, 32.14, 38.78, 44.59, 51.3, 57.07, 63.37,
                     68.61, 74.36, 6.5, 15.66, 28.71, 45.23, 66.49, 85.18,
                     118.99, 151.38, 186.54, 226.34, 269.53, 8.18, 17.95,
                     31.52, 48.28, 70.6, 90.39, 124.25, 157.11, 192.84, 232.49,
                     277.39, 11.65, 23.52, 37.22, 55.43, 78.87, 104.2, 136.06,
                     168.92, 204.79, 246.27, 292.65), c(11, 3, 2))

    lambdas <- vecm(diff(X), X[-N,], r = 1, diag(P), diag(P), dt = 1,
                    intercept = T, normalize = F)[P+5,]
    lrts1 <- -N*rev(cumsum(log(1-rev(lambdas))))
    lrts2 <- -N*log(1-lambdas)

    r <- c(NA, NA)
    names(r) <- c("trace", "max")
    r["trace"] <- which(c(lrts1, 0) <= c(cvals[0:P, conf.level == alpha.levels, 2], 1))[1] - 1
    r["max"] <- which(c(lrts2, 0) <= c(cvals[0:P, conf.level == alpha.levels, 1], 1))[1] - 1
    return(list(r = r[type],
         lambdas = lambdas,
         lrts = data.frame(trace = lrts1, max = lrts2,
                           row.names = paste("r =", 0:(P-1)))))
  }
}

#' Determine the cointegration rank by the method of Bunea et al.
#'
#' This function loads data matrix X and determines the cointegration rank by
#' the method proposed by Bunea et al. using a rank penalty
#'
#' @param X Matrix of lagged values of the process
#' @param figure logical value indicating if an illustrative figure should be plotted
#' @return A list containing the following components:
#' @slot r cointegration rank
#' @slot lambdas eigenvalues
#' @slot mu threshold of statistically significant eigenvalues
#' @export
rank.bunea <- function(X, figure = FALSE){
  p <- ncol(X)
  N <- nrow(X)
  l <- rankMatrix(X)[1]

  Delta_phi <- diff(X, differences = 1)
  R0 <- lm(Delta_phi ~ 1)$residuals
  # mu_hat <- matrix(1, ncol = 1, nrow = dim(Delta_phi)[1]) %*% apply(Delta_phi, 2, mean)
  Level_phi <- X[-N,] # removing the last observation to match the dimension with Delta_Y
  R1 <- lm(Level_phi ~ 1)$residuals
  M <- t(R1) %*% R1
  Q <- R1 %*% ginv(M) %*% t(R1)
  decomp <- Re(eigen(t(R0)%*% Q %*% (R0))$values)
  S <- sum((R0 - Q %*% R0)^2)/(N*p-p*l)
  mu <- 2*S*(p+l)

  if(figure){
    plot(decomp, type = "h", xlab = "r", ylab = "")
    abline(h = mu, lty = 2)
  }

  list(r = length(which(decomp >= mu)), lambdas = decomp, mu = mu)
}
