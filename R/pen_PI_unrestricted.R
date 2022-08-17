#' Fit VEC model with lasso penalty and restricted rank
#'
#' This function loads data matrix X and fits a VEC model,
#' where elements of the matrix Pi are penalized with lasso penalty
#' while at the same time the rank of Pi is restricted to a given value.
#'
#' @param yt Data matrix
#' @param crit Criterion for choosing optimal penalty value, the options are
#' \code{"fixed"} (given lambda), \code{"CV"} (crossvalidated error), \code{"AIC"} (Akaike information criterion),
#' \code{"BIC"} (Bayes information criterion) and \code{"HQ"} (Hannan-Quinn criterion).
#' @param lambda given lasso penalty when \code{crit="fixed"}
#' @param alpha vector of ridge penalties
#' @param n.lambda number of penalty values lambda to be used
#' @param n.cv number of repetitions of the crossvalidation procedure
#' @param thresh convergence parameter
#' @param maxit maximum number of iterations
#' @param dt timestep of the process
#' @return A list containing the following components:
#' @slot PI estimated matrix Pi, recalculated per time unit (Pi/\code{dt})
#' @slot MU estimated intercepts, recalculated per time unit (mu/\code{dt})
#' @slot OMEGA the covariance matrix of the random errors, recalculated per time unit (Omega/\code{dt})
#' @slot lambda the selected lasso penalty value
#' @slot alpha the selected ridge penalty value
#' @slot lambda.seq the sequence of penalty values
#' @slot PI.seq the sequence of Pi estimates for each penalty value in lambda.seq (not divided by \code{dt})
#' @export
pen.PI.unrestricted <- function(yt, crit = c("fixed", "CV", "AIC", "BIC", "HQ"),
                                lambda, alpha = 1, n.lambda = 100, n.cv = 10, thresh = 1e-12,
                                maxit = 1e6, dt = 1){
  Y <- diff(yt)
  N <- dim(Y)[1]
  p <- dim(Y)[2]
  Z <- yt[1:N,]

  # Standardize variables
  meanY <- apply(Y, 2, mean)
  sdY <- apply(Y, 2, sd)
  Ystd <- (Y - matrix(1, nrow = N, ncol = 1) %*% t(as.matrix(meanY))) %*% diag(1/sdY)

  meanZ <- apply(Z, 2, mean)
  sdZ <- apply(Z, 2, sd)
  Zstd <- (Z - matrix(1, nrow = N, ncol = 1) %*% t(as.matrix(meanZ))) %*% diag(1/sdZ)

  PI.Sparse <- matrix(NA, nrow = p, ncol = p)

  lambda.opt <- crit.opt <- rep(NA, length(alpha))

  for(a in 1:length(alpha)){
    alpha.i <- alpha[a]

    # Determine the lambda sequence
    lambda.max <- 0
    lambda.min <- NA
    for(i.r in 1:p){
      determine_lambdasequence <- glmnet(y = Ystd[,i.r], x = Zstd,
                                         intercept = F,
                                         family = "gaussian",
                                         alpha = alpha.i)
      lambda.max <- max(c(lambda.max, determine_lambdasequence$lambda))
      lambda.min <- min(c(lambda.min, determine_lambdasequence$lambda), na.rm = T)
    }

    lambda.seq <- exp(seq(log(lambda.max), log(lambda.min), length = n.lambda))

    if(crit=="fixed"){
      lambda.opt[a] <- lambda
    } else {

      pi.select <- array(NA, dim = c(p, p, n.lambda))
      if(crit=="AIC"){
        aic <- rep(NA, n.lambda)
        for(i.r in 1:p){
          penalty.factor <- rep(1, p)
          penalty.factor[i.r] <- 0
          determine_lambda <- glmnet(y = Ystd[,i.r], x = Zstd,
                                     intercept = F,
                                     family = "gaussian",
                                     lambda = lambda.seq,
                                     alpha = alpha.i)
          pi.select[,i.r,] <- matrix(determine_lambda$beta, nrow = p)
        }
        for(ii in 1:n.lambda){
          N <- length(Ystd[,i.r])
          k <- sum(pi.select[,,ii]!=0)
          res <- Ystd - Zstd %*% pi.select[,,ii]
          Omega.select <- (t(res) %*% res)/N
          Omega.inv <- solve(Omega.select)
          aic[ii] <- N*p*log(2*pi) + N*log(det(Omega.select)) + 2*k +
            sum(diag(res %*% Omega.inv %*% t(res)))
        }
        crit.opt[a] <- min(aic)
        lambda.opt[a] <- lambda.seq[which.min(aic)]
      }

      if(crit=="BIC"){
        bic <- rep(NA, n.lambda)
        for(i.r in 1:p){
          penalty.factor <- rep(1, p)
          penalty.factor[i.r] <- 0
          determine_lambda <- glmnet(y = Ystd[,i.r], x = Zstd,
                                     intercept = F,
                                     family = "gaussian",
                                     lambda = lambda.seq,
                                     alpha = alpha.i)
          pi.select[,i.r,] <- matrix(determine_lambda$beta, nrow = p)
        }
        for(ii in 1:n.lambda){
          N <- length(Ystd[,i.r])
          k <- sum(pi.select[,,ii]!=0)
          res <- Ystd - Zstd %*% pi.select[,,ii]
          Omega.select <- (t(res) %*% res)/N
          Omega.inv <- solve(Omega.select)
          bic[ii] <- N*p*log(2*pi) + N*log(det(Omega.select)) + k*log(N*p) +
            sum(diag(res %*% Omega.inv %*% t(res)))
        }
        crit.opt[a] <- min(bic)
        lambda.opt[a] <- lambda.seq[which.min(bic)]
      }

      if(crit=="HQ"){
        hq <- rep(NA, n.lambda)
        for(i.r in 1:p){
          penalty.factor <- rep(1, p)
          penalty.factor[i.r] <- 0
          determine_lambda <- glmnet(y = Ystd[,i.r], x = Zstd,
                                     intercept = F,
                                     family = "gaussian",
                                     lambda = lambda.seq,
                                     alpha = alpha.i)
          pi.select[,i.r,] <- matrix(determine_lambda$beta, nrow = p)
        }
        for(ii in 1:n.lambda){
          N <- length(Ystd[,i.r])
          k <- sum(pi.select[,,ii]!=0)
          res <- Ystd - Zstd %*% pi.select[,,ii]
          Omega.select <- (t(res) %*% res)/N
          Omega.inv <- solve(Omega.select)
          hq[ii] <- N*p*log(2*pi) + N*log(det(Omega.select)) + 2*k*log(log(N*p)) +
            sum(diag(res %*% Omega.inv %*% t(res)))
        }
        crit.opt <- min(hq)
        lambda.opt <- lambda.seq[which.min(hq)]
      }

      if(crit=="CV"){
        cv <- rep(0, n.lambda)
        for(i.cv in 1:n.cv){
          for(i.r in 1:p){
            penalty.factor <- rep(1, p)
            penalty.factor[i.r] <- 0
            determine_lambda <- cv.glmnet(y = Ystd[,i.r], x = Zstd,
                                          intercept = F,
                                          family = "gaussian",
                                          lambda = lambda.seq,
                                          alpha = alpha.i)
            cv <- cv + 1/(n.cv*p)*determine_lambda$cvm
          }
        }
        crit.opt <- min(cv)
        lambda.opt <- lambda.seq[which.min(cv)]
      }
    }
  }

  alpha.opt <- alpha[which.min(crit.opt)]
  lambda.final <- lambda.opt[which.min(crit.opt)]

  # Fit the final model
  PI.seq <- matrix(NA, nrow = length(lambda.seq), ncol = p^2)
  for(i.r in 1:p){
    penalty.factor <- rep(1, p)
    penalty.factor[i.r] <- 0
    LASSOfinal <- glmnet(y = Ystd[,i.r], x = Zstd,
                         standardize = F,
                         intercept = F,
                         alpha = alpha.opt,
                         family = "gaussian")
    for(i in 1:length(lambda.seq)){
      PI.seq[i, (i.r-1)*p+(1:p)] <- matrix(coef(LASSOfinal, s = lambda.seq[i]), nrow = 1)[-1] %*% diag(sdY[i.r]/sdZ)
    }
    PI.Sparse[i.r,] <- matrix(coef(LASSOfinal, s = lambda.final), nrow = 1)[-1] %*% diag(sdY[i.r]/sdZ)
  }

  mu.hat <- meanY - PI.Sparse %*% as.matrix(meanZ, ncol =1)
  res <- Y - matrix(1, nrow = N, ncol = 1) %*% t(mu.hat) - Z %*% t(PI.Sparse)
  Omega.hat <- (t(res) %*% res)/N

  return(list(PI = PI.Sparse/dt, MU = mu.hat/dt, OMEGA = Omega.hat/dt,
              lambda = lambda.final, alpha = alpha.opt, lambda.seq = lambda.seq,
              PI.seq = PI.seq))
}
