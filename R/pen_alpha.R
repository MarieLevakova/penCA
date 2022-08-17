#' Fit VEC model with penalty on alpha
#'
#' This function loads data matrix X and fits a cointegrated VEC model
#' with the cointegration rank r, where elements of the loading matrix
#' alpha are penalized with elastic net penalty.
#'
#' @param X Data matrix
#' @param r cointegration rank
#' @param dt time step
#' @param equal.penalty logical value indicating, whether the same penalty value lambda should be used for all rows of alpha
#' @param n.lambda number of penalty values lambda to be used
#' @param n.cv number of iterations of the crossvalidation procedures for picking optimal penalty value lambda
#' @param alpha
#' @return A list containing the following components:
#' @slot ALPHA estimated matrix alpha, recalculated per time unit (alpha/\code{dt})
#' @slot BETA estimated matrix beta
#' @slot PI \code{ALPHA*BETA}
#' @slot OMEGA the covariance matrix of the random errors, recalculated per time unit (Omega/\code{dt})
#' @slot MU estimated intercepts, recalculated per time unit (mu/\code{dt})
#' @slot lambda.lasso the selected lasso penalty
#' @slot lambda.ridge the selected ridge penalty
#' @slot lambda.seq the sequence of penalty values
#' @slot ALPHA.seq the sequence of alpha estimates for each penalty value in lambda.seq (not divided by \code{dt})
#' @slot PI.seq the sequence of Pi estimates for each penalty value in lambda.seq (not divided by \code{dt})
#' @export
#'
pen.alpha <- function(X, r, dt = 1, equal.penalty = F, n.lambda = 100, n.cv = 10,
                      alpha = 1){

  ## START CODE
  Y <- diff(X)
  N <- dim(Y)[1]
  p <- dim(Y)[2]
  Z <- X[1:N,]

  # Centering variables, to remove the effect of intercept
  Ystd <- Y
  meanY <- apply(Y, 2, mean)
  sdY <- apply(Y, 2, sd)
  Ystd <- (Y - matrix(1, nrow = N, ncol = 1) %*% t(as.matrix(meanY)))

  meanZ <- apply(Z, 2, mean)
  sdZ <- apply(Z, 2, sd)
  Zstd <- (Z - matrix(1, nrow = N, ncol = 1) %*% t(as.matrix(meanZ)))

  # Fit the usual cointegration model
  fit0 <- PI_johansen(Y = Ystd, Z = Zstd, r = r, dt = dt, intercept = F,
                             normalize = F)

  # Extract beta and create new set of predictors
  BETA <- matrix(fit0$beta, nrow = p, ncol = r)
  Z.r <- matrix(Zstd %*% BETA, nrow = N, ncol = r)

  # Estimate alpha with LASSO penalty

  ALPHA.Sparse <- matrix(NA, nrow = p, ncol = r)

  if(r == 1){
    ALPHA.Sparse  <- matrix(coef(lm(Ystd ~ Z.r - 1)), ncol = 1)
    lambda.final <- alpha.opt <- NA
    lambda.seq <- ALPHA.seq <- PI.seq <- NA
  } else {

    if(equal.penalty){
      lambda.opt <- crit.opt <- rep(NA, length(alpha))
      for(a in 1:length(alpha)){
        alpha.i <- alpha[a]
        # Determine the lambda sequence
        lambda.max <- 0
        lambda.min <- NA
        for(i.r in 1:p){
          determine_lambdasequence <- glmnet(y = Ystd[,i.r], x = Z.r, intercept = F,
                                             family = "gaussian",
                                             alpha = alpha.i)
          lambda.max <- max(c(lambda.max, determine_lambdasequence$lambda))
          lambda.min <- min(c(lambda.min, determine_lambdasequence$lambda), na.rm = T)
        }

        lambda.seq <- exp(seq(log(lambda.max), log(lambda.min), length = n.lambda))
        cv <- rep(0, n.lambda)
        for(i.cv in 1:n.cv){
          for(i.r in 1:p){
            determine_lambda <- cv.glmnet(y = Ystd[,i.r], x = Z.r, intercept = F,
                                          family = "gaussian", lambda = lambda.seq,
                                          alpha = alpha.i)
            cv <- cv + 1/(n.cv*r)*determine_lambda$cvm
          }
        }
        lambda.opt[a] <- lambda.seq[which.min(cv)]
        crit.opt[a] <- min(cv)
      }

      # Fit the final model
      alpha.opt <- alpha[which.min(crit.opt)]
      lambda.final <- lambda.opt[which.min(crit.opt)]

      for(i.r in 1:p){
        LASSOfinal <- glmnet(y = Ystd[,i.r], x = Z.r, intercept = F,
                             lambda = lambda.seq, family = "gaussian",
                             alpha = alpha.opt)
        ALPHA.Sparse[i.r,] <- matrix(coef(LASSOfinal, s = lambda.final), nrow = 1)[-1]
      }

      ALPHA.seq <- matrix(NA, nrow = length(lambda.seq), ncol = p*r)
      PI.seq <- matrix(NA, nrow = length(lambda.seq), ncol = p*p)
      for(i in 1:length(lambda.seq)){
        ALPHA.aux <- matrix(NA, nrow = p, ncol = r)
        for(i.r in 1:p){
          LASSOaux <- glmnet(y = Ystd[,i.r], x = Z.r, intercept = F,
                             lambda = lambda.seq, family = "gaussian",
                             alpha = alpha.opt)
          ALPHA.aux[i.r,] <- matrix(coef(LASSOaux, s = lambda.seq[i]), nrow = 1)[-1]
        }
        PI.aux <- ALPHA.aux %*% t(BETA)
        ALPHA.seq[i,] <- matrix(ALPHA.aux, ncol = 1)
        PI.seq[i,] <- matrix(PI.aux, ncol = 1)
      }

    } else {
      alpha.opt <- lambda.final <- rep(NA, p)
      for(i.r in 1:p){
        for(a in 1:length(alpha)){
          alpha.i <- alpha[a]
          lambda.opt <- crit.opt <- rep(NA, length(alpha))

          determine_lambdasequence <- glmnet(y = Ystd[,i.r], x = Z.r, intercept = F,
                                             family = "gaussian",
                                             alpha = alpha.i)
          lambda.max <- max(c(lambda.max, determine_lambdasequence$lambda))
          lambda.min <- min(c(lambda.min, determine_lambdasequence$lambda), na.rm = T)

          cv <- rep(0, n.lambda)
          for(i.cv in 1:n.cv){
            determine_lambda <- cv.glmnet(y = Ystd[,i.r], x = Z.r, intercept = F,
                                          family = "gaussian", lambda = lambda.seq,
                                          alpha = alpha.i)
            cv <- cv + 1/(n.cv*r)*determine_lambda$cvm
          }
          lambda.opt[a] <- determine_lambda$lambda[which.min(cv)]
          crit.opt[a] <- min(cv)
        }
        alpha.opt[i.r] <- alpha[which.min(crit.opt)]
        lambda.final[i.r] <- lambda.opt[which.min(crit.opt)]

        # Fit the final model
        LASSOfinal <- glmnet(y = Ystd[,i.r], x = Z.r, intercept = F,
                             nlambda = n.lambda, family = "gaussian",
                             alpha = alpha.opt[i.r])
        ALPHA.Sparse[i.r,] <- matrix(coef(LASSOfinal, s = lambda.final[i.r]), nrow = 1)[-1]
      }
    }
  }

  PI <- ALPHA.Sparse %*% t(BETA)
  MU <- meanY - PI %*% as.matrix(meanZ, ncol =1)
  res <- Y - matrix(1, nrow = N, ncol = 1) %*% t(MU) - Z %*% t(PI)
  OMEGA <- (t(res) %*% res)/N

  return(list(ALPHA = ALPHA.Sparse/dt, BETA = BETA, PI = PI/dt, OMEGA = OMEGA/dt,
              MU = MU/dt, lambda.lasso = lambda.final, lambda.ridge = alpha.opt,
              lambda.seq = lambda.seq, ALPHA.seq = ALPHA.seq, PI.seq = PI.seq))
}
