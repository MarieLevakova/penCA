#' Fit VEC model with penalty on beta
#'
#' This function loads data matrix X and fits a cointegrated VEC model
#' with the cointegration rank r, where elements of the cointegration matrix
#' beta are penalized with lasso penalty.
#'
#' @param X Data matrix
#' @param r cointegration rank
#' @param max.iter maximum number of iterations of the algorithm
#' @param conv convergence parameter
#' @param nlambda number of penalty values lambda to be used
#' @param n.cv number of iterations of the crossvalidation procedures for picking optimal penalty value lambda
#' @param glmnetthresh tolerance parameter for the interior call of glmnet function
#' @param dt timestep
#' @param equal.penalty logical value indicating, whether the same penalty value lambda should be used for all columns of beta
#' @return A list containing the following components:
#' @slot BETA estimated matrix beta
#' @slot ALPHA estimated matrix alpha, recalculated per time unit (alpha/\code{dt})
#' @slot OMEGA the covariance matrix of the random errors, recalculated per time unit (Omega/\code{dt})
#' @slot PI \code{ALPHA*BETA}
#' @slot MU estimated intercepts, recalculated per time unit (mu/\code{dt})
#' @slot it number of iterations
#' @slot value.obj values of the objective functions throughout iterations
#' @slot Pi.it estimate of Pi throughout iterations (vectorized form)
#' @export
#'
pen.beta <- function(X, r, max.iter = 10, conv = 10^-2, nlambda = 100, n.cv = 20,
                     glmnetthresh = 1e-04, dt = 1, equal.penalty = F){

  ## OUTPUT
  # BETA: estimate of cointegrating vectors
  # ALPHA: estimate of adjustment coefficients
  # OMEGA: estimate of covariance matrix
  # PI: estimate of Pi
  # MU: estimate of the intercept
  # it: number of iterations
  # value.obj: value of the objective function at each iteration
  # Pi.it: estimate of Pi at each iteration

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

  # Starting value
  # fit0 <- johansenAggregated(Y = Ystd, Z = Zstd, r = r, dt = dt, intercept = F, normalize = F)
  # beta.init <- fit0$beta
  mm <- diag(as.vector(1/sdZ^2)) %*% cov(Zstd, Ystd) %*% diag(as.vector(1/sdY^2)) %*% cov(Ystd, Zstd)
  beta.init <- matrix(eigen(mm)$vectors[,1:r], ncol = r)
  alpha.init <- t(coef(lm(Ystd ~ I(Zstd %*% beta.init)-1)))
  Pi.init <- matrix(alpha.init, ncol = r) %*% t(matrix(beta.init, ncol = r))
  mu.init <- meanY - Pi.init %*% as.matrix(meanZ, ncol = 1)
  RESID <- Y - matrix(1, nrow = N, ncol = 1) %*% t(mu.init) - Z %*% t(Pi.init)
  Omega.init <- (t(RESID) %*% RESID)/N

  # Convergence parameters: initialization
  it <- 1
  diff.obj <- 10*conv
  value.obj <- matrix(NA, ncol = 1, nrow = max.iter+1)
  Pi.it <- matrix(NA, ncol = p*p, nrow = max.iter + 1)

  RESID <- Y - Z %*% t(Pi.init) - matrix(1, nrow = N, ncol = 1) %*% matrix(mu.init, nrow = 1)
  value.obj[1,] <- (1/N)*sum(diag(RESID %*% solve(Omega.init) %*% t(RESID))) + log(det(Omega.init))
  Pi.it[1,] <- matrix(Pi.init, nrow = 1)

  while((it < max.iter) &  (diff.obj > conv)){
    # Obtain Alpha
    FIT2 <- NTS.ALPHA.Procrusted(Y = Ystd, Z = Zstd,
                                 r = r, Omega = Omega.init, beta = beta.init)
    # Obtain Beta and Omega
    FIT3 <- NTS.BETA(Y = Ystd, Z = Zstd, r = r, Omega = Omega.init,
                     P = FIT2$P, alpha = FIT2$ALPHA, alphastar = FIT2$ALPHAstar,
                     nlambda = nlambda, n.cv = n.cv, glmnetthresh = glmnetthresh,
                     equal.penalty = equal.penalty)

    # Check convergence
    beta.new <- FIT3$BETA
    beta.init <- matrix(beta.new, nrow = p, ncol = r)
    alpha.init <- FIT2$ALPHA
    Pi.init <- alpha.init%*%t(beta.init)
    mu.init <- meanY - Pi.init %*% as.matrix(meanZ, ncol = 1)
    RESID <- Y - matrix(1, nrow = N, ncol = 1) %*% t(mu.init) - Z %*% t(Pi.init)
    Omega.init <- (t(RESID) %*% RESID)/N

    value.obj[1+it,] <- (1/N)*sum(diag((RESID) %*% solve(Omega.init) %*% t(RESID))) + log(det(Omega.init))
    diff.obj <- abs(value.obj[1+it,] - value.obj[it,])/abs(value.obj[it,])
    Pi.it[1+it,] <- matrix(Pi.init, nrow = 1)
    it <- it + 1
  }

  out <- list(BETA = beta.init, ALPHA = alpha.init/dt, OMEGA = Omega.init/dt,
              PI = Pi.init/dt, MU = mu.init/dt, it = it, value.obj = value.obj,
              Pi.it = Pi.it)

  return(out)
}

# Auxiliary functions for solving the problem of penalized beta
NTS.ALPHA.Procrusted <- function(Y, Z, r, Omega, beta){
  ### FUNCTION TO ESTIMATE ALPHA ###

  ## INPUT
  # Y: Response Time Series
  # Z: Predictor time Series in Levels
  # r: cointegration rank
  # Omega: estimate of error covariance matrix
  # beta: estimate of cointegrating vector

  ## OUTPUT
  # ALPHA: estimate of adjustment coefficients
  # ALPHAstar: estimate of transformed adjustment coefficients
  # P: transformation matrix

  ## START CODE
  # Data matrices
  Xmatrix <- Z%*%beta
  decomp <- eigen(Omega)
  P <- (decomp$vectors) %*% diag(1/sqrt(decomp$values)) %*% solve(decomp$vectors)
  Pmin <- (decomp$vectors) %*% diag(sqrt(decomp$values)) %*% solve(decomp$vectors)

  # Singular Value Decomposition to compute ALPHA
  SingDecomp <- svd(t(Xmatrix) %*% Y %*% P)
  ALPHAstar <- t(SingDecomp$u %*% t(SingDecomp$v))
  if (r==1){
    ALPHAstar <- matrix(ALPHAstar, ncol = 1)
  }
  ALPHA <- Pmin %*% ALPHAstar

  out <- list(ALPHA = ALPHA, ALPHAstar = ALPHAstar, P = P)
  return(out)
}

NTS.BETA <- function(Y, Z, r, Omega, P, alpha, alphastar,
                     nlambda, n.cv, glmnetthresh = 1e-04, equal.penalty){
  ### FUNCTIONS TO ESTIMATE BETAs ###
  # First step: Determine Beta conditional on Gamma, alpha and Omega using Lasso Regression
  # Second step: Determine Omega conditional on Gamma, alpha and beta using GLasso

  ## INPUT
  # Y: Response Time Series
  # Z: Time Series in Levels
  # r: cointegration rank
  # Omega: estimate of error covariance matrix
  # P: transformation matrix P derived from Omega
  # alpha: estimate of adjustment coefficients
  # alphastar: estimate of transformed adjustment coefficients
  # nlambda: number of tuning parameter values to try
  # n.cv: number of crossvalidation iterations
  # glmnetthresh: tolerance parameter glmnet function
  # equal.penalty: logical, whether the same penalty value lambda should be used for all columns of beta

  ## OUTPUT
  # BETA: estimate of cointegrating vectors

  # Data matrices
  Ymatrix <- Y %*% t(P) %*% alphastar
  Xmatrix <- Z
  n <- nrow(Ymatrix)

  # Store Estimates of cointegrating vector
  BETA.Sparse <- matrix(NA, ncol = r, nrow = ncol(Y))

  # Standardized Response
  Ymatrixsd <- Ymatrix
  for (i.r in 1:r){
    Ymatrixsd[,i.r] <- Ymatrix[,i.r]/sd(Ymatrix[,i.r])
  }

  if(equal.penalty){

    # Determine the lambda sequence
    lambda.max <- 0
    lambda.min <- NA
    for(i.r in 1:r){
      determine_lambdasequence <- glmnet(y = Ymatrixsd[,i.r], x = Xmatrix,
                                         intercept = F,
                                         family = "gaussian")
      lambda.max <- max(c(lambda.max, determine_lambdasequence$lambda))
      lambda.min <- min(c(lambda.min, determine_lambdasequence$lambda), na.rm = T)
    }

    lambda_restricted <- exp(seq(log(lambda.max), log(lambda.min), length = nlambda))

    cv <- rep(0, nlambda)
    for(i in 1:n.cv){
      for (i.r in 1:r){
        # Calculate crossvalidated error
        BETA.scaled <- cv.glmnet(y = Ymatrixsd[,i.r], x = Xmatrix, lambda = lambda_restricted,
                                 standardize = F, intercept = F, family = "gaussian",
                                 thresh = glmnetthresh)
        cv <- cv + 1/(n.cv*r)*BETA.scaled$cvm
      }
    }

    # Lasso with fixed lambda
    lambda.opt <- lambda_restricted[which.min(cv)]

    for(i.r in 1:r){
      LASSOfinal <- glmnet(y = Ymatrixsd[,i.r], x = Xmatrix, standardize = F, intercept = F,
                           lambda = lambda_restricted, family = "gaussian", thresh = glmnetthresh)
      BETAscaled <- matrix(coef(LASSOfinal, s = lambda.opt), nrow = 1)[-1]
      BETA.Sparse[,i.r] <- BETAscaled*sd(Ymatrix[,i.r])
    }
  } else {

    for(i.r in 1:r){
      # Determine lambda
      cv <- rep(0, nlambda)
      for(i in 1:n.cv){
        BETA.scaled <- cv.glmnet(y = Ymatrixsd[,i.r], x = Xmatrix, nlambda = nlambda,
                                 standardize = F, intercept = F, family = "gaussian",
                                 thresh = glmnetthresh)
        cv <- cv + 1/(n.cv*r)*BETA.scaled$cvm
      }
      # LASSO with varying lambdas
      lambda.opt <- BETA.scaled$lambda[which.min(cv)]
      LASSOfinal <- glmnet(y = Ymatrixsd[,i.r], x = Xmatrix, standardize = F, intercept = F,
                           nlambda = nlambda, family = "gaussian", thresh = glmnetthresh)
      BETAscaled <- matrix(coef(LASSOfinal, s = lambda.opt), nrow = 1)[-1]
      BETA.Sparse[,i.r] <- BETAscaled*sd(Ymatrix[,i.r])
    }
  }

  out <- list(BETA = BETA.Sparse)
}


