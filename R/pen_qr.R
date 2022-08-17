#' Fit VEC model with penalty imposed on QR decomposition of Pi
#'
#' This function loads data matrix X and fits a VEC model,
#' where matrix Pi is factorized by QR decomposition and a row-lasso penalty
#' is imposed on the matrix R
#'
#' @param X Data matrix
#' @param crit Criterion for choosing optimal penalty value, the options are
#' \code{"fixed"} (given lambda), \code{"CV"} (crossvalidated error), \code{"AIC"} (Akaike information criterion),
#' \code{"BIC"} (Bayes information criterion) and \code{"HQ"} (Hannan-Quinn criterion).
#' @param lambda given penalty value when \code{crit="fixed"}
#' @param n.lambda number of penalty values lambda to be used
#' @param n.cv number of repetitions of the crossvalidation procedure
#' @param thresh convergence parameter
#' @param maxit maximum number of iterations
#' @param dt timestep of the process
#' @param psi power of the weights in the objective function
#' @return A list containing the following components:
#' @slot PI estimated matrix Pi, recalculated per time unit (Pi/\code{dt})
#' @slot MU estimated intercepts, recalculated per time unit (mu/\code{dt})
#' @slot OMEGA the covariance matrix of the random errors, recalculated per time unit (Omega/\code{dt})
#' @slot lambda the selected lasso penalty value
#' @slot S the matrix Q
#' @slot R the matrix R
#' @slot lambda.seq the sequence of penalty values
#' @slot PI.seq the sequence of Pi estimates for each penalty value in lambda.seq (not divided by \code{dt})
#' @slot R.seq the sequence of R estimates for each penalty value in lambda.seq
#' @export
pen.qr <- function(X, crit = c("fixed", "CV", "AIC", "BIC", "HQ"),
                   lambda, n.lambda = 100, n.cv = 10, thresh = 1e-12,
                   maxit = 1e6, dt = 1, psi = 2){
  Y <- diff(X)
  N <- dim(Y)[1]
  p <- dim(Y)[2]
  Z <- X[1:N,]

  # Standardize variables
  meanY <- apply(Y, 2, mean)
  sdY <- apply(Y, 2, sd)
  Ystd <- (Y - matrix(1, nrow = N, ncol = 1) %*% t(as.matrix(meanY))) #%*% diag(1/sdY)

  meanZ <- apply(Z, 2, mean)
  sdZ <- apply(Z, 2, sd)
  Zstd <- (Z - matrix(1, nrow = N, ncol = 1) %*% t(as.matrix(meanZ))) #%*% diag(1/sdZ)

  # Preliminary OLS estimate
  Pi0 <- PI.OLS(Ystd, Zstd, dt = dt, intercept = F)$PI
  QR.Pi0 <- qr(t(Pi0))
  Smat <- qr.Q(QR.Pi0)
  Rmat <- qr.R(QR.Pi0)

  # Construct a new predictor matrix
  ZS <- Zstd %*%Smat

  R.Sparse <- matrix(NA, nrow = p, ncol = p)

  # Determine the lambda sequence
  penalty.factor <- rep(NA, p)
  for(i in 1:p){
    penalty.factor[i] <- sqrt(sum((Rmat[i,i:p])^2))
  }
  determine_lambdasequence <- glmnet(y = Ystd, x = ZS,
                                     intercept = F,
                                     penalty.factor = 1/(penalty.factor^psi),
                                     family = "mgaussian",
                                     nlambda = n.lambda)

  lambda.seq <- determine_lambdasequence$lambda
  n.lambda <- length(lambda.seq)

  if(crit=="fixed"){
    lambda.opt[i.r] <- lambda
  } else {


    if(crit!="CV"){
      R.select <- array(0, dim = c(p, p, n.lambda))
      for(i in 1:p){
        R.select[,i,] <- as.matrix(determine_lambdasequence$beta[[i]])
      }

      aic <- rep(NA, n.lambda)
      bic <- rep(NA, n.lambda)
      hq <- rep(NA, n.lambda)
      for(ii in 1:n.lambda){
        N <- dim(Ystd)[1]
        k <- sum((Smat %*% R.select[,,ii])!=0)
        res <- Ystd - ZS %*% R.select[,,ii]
        Omega.select <- (t(res) %*% res)/N
        Omega.inv <- solve(Omega.select)
        aic[ii] <- N*p*log(2*pi) + N*log(det(Omega.select)) + 2*k +
          sum(diag(res %*% Omega.inv %*% t(res)))
        bic[ii] <- N*p*log(2*pi) + N*log(det(Omega.select)) + k*log(N*p) +
          sum(diag(res %*% Omega.inv %*% t(res)))
        hq[ii] <- N*p*log(2*pi) + N*log(det(Omega.select)) + 2*k*log(log(N*p)) +
          sum(diag(res %*% Omega.inv %*% t(res)))
      }
      if(crit=="AIC"){
        lambda.opt <- lambda.seq[which.min(aic)]
      }
      if(crit=="BIC"){
        lambda.opt <- lambda.seq[which.min(bic)]
      }
      if(crit=="HQ"){
        lambda.opt <- lambda.seq[which.min(hq)]
      }
    }
    if(crit=="CV"){
      cv <- rep(0, n.lambda)
      for(i.cv in 1:n.cv){
        determine_lambda <- cv.glmnet(y = Ystd, x = ZS,
                                      intercept = F,
                                      penalty.factor = 1/(penalty.factor^psi),
                                      lambda = lambda.seq,
                                      family = "mgaussian")
        cv <- cv + 1/(n.cv*p)*determine_lambda$cvm
      }
      lambda.opt <- lambda.seq[which.min(cv)]
    }
  }

  # Fit the final model
  LASSOfinal <- glmnet(y = Ystd, x = ZS,
                       intercept = F,
                       penalty.factor = 1/(penalty.factor^psi),
                       family = "mgaussian")
  R.Sparse <- matrix(0, nrow = p, ncol = p)
  for(i in 1:p){
    R.Sparse[,i] <- as.matrix(coef(LASSOfinal, s = lambda.opt)[[i]])[-1]
  }
  PI.Sparse <- t(Smat %*% R.Sparse)
  R.seq <- matrix(NA, nrow = length(lambda.seq), ncol = p^2)
  PI.seq <- matrix(NA, nrow = length(lambda.seq), ncol = p^2)
  for(i in 1:length(lambda.seq)){
    R.aux <- matrix(NA, nrow = p, ncol = p)
    for(i.p in 1:p){
      R.aux[,i.p] <- as.matrix(coef(LASSOfinal, s = lambda.seq[i])[[i.p]])[-1]
    }
    Pi.aux <- t(Smat %*% R.aux)
    R.seq[i,] <- matrix(R.aux, ncol = 1)
    PI.seq[i,] <- matrix(Pi.aux, ncol = 1)
  }
  # for(i in 1:p){
  #   PI.Sparse[,i] <- PI.Sparse[,i] %*% diag(sdY[i]/sdZ)
  # }

  mu.hat <- meanY - PI.Sparse %*% as.matrix(meanZ, ncol =1)
  res <- Y - matrix(1, nrow = N, ncol = 1) %*% t(mu.hat) - Z %*% t(PI.Sparse)
  Omega.hat <- (t(res) %*% res)/N

  return(list(PI = PI.Sparse/dt, MU = mu.hat/dt, OMEGA = Omega.hat/dt,
              lambda = lambda.opt, S = Smat, R = R.Sparse,
              lambda.seq = lambda.seq, PI.seq = PI.seq, R.seq = R.seq))
}
