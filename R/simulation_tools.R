#' Generate matrix Pi
#'
#' This function generates a random matrix Pi that defines a cointegrated VEC model
#' with a specified cointegration rank through Jordan normal form Pi = P^{-1}*J*P.
#'
#' @param p dimension of the matrix
#' @param r cointegration rank
#' @param sparse logical value indicating if the matrix P should have a sparse form
#' @return pxp matrix Pi
#' @export
generate_Pi <- function(p, r, sparse = TRUE){
  j.diag <- c(-runif(r), rep(0, p-r))
  J <- diag(j.diag)

  if(sparse){
    P <- matrix(1/sqrt(rep(1:p, each = p)), nrow = p, ncol = p)
    for(i1 in 2:p){
      for(i2 in 1:(i1-1)){
        P[i1,i2] <- 0
      }
    }
    P <- P[,sample(1:p)]
  } else {
    repeat{
      P <- matrix(rnorm(p^2), nrow = p)
      if (rankMatrix(P)[1]==p) break
    }
  }

  P.inv <- solve(P)
  Pi <- P.inv %*% J %*% P
  return(Pi)
}

#' Calculate the intercept
#'
#' This function calculates the intercept for a given matrix Pi and data matrix.
#'
#' @param yt N times p matrix of observations of the process
#' @param Pi matrix Pi
#' @param dt timestep of the process
#' @return p-dimensional vector of intercepts
#' @export
estimate.mu <- function(yt, Pi, dt = 1){
  Y <- diff(yt)
  N <- dim(Y)[1]
  P <- dim(Y)[2]
  Z <- as.matrix(yt[-N,])
  apply(Y - Z %*% t(dt *Pi), 2, mean)/dt
}

#' Performance metrics of the estimator
#'
#' This function calculates performance metrics of an estimator
#' with respect to the true matrix Pi and forecasting of independent data.
#'
#' @param Pi.est estimate of Pi
#' @param Pi.true true matrix Pi
#' @param y1 differences of the process (training sample)
#' @param x1 lagged values of the process (training sample)
#' @param y2 differences of the process (testing sample)
#' @param x2 lagged values of the process (testing sample)
#' @return p-dimensional vector of intercepts
#' @export A list with the following components
#' @slot MSPE mean square error of one step prediction of the process
#' @slot MSE mean square error of Pi
#' @slot r.norm l2 norm of Pi.est divided by l2 norm of Pi.true
#' @slot angle principal angle of the space spanned by Pi.est and Pi.true
#' @slot rank rank of Pi.est
#' @slot zeros number of zero elements in Pi.est
#' @slot correct.zeros proportion of correctly identified zeros in Pi.est
#' @slot correct.nonzeros proportion of correctly identified nonzeros in Pi.est
#' @slot SSE sum of squared residuals
perf.metrics <- function(Pi.est, Pi.true, y1, x1, y2, x2){
  tse <- sum(y1^2)

  MSPE <-  sum((y2 - x2 %*% t(Pi.est))^2)/(prod(dim(y2))) # mean square error of 1-step prediction
  MSE <- matrix.mse(Pi.est, Pi.true)
  r.norm <- rel.norm(Pi.est, Pi.true)
  angle <- matrix.distance(Pi.est, Pi.true)
  rank <- rankMatrix(Pi.est)[1]
  zeros <- sum(Pi.est == 0)/prod(dim(Pi.est))
  correct.zeros <- correct.zeros(Pi.est, Pi.true)
  correct.nonzeros <- correct.nonzeros(Pi.est, Pi.true)
  SSE <- sum((y1- x1 %*% t(Pi.est))^2)/tse

  list(MSPE = MSPE, MSE = MSE, r.norm = r.norm, angle = angle, rank = rank,
       zeros = zeros, correct.zeros = correct.zeros, correct.nonzeros = correct.nonzeros,
       SSE = SSE)
}

matrix.mse <- function(A, B){
  sqrt(sum((A-B)^2))/sqrt(sum(B^2))
}

rel.norm <- function(A, B){
  sqrt(sum((A)^2))/sqrt(sum(B^2))
}

matrix.distance <- function(A, B){
  ifelse(all(A==0) | all(B==0), pi/2,
         acos(sum(diag(t(A) %*% B))/(sqrt( sum(diag(t(A) %*% A)) * sum(diag(t(B) %*% B)))))
  )
}

principal_angles <- function(a, b){
  # AUXILIARY FUNCTION
  # Calculate minimal angle between subspace a and b
  # INPUT: subspace a and b
  # OUTPUT: first of the principal angles
  angles <- matrix(0, ncol = ncol(a), nrow = 1)
  qa <- qr.Q(qr(a))
  qb <-  qr.Q(qr(b))
  rkA <- qr(a)$rank
  if(rkA==0){
    return(NA)
  } else {
    rkB <- qr(b)$rank
    C <- svd(t(qa)%*%qb)$d

    if (rkA <= rkB){
      B <- qb - qa%*%(t(qa)%*%qb);
    } else {B <- qa - qb%*%(t(qb)%*%qa)}
    S <- svd(B)$d
    S <- sort(S)

    for (i in 1:min(rkA,rkB)){
      if (C[i]^2 < 0.5) {angles[1,i]  <- acos(C[i])}
      else if (S[i]^2 <= 0.5) {angles[1,i] = asin(S[i])}
    }
    angles <- t(angles)

    return(angles[1])
  }
}

correct.zeros <- function(A, B){
  sum(A[B==0]==0)/sum(B==0)
}

correct.nonzeros <- function(A, B){
  sum(A[B!=0]!=0)/sum(B!=0)
}

