#' Fit VEC model with OLS
#'
#' This function loads data matrix X and fits a VEC model
#' using the ordinary least squares method.
#'
#' @param Y matrix of first differences of the process
#' @param Z matrix of lagged values of the process
#' @param dt timestep of the process
#' @param intercept logical value indicating whether the model includes the intercept or not
#' @return A list containing the following components:
#' @slot PI estimated matrix Pi, recalculated per time unit (Pi/\code{dt})
#' @slot MU estimated intercepts, recalculated per time unit (mu/\code{dt})
#' @slot OMEGA the covariance matrix of the random errors, recalculated per time unit (Omega/\code{dt})
#' @export
PI.OLS <- function(Y, Z, dt = 1, intercept = T){
  N <- dim(Y)[1]
  P <- dim(Y)[2]

  if(intercept){
    m0 <- lm(Y ~ Z)
    est.coef <- m0$coefficients
    est.Omega <- (t(m0$residuals) %*% m0$residuals)/N
    out <- list(PI = t(est.coef[-1,])/dt, MU = est.coef[1,]/dt, OMEGA = est.Omega/dt)
  } else {
    m0 <- lm(Y ~ Z-1)
    est.coef <- m0$coefficients
    est.Omega <- (t(m0$residuals) %*% m0$residuals)/N
    out <- list(PI = t(est.coef)/dt, MU = est.coef[1,]/dt, OMEGA = est.Omega/dt)
  }
}

#' Fit VEC model with Johansen procedure
#'
#' This function loads data matrix X and fits a VEC model
#' with a restricted cointegration rank using the likelihood-based approach
#' of Søren Johansen.
#'
#' @param Y matrix of first differences of the process
#' @param Z matrix of lagged values of the process
#' @param r cointegration rank
#' @param A matrix of linear restrictions on the space spanned by the loading matrix
#' @param H matrix of linear restriction on the space spanned by the cointegration matrix
#' @param dt timestep of the process
#' @param intercept logical value indicating whether the model includes the intercept or not
#' @param normalize logical value indicating whether the cointegration matrix should be normalized to contain a unit matrix as the upper rxr block
#' @return A list containing the following components:
#' @slot N number of observations of the process
#' @slot p dimension of the process
#' @slot r cointegration rank
#' @slot alpha estimate of the loading matrix, recalculated per time unit (alpha/\code{dt})
#' @slot beta estimate of the cointegration matrix
#' @slot Psi estimated intercepts, recalculated per time unit (mu/\code{dt})
#' @slot Omega the covariance matrix of the random errors, recalculated per time unit (Omega/\code{dt})
#' @slot test test statistics of the trace test
#' @slot lambda eigenvalues resulting from calculating the estimate of the cointegration matrix beta
#' @slot A linear restrictions imposed on the space spanned by the loading matrix alpha
#' @slot H linear restrictions imposed on the space spanned by the cointegration matrix beta
#' @slot df degrees of freedom
#' @slot dt timestep of the process
#' @slot r2 proportion of the variability explained by the model
#' @export
PI_johansen <- function(Y = NULL, Z = NULL, r = 1, A = NULL, H = NULL,
                        dt = 0.1, intercept = TRUE, normalize = TRUE){
  N = nrow(Y)
  p = ncol(Y)

  if(is.null(A)){
    A = as.matrix(diag(p)) # produces a unit matrix
  }
  if(is.null(H)){
    H = as.matrix(diag(p))
  }
  df = r*(p-ncol(A))+r*(p-ncol(H))

  out = johansenCpp(Y, Z, r, as.matrix(A), as.matrix(H), dt, intercept, normalize)

  if(r > 0){
    a.hat = matrix(out[1:r,],nr=p,byrow = TRUE)
    b.hat = matrix(out[(r+1):(2*r),],nr=p,byrow = TRUE)
    rownames(a.hat) = paste0("x",1:p)
    colnames(a.hat) = paste0("r",1:r)
    rownames(b.hat) = paste0("x",1:p)
    colnames(b.hat) = paste0("r",1:r)
  }
  P.hat = out[2*r+1,]
  O.hat = out[(2*r+2):(2*r+1+p),]
  test  = out[(2*r+2+p),]
  eigs  = out[(2*r+2+p+1),]
  res   = out[(2*r+2+p+2):(2*r+2+p+N+1),]

  meanY <- apply(Y, 2, mean)

  res0  = Y - matrix(meanY, nrow = N, ncol = p, byrow = T)
  R2    = 1 - sum(res^2)/sum(res0^2)

  names(P.hat)    = paste0("x",1:p)
  rownames(O.hat) = paste0("x",1:p)
  colnames(O.hat) = paste0("x",1:p)
  names(test)     = paste0("r=",0:(p-1))
  names(eigs)     = paste0("l",1:p)
  colnames(res)   = paste0("x",1:p)

  if(r==0){
    a.hat = b.hat = rep(0, p)
  }
  return(list(N = N, p = p, r = r, alpha = a.hat, beta = b.hat, Psi = P.hat, Omega = O.hat,
              test = test, lambda = eigs, A = A, H = H, df = df, dt = dt, r2 = R2))
}
