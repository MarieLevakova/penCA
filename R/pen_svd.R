#' Fit VEC model with penalty imposed on SVD of Pi
#'
#' This function loads data matrix X and fits a VEC model,
#' where penalty is imposed on singular values and vectors
#' in the singular vector decomposition of Pi
#'
#' @param Y Matrix of differences of the process
#' @param X Matrix of lagged values of the process
#' @param nrank cointegration rank
#' @param ic Criterion for choosing optimal penalty value, the options are
#' \code{"CV"} (crossvalidated error), \code{"AIC"} (Akaike information criterion),
#' \code{"BIC"} (Bayes information criterion) and \code{"HQ"} (Hannan-Quinn criterion).
#' @param control parameters specifying the iterative algorithm
#' @param n.cv number of repetitions of the crossvalidation procedure
#' @return A list containing the following components:
#' @slot Y matrix of differences of the process
#' @slot X matrix of lagged values of the process
#' @slot U.path the regularization path of the matrix of left single vectors
#' @slot V.path the regularization path of the matrix of right single vectors
#' @slot D.path the regularization path of the singular values
#' @slot ic.path the regularization path of the objective function
#' @slot U the final left singular vectors
#' @slot V the final right singular vectors
#' @slot D the final singular values
#' @slot rank the final cointegration rank
#' @slot PI the final estimate of the matrix Pi
#' @export
pen.svd <- function(Y, X, nrank, ic.type = c("AIC", "BIC", "HQ", "CV"),
                  control = list(), n.cv = 20) {

  control <- do.call("pen.svd.control", control)
  tol <- control$df.tol
  ic <- match.arg(ic.type)

  p <- ncol(X)
  q <- ncol(Y)
  n <- nrow(Y)

  if (n != nrow(X))
    stop("'nrow(X)' has to be equal to 'nrow(Y)'.")

  ini <- NULL
  ranknow <- nrank
  iter.eea <- 1
  Ufit <- matrix(nrow = p, ncol = ranknow, 0)
  Vfit <- matrix(nrow = q, ncol = ranknow, 0)
  Dfit <- rep(0, ranknow)

  while (iter.eea <= control$niter.eea) {
    fiti <-  pen.svd.eea(Y, X, nrank = ranknow, ic.type = ic, ini = ini,
                       control = control, n.cv = n.cv)

    ## Current solution
    dorder <- order(fiti$D, decreasing = TRUE)
    Ufit[, 1:ranknow] <- fiti$U[, dorder]
    Vfit[, 1:ranknow] <- fiti$V[, dorder]
    Dfit[1:ranknow] <- fiti$D[dorder]

    ## New rank
    ranknew <- sum(fiti$D > tol)

    ## if zero solution happens
    if (is.na(ranknew) | ranknew == 0) {
      iter.eea <- control$niter.eea
    } else {
      Ufit <- matrix(Ufit[, 1:ranknew], nrow = p, ncol = ranknew)
      Vfit <- matrix(Vfit[, 1:ranknew], nrow = q, ncol = ranknew)
      Dfit <- Dfit[1:ranknew]
    }

    # if ranknew == 1, then only need to iterate once if ranknow!=ranknew.
    if (ranknew == 1 & ranknow != ranknew) {
      ini <- pen.svd.init(Y, X, nrank = ranknew, control = control)
      ini <- pen.svd.init(Y, X, nrank = ranknew, U0 = Ufit, V0 = Vfit, D0 = Dfit,
                        Uw = ini$U0, Vw = ini$V0, control = control)
      ranknow <- ranknew
      # out of loop in the next step. one more iteration
      iter.eea <- control$niter.eea - 1
    }

    if (ranknew > 1 & iter.eea < control$niter.eea) {
      # update weights
      if (ranknow != ranknew) ini <- pen.svd.init(Y, X, nrank = ranknew, control = control)
      ini <- pen.svd.init(Y, X, nrank = ranknew, U0 = Ufit, V0 = Vfit, D0 = Dfit,
                        Uw = ini$U0, Vw = ini$V0, control = control)
      ranknow <- ranknew
    }

    iter.eea <- iter.eea + 1
  }

  U <- Ufit
  V <- Vfit

  if(ranknew != 0){
    Dmat <- matrix(0, nrow = ranknew, ncol = ranknew)
    diag(Dmat) <- Dfit

    PI <- V %*% Dmat %*% t(U)  # PI is transposed
  } else {
    PI <- matrix(0, nrow = p, ncol = p)
  }

  ret <- list(Y = Y, X = X, U.path = fiti$U.path, V.path = fiti$V.path, D.path = fiti$D.path,
              ic.path = fiti$ic.path, U = U, V = V, D = Dfit, rank = ranknew, PI = PI)
  ret
}


pen.svd.eea <- function(Y, X, nrank, minLambda = 1e-04, ic.type = c("AIC", "BIC", "HQ", "CV"),
                      ini = NULL, control = list(), n.cv){

  control <- do.call("pen.svd.control", control)
  ic <- match.arg(ic.type)
  IC <- match(ic, c("AIC", "BIC", "HQ", "CV"))

  p <- ncol(X)
  q <- ncol(Y)
  n <- nrow(Y)

  if (n != nrow(X))
    stop("'nrow(X)' has to be equal to 'nrow(Y)'.")

  if (is.null(ini)) {
    ini <- pen.svd.init(Y, X, nrank = nrank, control = control)
  }

  ## Information criteria
  ## AIC, BIC, HQ
  ## Default is AIC
  nlambda <- control$nlambda
  IC_path <- array(dim = c(nlambda, nrank, 4))

  ## Solution path
  U_path <- array(dim = c(p, nlambda, nrank))
  V_path <- array(dim = c(q, nlambda, nrank))
  D_path <- matrix(nrow = nlambda, ncol = nrank)

  ## Final solution selected by IC
  U <- matrix(nrow = p, ncol = nrank)
  V <- matrix(nrow = q, ncol = nrank)
  D <- rep(NA, nrank)
  Lam <- rep(NA, nrank)

  for (k in 1:nrank) {
    lamSca <- log(ini$lamMax[k]) + log(c(0.9, control$minLambda))
    lambda <- c(exp(seq(lamSca[1], lamSca[2], length = nlambda - 1)), 0)

    UU <- matrix(ncol = nlambda, nrow = p)
    VV <- matrix(ncol = nlambda, nrow = q)
    DD <- rep(NA, nlambda)

    # initial nonzero values
    Uini <- ini$Uini[, k]
    Vini <- ini$Vini[, k]
    Dini <- 1

    WU <- ini$Wuel[, k]
    WV <- ini$Wvel[, k]
    Yel <- ini$Yel[, , k]
    Xtran <- X %*% diag(WU)

    # Fitting
    for (i in 1:nlambda) {
      fit <- surr_fit_Rcpp(Yel, X, lambda[i], Uini, Vini, WU, WV, Xtran, control,
                           n.cv)
      IC_path[i, k, ] <- fit$ic

      VV[, i] <- fit$V
      UU[, i] <- fit$U
      DD[i] <- fit$D
      if (fit$D != 0 & !is.nan(fit$D)) {
        Uini <- fit$U
        Vini <- fit$V
        Dini <- fit$D
      }
    }

    logSY <- sum(Yel ^ 2)
    # Choose lambda
    minid <- apply(IC_path[, k, ], 2, which.min)

    ## Final solution based on IC
    if (logSY <= IC_path[minid[IC], k, IC]) {
      U[, k] <- 0
      V[, k] <- 0
      D[k] <- 0
      Lam[k] <- ini$lamMax[k]
    } else {
      U[, k] <- UU[, minid[IC]]
      V[, k] <- VV[, minid[IC]]
      D[k] <- DD[minid[IC]]
      Lam[k] <- lambda[minid[IC]]
    }

    U_path[, , k] <- UU
    V_path[, , k] <- VV
    D_path[, k] <- DD
  }

  ret <- list(U.path = U_path, V.path = V_path, D.path = D_path, ic.path = IC_path,
    U = U, V = V, D = D)

  ret
}


### internal functions =========================================================

## Internal function for specifying computation parameters
##
## a list of internal computational parameters controlling optimization
##
## @param maxit maximum number of iterations
## @param epsilon convergence tolerance
## @param innerMaxit maximum number of iterations for inner steps
## @param innerEpsilon convergence tolerance for inner steps
## @param nlambda number of lambda values
## @param adaptive If Ture, use adaptive penalization
## @param gamma0 power parameter for constructing adaptive weights
## @param minLambda multiplicate factor to determine the minimum lambda
## @param niter.eea the number of iterations in the iterative exclusive
##     extraction algorithm
## @param df.tol tolerance
##
## @return a list of computational parameters.
pen.svd.control <- function(maxit = 100L,
                          epsilon = 1e-4,
                          innerMaxit = 50L,
                          innerEpsilon = 1e-4,
                          nlambda = 100,
                          adaptive = TRUE,
                          gamma0 = 2,
                          minLambda = 1e-5,
                          niter.eea = 1,
                          df.tol = 1e-8) {
  maxit <- as.integer(maxit)
  innerMaxit <- as.integer(innerMaxit)
  nlambda <- as.integer(nlambda)
  list(
    maxit = maxit,
    epsilon = epsilon,
    innerMaxit = innerMaxit,
    innerEpsilon = innerEpsilon,
    nlambda = nlambda,
    adaptive = adaptive,
    gamma0 = gamma0,
    minLambda = minLambda,
    niter.eea = niter.eea,
    df.tol = df.tol
  )
}


## Initialization for EEA
pen.svd.init <- function(Y, X, U0 = NULL, V0 = NULL, D0 = NULL,
                       Uw = U0, Vw = V0, nrank = 1, control = list()) {
  p <- ncol(X)
  q <- ncol(Y)
  n <- nrow(Y)

  if (is.null(U0) || is.null(V0) || is.null(D0)) {
    ## Calculate LS estimator and RR estimator
    ini <- rrr.fit(Y, X, nrank = nrank, coefSVD = TRUE)
    C_ls <- ini$coef.ls
    C_rr <- ini$coef
    U0 <- ini$coefSVD$u
    V0 <- ini$coefSVD$v
    D0 <- ini$coefSVD$d
  } else {
    C_rr <- U0 %*% (D0 * t(V0))
    C_ls <- NULL
  }

  U0 <- matrix(U0, nrow = p, ncol = nrank)
  V0 <- matrix(V0, nrow = q, ncol = nrank)

  if (is.null(Uw))
    Uw <- U0
  if (is.null(Vw))
    Vw <- V0

  ## Compute exclusive layers
  Yel <- array(dim = c(n, q, nrank))
  Wuel <- matrix(1, p, nrank)
  Wvel <- matrix(1, q, nrank)

  ## Initial nonzero singular vectors
  Uini <- matrix(0, p, nrank)
  Vini <- matrix(0, q, nrank)

  control <- do.call("pen.svd.control", control)
  ## Maximum lambda
  lamMax <- double(nrank)
  tol <- control$df.tol
  gamma0 <- control$gamma0
  for (k in 1:nrank) {
    Ck <- tcrossprod(U0[, k], V0[, k]) * D0[k]
    Cmk <- C_rr - Ck
    Yel[, , k] <- Y - X %*% Cmk

    if (control$adaptive) {
      if (all(Uw[, k] == 0)) {
        Wuel[, k] <- 1
      } else {
        Wuel[, k] <- abs(Uw[, k]) ^ gamma0
        Wuel[, k][Wuel[, k] == 0] <- tol
      }
      if (all(Vw[, k] == 0)) {
        Wvel[, k] <- 1
      } else {
        Wvel[, k] <- abs(Vw[, k]) ^ gamma0
        Wvel[, k][Wvel[, k] == 0] <- tol
      }
    }

    ww <- outer(Wuel[, k], Wvel[, k])
    xyk <-  crossprod(X, Yel[, , k])
    lam <- abs(ww * xyk)

    imax <- which.max(lam)
    b <- round(imax / p - 0.50001) + 1
    a <- imax - (b - 1) * p
    lamMax[k] <- lam[a, b]
    Uini[as.integer(a), k] <- 1
    Vini[as.integer(b), k] <- 1
  }

  list(
    Yel = Yel,
    C_ls = C_ls,
    C_rr = C_rr,
    Wuel = Wuel,
    Wvel = Wvel,
    U0 = matrix(U0, ncol = nrank),
    V0 = matrix(V0, ncol = nrank),
    D0 = D0,
    Uini = Uini,
    Vini = Vini,
    lamMax = lamMax,
    Mdim = c(p, q, n)
  )
}

## Fitting reduced-rank regression with a specific rank
##
## Given a response matrix and a covariate matrix, this function fits reduced
## rank regression for a specified rank. It reduces to singular value
## decomposition if the covariate matrix is the identity matrix.
##
## @usage
## rrr.fit(Y, X, nrank = 1, weight = NULL, coefSVD = FALSE)
##
## @param Y a matrix of response (n by q)
## @param X a matrix of covariate (n by p)
## @param nrank an integer specifying the desired rank
## @param weight a square matrix of weight (q by q); The default is the
##     identity matrix
## @param coefSVD logical indicating the need for SVD for the coefficient matrix
##     in the output; used in ssvd estimation
## @return S3 \code{rrr} object, a list consisting of \item{coef}{coefficient
##     of rrr} \item{coef.ls}{coefficient of least square} \item{fitted}{fitted
##     value of rrr} \item{fitted.ls}{fitted value of least square}
##     \item{A}{right singular matrix} \item{Ad}{a vector of sigular values}
##     \item{rank}{rank of the fitted rrr}
## @examples
## Y <- matrix(rnorm(400), 100, 4)
## X <- matrix(rnorm(800), 100, 8)
## rfit <- rrr.fit(Y, X, nrank = 2)
## coef(rfit)
## @importFrom MASS ginv


rrr.fit <- function(Y, X, nrank = 1, coefSVD = FALSE){

  q <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)
  stopifnot(n == nrow(X))

  S_yx <- crossprod(Y, X)
  S_xx <- crossprod(X)

  ## FIXME: 0.01 is too arbitrary
  S_xx_inv <- tryCatch(
    ginv(S_xx),
    error = function(e)
      solve(S_xx + 0.01 * diag(p))
  )

  C_ls <- tcrossprod(S_xx_inv, S_yx)

  XC <- X %*% C_ls
  svdXC <- svd(XC, nrank, nrank)
  A <- svdXC$v[, 1:nrank]
  Ad <- (svdXC$d[1:nrank]) ^ 2
  AA <- tcrossprod(A)
  C_rr <- C_ls %*% AA

  ret <- list(
    coef = C_rr,
    coef.ls = C_ls,
    fitted = X %*% C_rr,
    fitted.ls = XC,
    A = A,
    Ad = Ad,
    rank = nrank
  )

  if (coefSVD) {
    coefSVD <- svd(C_rr, nrank, nrank)
    coefSVD$d <- coefSVD$d[1:nrank]
    coefSVD$u <- coefSVD$u[, 1:nrank, drop = FALSE]
    coefSVD$v <- coefSVD$v[, 1:nrank, drop = FALSE]
    ret <- c(ret, list(coefSVD = coefSVD))
  }

  ret
}
