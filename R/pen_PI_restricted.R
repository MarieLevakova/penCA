#' Fit VEC model with lasso penalty and restricted rank
#'
#' This function loads data matrix X and fits a VEC model,
#' where elements of the matrix Pi are penalized with lasso penalty
#' while at the same time the rank of Pi is restricted to a given value.
#'
#' @param X Data matrix
#' @param n.lambda number of penalty values lambda to be used
#' @param lambda.min the smallest allowed penalty value
#' @param lambda.max the largest allowed penalty value (not used, calculated automatically)
#' @param r rank of Pi
#' @param maxiter maximum number of iterations
#' @param dt timestep of the process
#' @param crit Criterion for choosing optimal penalty value, the options are
#' \code{"CV"} (crossvalidated error), \code{"AIC"} (Akaike information criterion),
#' \code{"BIC"} (Bayes information criterion) and \code{"HQ"} (Hannan-Quinn criterion).
#' @param w.auto logical value, whether the sequence of timesteps alpha_k of the algorithm should be determined by an in-built procedure
#' @param n.cv number of repetitions of the crossvalidation procedure
#' @param q timestep parameter when timesteps of the algorithm are chosen automatically
#' @param weights pxp matrix of weights (default: identical weights)
#' @return A list containing the following components:
#' @slot PI estimated matrix Pi, recalculated per time unit (Pi/\code{dt})
#' @slot MU estimated intercepts, recalculated per time unit (mu/\code{dt})
#' @slot OMEGA the covariance matrix of the random errors, recalculated per time unit (Omega/\code{dt})
#' @slot lambda the selected penalty value
#' @slot PI.iter matrix with estimates of Pi throughout iterations of the algorithm, each row contains vectorized Pi
#' @slot objective.iter values of the objective function throughout iterations
#' @slot w the sequence of timesteps of the algorithm
#' @slot lambda.seq the sequence of penalty values
#' @slot crit the type of criterion to select optimal lambda
#' @slot crit.seq the sequence of values of the criterion throughout iterations
#' @slot Pi.seq the sequence of Pi estimates for each penalty value in lambda.seq (not divided by \code{dt})
#' @export
pen.PI.restricted <- function(X, n.lambda, lambda.min = 1e-12, lambda.max = 0, r,
                              maxiter = 100,
                              dt = 1, crit = "CV", w.auto = TRUE, n.cv = 100, q = 2,
                              weights = matrix(1, nrow = ncol(X), ncol = ncol(X))){

  p <- dim(X)[2]

  # Translating crit from a string to an integer value
  crit.vector <- c("CV", "AIC", "BIC", "HQ")
  crit.int <- which(crit.vector ==  crit) - 1

  # Calling cpp function
  out <- penPiCpp(X, n.lambda, lambda.min, r, maxiter, crit.int, dt,
                  w.auto, n.cv, q, t(weights), lambda.max)

  # Extracting results from the cpp output
  PI <- out[1:p,1:p]
  MU <- out[1:p,p+1]
  OMEGA <- out[1:p,(p+2):(2*p+1)]
  lambda <- out[1,2*p+2]
  PI.iter <- out[(p+1):(p+maxiter-1),1:(p^2)]
  objective.iter <- out[(p+1):(p+maxiter-1),p^2+1]
  w <- out[(p+1):(p+maxiter),p^2+2]
  lambda.seq <- out[(p+maxiter+1):(p+maxiter+n.lambda),1]
  crit.seq <- out[(p+maxiter+1):(p+maxiter+n.lambda),2]
  Pi.seq <- out[(p+maxiter+1):(p+maxiter+n.lambda),3:(p^2+2)]

  return(list(PI = PI, MU = MU, OMEGA = OMEGA, lambda = lambda,
              PI.iter = PI.iter, objective.iter = objective.iter, w = w,
              lambda.seq = lambda.seq, crit = crit, crit.seq = rev(crit.seq),
              Pi.seq = Pi.seq))
}
