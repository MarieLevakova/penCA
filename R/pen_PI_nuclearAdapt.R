#' Fit VEC model with adaptive nuclear-norm penalty
#'
#' This function loads data matrix X and fits a VEC model,
#' where elements of the matrix Pi are penalized with adaptive nuclear norm penalty.
#'
#' @param X Data matrix
#' @param n.lambda number of penalty values lambda to be used
#' @param lambda.min the smallest allowed penalty value
#' @param dt timestep
#' @param crit Criterion for choosing optimal penalty value, the options are
#' \code{"CV"} (crossvalidated error), \code{"AIC"} (Akaike information criterion),
#' \code{"BIC"} (Bayes information criterion) and \code{"HQ"} (Hannan-Quinn criterion).
#' @param n.cv number of repetitions of the crossvalidation procedure
#' @param w.gamma parameter of the weights
#' @param alpha proportion of the ridge penalty
#' @return A list containing the following components:
#' @slot PI estimated matrix Pi, recalculated per time unit (Pi/\code{dt})
#' @slot MU estimated intercepts, recalculated per time unit (mu/\code{dt})
#' @slot OMEGA the covariance matrix of the random errors, recalculated per time unit (Omega/\code{dt})
#' @slot lambda the selected penalty value
#' @slot lambda.seq the sequence of penalty values
#' @slot crit the type of criterion to select optimal lambda
#' @slot crit.seq the sequence of values of the criterion throughout iterations
#' @slot Pi.seq the sequence of Pi estimates for each penalty value in lambda.seq (not divided by \code{dt})
#' @export
pen.PI.nuclearAdapt <- function(X, n.lambda, lambda.min = 1e-12, dt = 1, crit = "CV",
                                n.cv = 100, w.gamma = 0, alpha = 0.01){

  p <- dim(X)[2]

  # Translating crit from a string to an integer value
  crit.vector <- c("CV", "AIC", "BIC", "HQ")
  crit.int <- which(crit.vector ==  crit) - 1

  # Calling cpp function
  out <- penAdaptNuclearCpp(X, n.lambda, lambda.min, crit.int, dt, n.cv, w.gamma, alpha)

  # Extracting results from the cpp output
  PI <- out[1:p,1:p]
  MU <- out[1:p,p+1]
  OMEGA <- out[1:p,(p+2):(2*p+1)]
  lambda <- out[1,2*p+2]
  lambda.seq <- out[(p+1):(p+n.lambda), 1]
  crit.seq <- out[(p+1):(p+n.lambda), 2]
  PI.seq <- out[(p+n.lambda):(p+1), 3:(p^2+2)] # reverse the order

  return(list(PI = PI, MU = MU, OMEGA = OMEGA, lambda = lambda,
              lambda.seq = rev(lambda.seq), crit = crit, crit.seq = rev(crit.seq),
              PI.seq = PI.seq))
}
