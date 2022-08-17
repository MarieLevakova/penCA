#' Fit VEC model with nuclear-norm penalty
#'
#' This function loads data matrix X and fits a VEC model,
#' where elements of the matrix Pi are penalized with nuclear norm penalty.
#'
#' @param X Data matrix
#' @param n.lambda number of penalty values lambda to be used
#' @param lambda.min the smallest allowed penalty value
#' @param miniter minimum number of iterations
#' @param maxiter maximum number of iterations
#' @param dt timestep
#' @param crit Criterion for choosing optimal penalty value, the options are
#' \code{"CV"} (crossvalidated error), \code{"AIC"} (Akaike information criterion),
#' \code{"BIC"} (Bayes information criterion) and \code{"HQ"} (Hannan-Quinn criterion).
#' @param n.cv number of repetitions of the crossvalidation procedure
#' @param thresh convergence parameter
#' @param alpha proportion of the ridge penalty
#' @return A list containing the following components:
#' @slot PI estimated matrix Pi, recalculated per time unit (Pi/\code{dt})
#' @slot MU estimated intercepts, recalculated per time unit (mu/\code{dt})
#' @slot OMEGA the covariance matrix of the random errors, recalculated per time unit (Omega/\code{dt})
#' @slot lambda the selected penalty value
#' @slot PI.iter estimates of Pi throughout iterations of the algorithm, each row contains a vectorized form of matrix Pi
#' @slot objective.iter vector of values of the objective function throughout iterations
#' @slot lambda.seq the sequence of penalty values
#' @slot crit the type of criterion to select optimal lambda
#' @slot crit.seq the sequence of values of the criterion throughout iterations
#' @slot Pi.seq the sequence of Pi estimates for each penalty value in lambda.seq (not divided by \code{dt})
#' @slot iter.seq number of iterations for each penalty value in lambda.seq (not divided by \code{dt})
#' @export
pen.PI.nuclear <- function(X, n.lambda, lambda.min = 1e-12,
                           miniter = 10, maxiter = 100,
                           dt = 1, crit = "CV", n.cv = 100, thresh = 1e-6, alpha = 0.01){

  p <- dim(X)[2]

  # Translating crit from a string to an integer value
  crit.vector <- c("CV", "AIC", "BIC", "HQ")
  crit.int <- which(crit.vector ==  crit) - 1

  # Calling cpp function
  out <- penNuclearCpp(X, n.lambda, lambda.min, miniter, maxiter, crit.int,
                       dt, n.cv, thresh, alpha)

  # Extracting results from the cpp output
  PI <- out[1:p,1:p]
  MU <- out[1:p,p+1]
  OMEGA <- out[1:p,(p+2):(2*p+1)]
  lambda <- out[1,2*p+2]
  PI.iter <- out[(p+1):(p+maxiter-1),1:(p^2)]
  objective.iter <- out[(p+1):(p+maxiter-1),p^2+1]
  lambda.seq <- out[(p+maxiter+1):(p+maxiter+n.lambda),1]
  crit.seq <- out[(p+maxiter+1):(p+maxiter+n.lambda),2]
  Pi.seq <- out[(p+maxiter+1):(p+maxiter+n.lambda),3:(p^2+2)]
  iter.seq <- out[(p+maxiter+1):(p+maxiter+n.lambda),p^2+3]

  return(list(PI = PI, MU = MU, OMEGA = OMEGA, lambda = lambda,
              PI.iter = PI.iter, objective.iter = objective.iter,
              lambda.seq = rev(lambda.seq), crit = crit, crit.seq = rev(crit.seq),
              Pi.seq = Pi.seq, iter.seq = iter.seq))
}
