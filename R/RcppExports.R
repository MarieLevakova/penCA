# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

vecm <- function(Z0, Z1, r, A, B, dt, intercept, normalize) {
    .Call('_penCA_vecm', PACKAGE = 'penCA', Z0, Z1, r, A, B, dt, intercept, normalize)
}

johansenCpp <- function(Z0, Z1, r, A, B, dt, intercept, normalize) {
    .Call('_penCA_johansenCpp', PACKAGE = 'penCA', Z0, Z1, r, A, B, dt, intercept, normalize)
}

penreg_Rcpp <- function(Y, X, lambda, beta0, control) {
    .Call('_penCA_penreg_Rcpp', PACKAGE = 'penCA', Y, X, lambda, beta0, control)
}

penreg_Rcpp_XY <- function(XY, XX, lambda, beta0, control) {
    .Call('_penCA_penreg_Rcpp_XY', PACKAGE = 'penCA', XY, XX, lambda, beta0, control)
}

penAdaptNuclearCpp <- function(X, n_lambda, lambda_min, crit, dt, n_cv, w_gamma, alpha) {
    .Call('_penCA_penAdaptNuclearCpp', PACKAGE = 'penCA', X, n_lambda, lambda_min, crit, dt, n_cv, w_gamma, alpha)
}

penNuclearCpp <- function(X, n_lambda, lambda_min, miniter, maxiter, crit, dt, n_cv, thresh, alpha) {
    .Call('_penCA_penNuclearCpp', PACKAGE = 'penCA', X, n_lambda, lambda_min, miniter, maxiter, crit, dt, n_cv, thresh, alpha)
}

penRankCpp <- function(X, n_lambda, lambda_min, crit, dt, n_cv, alpha) {
    .Call('_penCA_penRankCpp', PACKAGE = 'penCA', X, n_lambda, lambda_min, crit, dt, n_cv, alpha)
}

accu2 <- function(obj) {
    .Call('_penCA_accu2', PACKAGE = 'penCA', obj)
}

penPiCpp <- function(X, n_lambda, lambda_min, r, maxiter, crit, dt, w_auto, n_cv, q, weights, lambda_max) {
    .Call('_penCA_penPiCpp', PACKAGE = 'penCA', X, n_lambda, lambda_min, r, maxiter, crit, dt, w_auto, n_cv, q, weights, lambda_max)
}

surr_fit_Rcpp <- function(Y, X, lambda, U0, V0, WU, WV, Xtran, control, n_cv) {
    .Call('_penCA_surr_fit_Rcpp', PACKAGE = 'penCA', Y, X, lambda, U0, V0, WU, WV, Xtran, control, n_cv)
}

