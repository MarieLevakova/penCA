
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <stdio.h>
#include <Rcpp.h>
#include <math.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

// Soft-thresholding operator
double softThres(double x, double lambda) {
  return((x > lambda) ? x - lambda :
           (x < -lambda) ? x + lambda : 0.);
}

// Penalized univariate regression with lasso penalty
//
// @param Y response vector
// @param X covariate matrix
// @param lambda scalar penalty parameter
// @param beta0 initial value of regression coefficient vector
// @param control a list of parameters controlling the fitting process
// @return estimated regression coefficient vector

// [[Rcpp::export]]
arma::vec penreg_Rcpp(arma::vec Y, arma::mat X, double lambda, arma::vec beta0,
                      Rcpp::List control) {

  int p = X.n_cols;
  int n = X.n_rows;

  int maxit = control["maxit"];
  double epsilon = control["epsilon"];

  arma::vec ss = arma::conv_to<arma::vec>::from(sum(X % X));

  // for (int i=0; i<p; i++) X.col(i) /= ss(i);
  arma::vec vlambda = lambda * n * pow(ss, -1.);

  beta0 = beta0;
  arma::vec r = Y - X*beta0;
  arma::vec beta(beta0);

  for (int i=0; i<maxit; i++) {
    arma::vec bchg(p);
    for (int j=0; j<p; j++) {
      double bj = sum(X.col(j) % r)/ss(j) + beta(j);
      bj = softThres(bj, vlambda(j));
      beta(j) = bj;
      bchg(j) = beta(j) - beta0(j);
      if (bchg(j) != 0) r -= bchg(j) * X.col(j);
    }
    if (sqrt(sum(bchg % bchg) / sum(beta0 % beta0)) < epsilon) break;
    beta0 = beta;
  }

  return(beta);
}


// Penalized regression with lasso penalty
//
// @param XY product moment matrix of X and Y (i.e., (1/n) * X'Y)
// @param XX Gram matrix (i.e., (1/n) * X'X)
// @param lambda scalar penalty parameter
// @param beta0 initial value of regression coefficient vector
// @param control a list of parameters controling the fitting process
// @return estimated regression coefficient vector
// [[Rcpp::export]]
arma::vec penreg_Rcpp_XY(arma::vec XY, arma::mat XX, double lambda, arma::vec beta0,
                         Rcpp::List control) {
  int p = XX.n_cols;
  arma::vec vlambda = lambda / diagvec(XX);

  int maxit = control["maxit"];
  double epsilon = control["epsilon"];

  arma::vec beta(beta0);
  for (int i=0; i<maxit; i++) {
    arma::vec bchg(p);
    for (int j=0; j<p; j++) {
      double bj = (XY(j) - sum(XX.col(j) % beta))/XX(j,j) + beta(j);
      bj = softThres(bj, vlambda(j));
      beta(j) = bj;
      bchg(j) = beta(j) - beta0(j);
    }
    if (sqrt(sum(bchg % bchg) / sum(beta0 % beta0)) < epsilon) break;
    beta0 = beta;
  }
  return(beta);
}
