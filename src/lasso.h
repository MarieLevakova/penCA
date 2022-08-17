// This is start of the header guard.  ADD_H can be any unique name.  By convention, we use the name of the header file.
#ifndef lasso_H
#define lasso_H

using namespace arma;
// This is the content of the .h file, which is where the declarations go
arma::vec penreg_Rcpp(arma::vec Y, arma::mat X, double lambda, arma::vec beta0, Rcpp::List control);
arma::vec penreg_Rcpp_XY(arma::vec XY, arma::mat XX, double lambda, arma::vec beta0, Rcpp::List control);

// This is the end of the header guard
#endif
