// This is start of the header guard.  ADD_H can be any unique name.  By convention, we use the name of the header file.
#ifndef johansen_H
#define johansen_H

using namespace arma;
// This is the content of the .h file, which is where the declarations go
arma::mat vecmAggregated(arma::mat Z0, arma::mat Z1, int r, arma::mat A, arma::mat B, double dt, bool intercept, bool normalize);
arma::mat varAggregated(arma::mat Z0, arma:: mat Z1, double dt, bool intercept, bool normalize);
arma::mat johansenCppAggregated(arma::mat Z0, arma::mat Z1, int r, arma::mat A, arma::mat B, double dt, bool intercept, bool normalize);


// This is the end of the header guard
#endif
