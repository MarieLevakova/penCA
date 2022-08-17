// Estimates the matrix Pi of a restricted rank r in a VECM model.
// Algorithm is from Lian, Zhao:
// Parametric and semiparametric reduced-rank regression with flexible sparsity (2015)

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>     // Cpp linear algebra library
#include <stdio.h>
#include <math.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

arma::mat penNuclearLoop(arma::mat Ystd, arma::mat Zstd, arma::mat Pi_init,
                         double lambda, int miniter, int maxiter, double thresh,
                         double alpha){

  int N = Ystd.n_rows;
  int p = Ystd.n_cols;

  Ystd = join_cols(Ystd, zeros<mat>(p,p));
  Zstd = join_cols(Zstd, sqrt(lambda*(1-alpha)/2)*eye<mat>(p,p));

  vec t_k = zeros<vec>(maxiter);
  vec r_k = zeros<vec>(maxiter);
  mat Pi_restricted = Pi_init;
  mat Pi_star = Pi_init;
  mat Pi_previous = Pi_init;

  t_k(0) = 1;
  double L = norm(Zstd.t()*Zstd, "fro");

  mat PI_iter = zeros<mat>(maxiter-1, p*p);
  vec objective_iter = zeros<vec>(maxiter-1);
  mat C_k = zeros<mat>(p,p);
  int ii = 2;

  do {
    C_k = Pi_star - (1/L)*(Zstd*Pi_star.t()-Ystd).t() * Zstd;

    mat svd_U;
    vec svd_d;
    mat svd_V;
    arma::svd(svd_U, svd_d, svd_V, C_k);

    vec d_soft = svd_d;
    for(int jj = 0; jj < p; jj++){
      if(d_soft(jj)>lambda*alpha/L) {
        d_soft(jj) = d_soft(jj) - lambda*alpha/L;
      } else {
        d_soft(jj) = 0;
      }
    }
    Pi_previous = Pi_restricted;
    Pi_restricted = svd_U * diagmat(d_soft) * svd_V.t();
    t_k(ii) = (1+sqrt(1 + 4*pow(t_k(ii-1),2)))/2;
    r_k(ii) = (t_k(ii-1)-1)/t_k(ii);

    Pi_star = Pi_restricted + r_k(ii)*(Pi_restricted - Pi_previous);

    PI_iter.row(ii-1) = reshape(Pi_restricted, 1, p*p);
    mat res = Ystd - Zstd * Pi_restricted.t();
    objective_iter[ii-1] = pow(norm(res, "fro"), 2)/(2*N) + lambda*norm(Pi_restricted, 1); //

    ii = ii + 1;
  } while ((ii < maxiter && abs(objective_iter[ii-2]-objective_iter[ii-3])>thresh) || ii < miniter); //

  // Assemble the output
  mat out = zeros<mat>(p+maxiter, p*p+1);

  // Insert PI in the output
  out(span(0, p-1), span(0, p-1)) = Pi_restricted;

  // Insert n_iter in the output
  out(0, p) = ii;

  // Insert PI_iter in the output
  out(span(p, p+maxiter-2), span(0, p*p-1)) = PI_iter;

  // Insert objective_iter in the output
  out(span(p, p+maxiter-2), p*p) = objective_iter;

  // return out;
  return out;
}

// Function to estimate Pi using nuclear norm penalty,
// lambda chosen by cross-validation method based on average error of 1-step forecast

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

arma::mat penNuclearCpp(arma::mat X, int n_lambda, double lambda_min,
                   int miniter, int maxiter, int crit, double dt, int n_cv,
                   double thresh, double alpha){

  int N = X.n_rows-1;   // Number of observations
  int p = X.n_cols;     // Dimension of the system
  mat Y = zeros<mat>(N,p);  // Matrix of zeros with p rows and N columns
  mat Z = zeros<mat>(N,p);
  for(int n=0;n<N;n++){
    Y.row(n) = X.row(n+1)-X.row(n); // Fills rows of Y with differences of X.
    Z.row(n) = X.row(n);            // Fills rows of Z with levels of X.
  }

  // Standardize variables
  mat meanY = mean(Y); // Column means of Y
  mat meanZ = mean(Z); // Column means of Z

  mat sdY = stddev(Y); // Column standard deviations of Y
  mat sdZ = stddev(Z); // Column standard deviations of Z

  mat Ystd = (Y - ones<mat>(N,1)*meanY); // *diagmat(1/sdY);
  mat Zstd = (Z - ones<mat>(N,1)*meanZ); // *diagmat(1/sdZ);

  // Calculate the sequence of lambdas
  double lambda_max = norm(Ystd.t() * Zstd, "fro")/alpha;
  vec lambda_seq = linspace(lambda_min, lambda_max, n_lambda);
  vec crit_value = zeros<vec>(n_lambda);

  // Choose optimal lambda
  double lambda;
  double lambda_opt;
  mat Pi_restricted;
  mat Pi_init = (Zstd.t()*Zstd).i()*Zstd.t()*Ystd; //zeros<mat>(p,p);

  mat  pen_out;
  if(crit != 0){ //lambda not chosen by crossvalidation
    for(int i=0; i<n_lambda; i++){
      lambda = lambda_seq(i);

      pen_out = penNuclearLoop(Ystd, Zstd, Pi_init, lambda, miniter, maxiter,
                               thresh, alpha);
      Pi_restricted = pen_out(span(0, p-1), span(0, p-1));
      Pi_init = Pi_restricted;

      int k = accu(conv_to<imat>::from(Pi_restricted!=zeros<mat>(p,p)));
      mat res = Ystd - Zstd * Pi_restricted.t();
      mat Omega_select = (res.t() * res)/N;
      mat Omega_inv = Omega_select.i();

      double logdet_Omega;
      double sign;

      log_det(logdet_Omega, sign, Omega_select);

      if(crit==1) { // AIC
        crit_value(i) = N*p*log(2*datum::pi) + N*logdet_Omega + 2*k + trace(res*Omega_inv*res.t());
      } else if(crit==2) { // BIC
        crit_value(i) = N*p*log(2*datum::pi) + N*logdet_Omega + k*log(N) + trace(res*Omega_inv*res.t());
      } else { // HQ
        crit_value(i) = N*p*log(2*datum::pi) + N*logdet_Omega + 2*k*log(log(N)) + trace(res*Omega_inv*res.t());
      }
    }
    lambda_opt = lambda_seq(crit_value.index_min());
  } else { // CV
    vec cv = zeros<vec>(n_lambda);

    // Run crossvalidation n_cv times
    for(int cv_run=0; cv_run<n_cv; cv_run++){
      // Divide data into 5 folds
      ivec folds = randi<ivec>(N, distr_param(1, 5));
      for(int ii=0; ii<5; ii++){
        mat Ystd_cv = Ystd.rows(find(folds!=ii));
        mat Zstd_cv = Zstd.rows(find(folds!=ii));
        Pi_init = (Zstd_cv.t()*Zstd_cv).i()*Zstd_cv.t()*Ystd_cv; //zeros<mat>(p,p); //
        for(int i=0; i<n_lambda; i++){
          lambda = lambda_seq(i);
          pen_out = penNuclearLoop(Ystd_cv, Zstd_cv, Pi_init, lambda,
                                   miniter, maxiter, thresh, alpha);
          Pi_restricted = pen_out(span(0,p-1), span(0,p-1));
          Pi_init = Pi_restricted;
          mat res = Ystd.rows(find(folds==ii)) - Zstd.rows(find(folds==ii))*Pi_restricted.t();
          cv(i) = cv(i) + trace(res*res.t())/(N*p*n_cv);
        }
      }
    }

    lambda_opt = lambda_seq(cv.index_min());
    crit_value = cv;
  }

  // Fit with an optimal lambda
  // (Go through the whole sequence of lambdas to achieve warm starts)
  Pi_init = (Zstd.t()*Zstd).i()*Zstd.t()*Ystd;//zeros<mat>(p,p); //
  mat PI_iter;
  mat objective_iter;
  mat Pi_lambda = zeros<mat>(n_lambda, p*p);
  vec iter_lambda = zeros<vec>(n_lambda);
  mat Pi_final = zeros<mat>(p,p);
  for(int i=0; i<n_lambda; i++){
    lambda = lambda_seq(i);

    pen_out = penNuclearLoop(Ystd, Zstd, Pi_init, lambda, miniter, maxiter,
                             thresh, alpha);
    Pi_restricted = pen_out(span(0, p-1), span(0, p-1));
    Pi_init = Pi_restricted;
    Pi_lambda.row(i) = reshape(Pi_restricted, 1, p*p);
    iter_lambda(n_lambda-i-1) = pen_out(0,p);
    if(lambda == lambda_opt){
      Pi_final = Pi_restricted;
      PI_iter = pen_out(span(p, p+maxiter-2), span(0, p*p-1));
      objective_iter = pen_out(span(p, p+maxiter-2), p*p);
    }
  }

  // Final unnormalization
  // for(int ir=0; ir<p; ir++){
  //   Pi_restricted.row(ir) = Pi_restricted.row(ir)*diagmat(sdY(ir)/sdZ);
  // }
  mat mu_hat = meanY.t() - Pi_final * meanZ.t();
  mat res = Y - ones<mat>(N, 1)*mu_hat.t() - Z*Pi_final.t();
  mat Omega_hat = (res.t()*res)/N;

  // Assemble the output
  mat mat_output = zeros<mat>(p+maxiter+n_lambda,p*p+3);

  mat_output(span(0, p-1), span(0, p-1)) = Pi_final/dt;
  mat_output(span(0, p-1), p) = mu_hat/dt;
  mat_output(span(0, p-1), span(p+1, 2*p)) = Omega_hat/dt;
  mat_output(0, 2*p+1) = lambda_opt;
  mat_output(span(p, p+maxiter-2), span(0, p*p-1)) = PI_iter;
  mat_output(span(p, p+maxiter-2), p*p) = objective_iter;
  mat_output(span(p+maxiter, p+maxiter+n_lambda-1), 0) = lambda_seq;
  mat_output(span(p+maxiter, p+maxiter+n_lambda-1), 1) = crit_value;
  mat_output(span(p+maxiter, p+maxiter+n_lambda-1), span(2, p*p+1)) = Pi_lambda;
  mat_output(span(p+maxiter, p+maxiter+n_lambda-1), p*p+2) = iter_lambda;

  return mat_output;
}
