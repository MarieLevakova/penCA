// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <stdio.h>
#include <Rcpp.h>
#include <math.h>
#include "lasso.h"

using namespace arma;
using namespace Rcpp;
using namespace std;

// Sparse unit-rank regression, Rcpp version

// [[Rcpp::export]]
Rcpp::List surr_fit_Rcpp(arma::mat Y, arma::mat X, double lambda,
                         arma::vec U0, arma::vec V0,
                         arma::vec WU, arma::vec WV,
                         arma::mat Xtran, Rcpp::List control,
                         int n_cv){

  int p = X.n_cols;
  int q = Y.n_cols;
  int n = Y.n_rows;

  arma::vec y = vectorise(Y);
  Rcpp::List out;

  arma::mat XXtran = Xtran.t()*Xtran;
  arma::mat XYtran = Xtran.t()*Y;

  int maxit = control["maxit"];
  double epsilon = control["epsilon"];
  Rcpp::List innerControl;
  innerControl["epsilon"] = control["innerEpsilon"];
  innerControl["maxit"] = control["innerMaxit"];

  // Matrices used in iterations
  arma::vec U = U0;
  arma::vec V = V0;
  double D = 0.0;

  arma::vec Up;
  arma::vec Vp;

  arma::mat Xu;
  arma::mat Xv;

  // Quantities used in iterations
  double norm_UV;
  double conv = 0.0;
  double sse, AIC, BIC, HQ;
  arma::vec res;
  bool converged;
  double dfu, dfv;
  bool flag = FALSE; // Logical value indicating if C_k is null

  // Repeat the same procedure for crossvalidation
  double CV = 0;
  for(int cv_run=0; cv_run<n_cv; cv_run++){ // Loop over repetitions

    // Divide data into 5 folds
    ivec folds = randi<ivec>(n, distr_param(1, 5));

    for(int ii=0; ii<5; ii++){ // Loop over folds

      arma::mat Xtran_i = Xtran.rows(find(folds!=ii));
      arma::mat Xtran_ii = Xtran.rows(find(folds==ii));
      arma::mat Y_i = Y.rows(find(folds!=ii));
      arma::mat XXtran_i = Xtran_i.t()*Xtran_i;
      arma::mat XYtran_i = Xtran_i.t()*Y_i;
      arma::vec y_i = vectorise(Y.rows(find(folds==ii)));

      arma::mat ww = (arma::conv_to<mat>::from(WU)) * (arma::conv_to<mat>::from(WV)).t();
      arma::mat xyk = (X.rows(find(folds!=ii))).t() * Y.rows(find(folds!=ii));
      double lam = (abs(ww % xyk)).max();

      U = U0;
      V = V0;
      D = 0.0;
      flag = FALSE;

      if(lam<lambda){

        U = U/WU;
        V = V/WV;
        U = U/norm(U,1);
        V = V/norm(V,1);

        for(int j = 1; j < maxit; j++){
          Up = U;
          Vp = V;
          norm_UV = norm(U*V.t()*as_scalar(D), 2); // Norm of C_k

          // Solve for (d,U) while fixing V
          U = penreg_Rcpp_XY(XYtran_i*(WV%V), XXtran_i*pow(norm(WV%V,2),2), lambda, Up, innerControl);
          D = norm(U,1);
          if (D < std::numeric_limits<double>::epsilon() && D > -std::numeric_limits<double>::epsilon()){
            D = 0.0;
            flag = TRUE;
            conv = epsilon;
            conv *= 2;
            break;
          }
          U = U/as_scalar(D);

          // Solve for (d,V) while fixing U
          Xv = kron(diagmat(WV), Xtran_ii*U);
          V = penreg_Rcpp_XY(WV%(XYtran_i.t()*U), diagmat(WV%WV*as_scalar(pow(norm(Xtran_i*U,2),2))),
                             lambda, Vp, innerControl);
          D = norm(V,1);
          if (D <  std::numeric_limits<double>::epsilon() && D > -std::numeric_limits<double>::epsilon()){
            D = 0.0;
            flag = TRUE;
            conv = epsilon;
            conv *= 2;
            break;
          }
          V = V/as_scalar(D);

          conv = norm(U*V.t() - Up*Vp.t(),2)/norm_UV;
          if(conv <= epsilon) break;
        }
      } else {
        flag = TRUE;
      }

      if(flag == TRUE){
        res = y_i;
      } else {
        res = y_i - Xv*V*as_scalar(D);
      }
      sse = pow(norm(res,2),2);
      CV += sse/(n_cv*5);

    } // End of loop over the folds
  } // End of loop over the CV runs

  // End of crossvalidation

  flag = FALSE;

  U = U0;
  V = V0;
  D = 0.0;

  U = U/WU;
  V = V/WV;
  U = U/norm(U,1);
  V = V/norm(V,1);

  for(int j = 1; j < maxit; j++){
    Up = U;
    Vp = V;
    norm_UV = norm(U*V.t()*as_scalar(D), 2); // Norm of C_k

    // Solve for (d,U) while fixing V
    U = penreg_Rcpp_XY(XYtran*(WV%V), XXtran*pow(norm(WV%V,2),2), lambda, Up, innerControl);
    D = norm(U,1);
    if (D < std::numeric_limits<double>::epsilon() && D > -std::numeric_limits<double>::epsilon()){
      D = 0.0;
      flag = TRUE;
      conv = epsilon;
      conv *= 2;
      break;
    }
    U = U/as_scalar(D);

    // Solve for (d,V) while fixing U
    Xv = kron(diagmat(WV), Xtran*U);
    V = penreg_Rcpp_XY(WV%(XYtran.t()*U), diagmat(WV%WV*as_scalar(pow(norm(Xtran*U,2),2))),
                       lambda, Vp, innerControl);
    D = norm(V,1);
    if (D <  std::numeric_limits<double>::epsilon() && D > -std::numeric_limits<double>::epsilon()){
      D = 0.0;
      flag = TRUE;
      conv = epsilon;
      conv *= 2;
      break;
    }
    V = V/as_scalar(D);

    conv = norm(U*V.t() - Up*Vp.t(),2)/norm_UV;
    if(conv <= epsilon) break;
  }

  if(flag == TRUE){
    res = y;
    sse = pow(norm(y,2),2);
    V = V.zeros();
    U = U.zeros();
    D = 0.0;
    dfu = 0.0;
    dfv = 0.0;
    AIC = sse;
    BIC = sse;
    HQ = sse;
  } else {
    res = y - Xv*V*as_scalar(D);
    sse = pow(norm(res,2),2);
    dfu =  accu(U != 0);
    dfv =  accu(V != 0);
    V = WV%V*as_scalar(D);
    U = WU%U;
    double dv = norm(V,2);
    double du = norm(U,2);
    D = dv*du;
    U = U/du;
    V = V/dv;
    AIC = sse + 2*(dfu+dfv-1);
    BIC = sse + std::log(static_cast<double>(q*n))*(dfu+dfv-1);
    HQ = sse + std::log(std::log(static_cast<double>(q*n)))*(dfu+dfv-1);
  }

  if(conv <= epsilon){
    converged = TRUE;
  } else {
    converged = FALSE;
  }
  arma::vec df(2);
  df(0) = dfu;
  df(1) = dfv;
  arma::vec ic(4);
  ic(0) = AIC;
  ic(1) = BIC;
  ic(2) = HQ;
  ic(3) = CV;

  out["sse"] = sse;
  out["df"] = df;
  out["ic"] = ic;
  out["U"] = U;
  out["V"] = V;
  out["D"] = D;
  out["converged"] = converged;
  return(out);
}
