#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include "CalCorr.hpp"

using namespace Rcpp;
using namespace arma;
using namespace std;




mat cal_blockcor(arma::umat& X){
  uword size1 = X.n_rows;
  uword size2 = X.n_cols;
  vec meanX(size2);
  vec sqrtsum(size2);
  mat Xnorm = conv_to<mat>::from(X);
  for (int i = 0; i < (int)(size2); i++) { //calculate the mean of the vector and sqrt sum
    meanX[i] = sum(Xnorm.col(i))*1.0 / size1;
    Xnorm.col(i) -= meanX[i];
    vec v_i = Xnorm.col(i);
    mat pd = v_i.t() * v_i;
    sqrtsum[i] = sqrt(pd.at(0));
    if (sqrtsum[i] > 0){
      Xnorm.col(i) /= sqrtsum[i];
    }
  }
  arma::mat corr(size2, size2);
  arma::mat eyeI(size2, size2);
  eyeI.eye();
  mat cor_ij(1, 1);
  corr.eye();
  for (int i = 0; i < (int)(size2); i++){
    for (int j = i + 1; j < (int)(size2); j++) {
      cor_ij = ((Xnorm.col(i)).t() * (Xnorm.col(j)));
      double value = cor_ij.at(0);
      corr(i, j) = value;
      corr(j, i) = value;
    }
  }

  return corr;
  
}


// [[Rcpp::export]]
mat Cal_block_SimR(umat block_inf, arma::umat &X, double lam){
  
  
  int p = X.n_cols;
  int nblocks = block_inf.n_rows;
  mat R = zeros(p, p);
  for(int j = 0; j < (int)(nblocks); j++){
    umat subX = X.cols(block_inf(j, 0), block_inf(j, 1));
    
    arma::mat corr = cal_blockcor(subX);
    int p1 = corr.n_rows;
    arma::mat eyeI(p1, p1);
    eyeI.eye();
    corr *= lam;
    corr += (1 - lam)*eyeI;
    R.submat(block_inf(j, 0), block_inf(j, 0), block_inf(j, 1), block_inf(j, 1)) = corr;
    
  }
  
  
  return R;
}


// [[Rcpp::depends(RcppArmadillo)]]
void fastLm(const arma::vec & y, const arma::mat & X, double &coefb, double &stdb) {
  
  int n = X.n_rows, k = X.n_cols;
  
  arma::colvec coef = arma::solve(X, y); 
  arma::colvec resid = y - X*coef; 
  
  double sig2 = arma::as_scalar(arma::trans(resid)*resid/(n-k));
  arma::colvec stderrest = 
    arma::sqrt(sig2 * arma::diagvec( arma::inv(arma::trans(X)*X)) );
  
  coefb = coef[1];
  stdb = stderrest[1];
}

// [[Rcpp::export]]
List fastSigLm(const arma::vec & y, const arma::mat & X) {
  
  // int n = X.n_rows, k = X.n_cols;
  int p = X.n_cols;int n = X.n_rows;
  arma::mat xx = zeros(p, 2);
  double coefb = 0;
  double stdb = 0;
  vec coef = zeros(p, 1);
  vec std = zeros(p, 1);
  
  for( int j = 0; j < p; j = j + 1 )
  {
    xx = join_rows(ones(n, 1), X.col(j));
    fastLm(y, xx, coefb, stdb);
    coef[j] = coefb;
    std[j] = stdb;
  }
  
  return List::create(Named("coef") = coef,
                      Named("std") = std);
}


//[[Rcpp::export]]
double normal_pdf(double x, double m, double s)
{
  static const double inv_sqrt_2pi = 0.3989422804014327;
  double a = (x - m) / s;
  
  return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
}

