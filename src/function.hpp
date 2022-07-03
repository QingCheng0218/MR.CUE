#ifndef function_hpp
#define function_hpp

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>

using namespace Rcpp;
using namespace arma;
using namespace std;

Rcpp::List fastSigLm(const arma::vec & y, const arma::mat & X);
double normal_pdf(double x, double m, double s);
arma::mat Cal_block_SimR(umat block_inf, arma::umat &X, double lam);
#endif /* function_hpp */