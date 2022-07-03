#ifndef truncatedNormal_hpp
#define truncatedNormal_hpp

#include "RcppArmadillo.h"
#include <iostream>
#include <sstream>
#include <wishart.h>
#include<boost/math/distributions.hpp>
using namespace Rcpp;
using namespace arma;
using namespace std;

arma::vec truncEstfun(arma::vec &a, arma::vec &b, arma::vec &x1, arma::vec &x2,
                      int maxIter, int burnin, int thin);

#endif 
