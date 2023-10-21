#include "RcppArmadillo.h"
#include <iostream>
#include <sstream>
#include <wishart.h>
#include <boost/math/distributions.hpp>
using namespace Rcpp;
using namespace arma;
using namespace std;





// [[Rcpp::depends(RcppArmadillo, BH, RcppDist)]]

//[[Rcpp::export]]
double phi(double x)
{
  // constants
  double a1 =  0.254829592;
  double a2 = -0.284496736;
  double a3 =  1.421413741;
  double a4 = -1.453152027;
  double a5 =  1.061405429;
  double p  =  0.3275911;
  
  // Save the sign of x
  int sign = 1;
  if (x < 0)
    sign = -1;
  x = fabs(x)/sqrt(2.0);
  
  // A&S formula 7.1.26
  double t = 1.0/(1.0 + p*x);
  double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
  
  return 0.5*(1.0 + sign*y);
}

//[[Rcpp::export]]
double cdfNormal(double x, double mean, double sd)
{
  boost::math::normal_distribution<>myNormal (mean, sd);
  return cdf(myNormal, x);
}

//[[Rcpp::export]]
arma::vec MulticdfNormal(arma::vec &x)
{
  int p = x.n_elem;
  arma::vec res = zeros(p, 1);
  for(int i = 0; i < p; i++){
    res[i] = cdfNormal(x[i], 0, 1);
  }
  return res;
}

//[[Rcpp::export]]
arma::vec multiphi(arma::vec x){
  int p = x.n_elem;
  arma::vec res = zeros(p, 1);
  for(int i = 0; i < p; i++){
    res[i] = phi(x[i]);
  }
  return res;
}

//[[Rcpp::export]]
double inverseNormal(double prob, double mean, double sd){
  boost::math::normal_distribution<>myNormal (mean, sd);
  return quantile(myNormal, prob);
}

//[[Rcpp::export]]
arma::vec MultiinverseNormal(arma::vec &prob){
  int p = prob.n_elem;
  arma::vec res = zeros(p, 1);
  for(int i = 0; i < p; i++){
    res[i] = inverseNormal(prob[i], 0, 1);
  }
  return res;
}

//[[Rcpp::export]]
arma::mat MatSum(arma::vec &y1, arma::vec &y2){
  arma::mat RES = zeros(2, 2);
  int p = y1.n_elem;
  for(int i = 0; i < p; i++){
    RES(0, 0) = RES(0, 0) + y1[i]*y1[i];
    RES(0, 1) = RES(0, 1) + y1[i]*y2[i];
    RES(1, 1) = RES(1, 1) + y2[i]*y2[i];
  }
  RES(1, 0) = RES(0, 1);
  return RES;
  
}


//[[Rcpp::export]]
arma::vec truncEstfun(arma::vec &a, arma::vec &b, arma::vec &x1, arma::vec &x2,
                      int maxIter, int burnin, int thin){
  double sig11, sig22, sig12, tau1, w1,  w2, phi2, phi3;
  vec phi1, phi21, y1, y2, tau2, phi22, phi23;
  uword p = 2;
  int n = x1.n_elem;
  double a1 = a[0], a2 = a[1];
  double b1 = b[0], b2 = b[1];

  arma::mat sig = diagmat(ones(p, 1));
  vec mu = zeros(p, 1);


  int numsave = maxIter / thin;
  vec sigres = ones(numsave, 1);
  int l = 0;

  for(int iter = 0; iter < (int)(maxIter + burnin); iter ++){
    sig11 = sig(0, 0);
    sig22 = sig(1, 1);
    sig12 = sig(0, 1);


    tau1 = mu(0);
    w1 = sqrt(sig11);
    
    vec tmpphi1 = (x1 - tau1)/w1;

    phi1 = MulticdfNormal(tmpphi1);
    phi2 = cdfNormal((a1 - tau1)/w1, 0, 1);
    phi3 = cdfNormal((b1 - tau1)/w1, 0, 1);
    
    vec tmpy1 = (phi1 - phi2)/(phi3 - phi2);

    y1 = tau1 + w1*MultiinverseNormal(tmpy1);
    
    tau2 = mu[1] + sig12*(1./sig11)*(y1 - mu[0]);
    w2 = sqrt(sig22 - sig12*sig12 / sig11);
    
    vec tmpphi21 = (x2 - tau2)/w2;
    vec tmpphi22 = (a2 - tau2)/w2;
    vec tmpphi23 = (b2 - tau2)/w2;
 
    phi21 = MulticdfNormal(tmpphi21);
    phi22 = MulticdfNormal(tmpphi22);
    phi23 = MulticdfNormal(tmpphi23);
    vec tmpy2 = (phi21 - phi22)/(phi23 - phi22);

    y2 = tau2 + w2*MultiinverseNormal(tmpy2);
    
    mat Sy = MatSum(y1, y2);
  
    sig = riwish(n, Sy);
    if(iter >= (int)burnin){
      if((iter - burnin) % thin ==0){
        sigres[l] = sig(0, 1)/sqrt(sig(0, 0)*sig(1, 1));
        l += 1;
      }
    }


  }

  return sigres;
}

