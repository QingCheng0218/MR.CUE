#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <thread>
#include <mutex>
#include <random>
#include "GibbsGamgamEta_ptr.hpp"


using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]


void paraBlock_GamgamEta::loop_by_block_gibbs_GamgamEta(int i){
  // ---------------------------------------------------------------- //
  // load the parameters in each block
  int eta = Eta[i];
  vec rRinsgmul = F4rRinsgmu(i, 0);
  vec rRinsGmutl = F4rRinsGmut(i, 0);
  vec rDinsGRinsG = F4rDinsGRinsG(i, 0);
  vec rDinsgRinsg = F4rDinsgRinsg(i, 0);
  mat rRinsl = F4rRins(i, 0);
  mat rRins2l = F4rRins2(i, 0);

  vec mutl = F4mut(i, 0);
  vec mul = F4mu(i, 0);
  vec se1l = F4se1(i, 0);
  vec se2l = F4se2(i, 0);

 
  vec rGinvsG2l = F4rGinvsG2(i, 0);
  vec rginvsg2l = F4rginvsg2(i, 0);
  vec rginvse12l = F4rginvse12(i, 0);
  vec rGinvse12l = F4rGinvse12(i, 0);

  // ---------------------------------------------------------------- //
  vec diagrRinsl = diagvec(rRinsl);
  vec diagrRins2l = diagvec(rRins2l);
  double invsgga2 = 1. / sgga2;
  double invtau02xi2 = 1. / (tau02xi2);
  double invtau12 = 1. / tau12;
  double b02 = beta0 * beta0;
  double b12 = beta1 * beta1;

  vec v20tl =  1. / (rDinsGRinsG + invtau02xi2);
  vec v21tl =  1. / (rDinsGRinsG + invtau12);

  vec v20l =  1. / (rDinsgRinsg + b02*invtau02xi2 + invsgga2);
  vec v21l =  1. / (rDinsgRinsg + b12*invtau12 + invsgga2);
  int pl = rRinsgmul.n_elem;
  // ------------------------
  // Update Gamma
  // ------------------------
  if(eta==1){
    vec tmp1;
    double rRinsGmutjj1, mut1;
    for(int j =0; j < pl; j++){
      tmp1 = rRinsGmutl - rRins2l.col(j)*mutl[j];
      rRinsGmutjj1 = rRinsGmutl[j] - diagrRins2l[j]*mutl[j];
      mut1 = (rGinvsG2l[j] - rRinsGmutjj1/se2l[j] - rginvse12l[j] + rho*rRinsgmul[j] / se2l[j] + beta1 * mul[j] * invtau12)*v21tl[j];

      mutl[j] = mut1 + randn()*sqrt(v21tl[j]); // mutl[j] = mut1;
      rRinsGmutl = tmp1 + rRins2l.col(j)*mutl[j];

    }

  }else{
    vec tmp1;
    double rRinsGmutjj0, mut0;
    for(int j =0; j < pl; j++){
      tmp1 = rRinsGmutl - rRins2l.col(j)*mutl[j];
      rRinsGmutjj0 = rRinsGmutl[j] - diagrRins2l[j]*mutl[j];

      mut0 = (rGinvsG2l[j] - rRinsGmutjj0/se2l[j] - rginvse12l[j] + rho*rRinsgmul[j] / se2l[j] + beta0 * mul[j] * invtau02xi2)*v20tl[j];

      mutl[j] = mut0 + randn()*sqrt(v20tl[j]); // mutl[j] = mut0;
      rRinsGmutl = tmp1 + rRins2l.col(j)*mutl[j];

    }
  }



  // ------------------------
  // Update gamma
  // ------------------------

  if(eta==1){
    vec tmp2;
    double rRinsgmujj1, mu1;
    for(int j =0; j < pl; j++){
      tmp2 = rRinsgmul - rRinsl.col(j)*mul[j];
      rRinsgmujj1 = rRinsgmul[j] - diagrRinsl[j]*mul[j];
      mu1 = (rginvsg2l[j] - rRinsgmujj1/se1l[j] - rGinvse12l[j] + rho*rRinsGmutl[j] / se1l[j] + beta1 * mutl[j] * invtau12)*v21l[j];

      mul[j] = mu1 + randn()*sqrt(v21l[j]); // mul[j] = mu1;
      rRinsgmul = tmp2 + rRinsl.col(j)*mul[j];
    }
  }else{
    vec tmp2;
    double rRinsgmujj0, mu0;
    for(int j =0; j < pl; j++){
      tmp2 = rRinsgmul - rRinsl.col(j)*mul[j];
      rRinsgmujj0 = rRinsgmul[j] - diagrRinsl[j]*mul[j];
      mu0 = (rginvsg2l[j] - rRinsgmujj0/se1l[j] - rGinvse12l[j] + rho*rRinsGmutl[j] / se1l[j] + beta0 * mutl[j] * invtau02xi2)*v20l[j];

      mul[j] = mu0 + randn()*sqrt(v20l[j]); // mul[j] = mu0;
      rRinsgmul = tmp2 + rRinsl.col(j)*mul[j];
    }
  }

  F4mu(i, 0) = mul;
  F4mut(i, 0) = mutl;
  F4rRinsgmu(i, 0) = rRinsgmul;
  F4rRinsGmut(i, 0) = rRinsGmutl;

  // ------------------------------
  // Update Eta
  // ------------------------------
  double prob0, prob;
  vec bmu0l = beta0*F4mu(i, 0);
  vec bmu1l = beta1*F4mu(i, 0);
  vec err1 = mutl - bmu1l;
  vec err0 = mutl - bmu0l;

  prob0 = logW - 0.5*invtau12*sum(err1%err1) - 0.5*pl*log(tau12) + 0.5*invtau02xi2*sum(err0%err0) + 0.5*pl*log(tau02xi2);
  prob = 1. / (1 + exp(-prob0));
  // cout<<"i:" << i << "prob" << prob << endl;
  
  eta = R::rbinom(1, prob);
  Eta[i] = eta;


  // -----------------------------------------------------------------------
  // for beta0 and beta1 iteration
  mu24bE0[i] = (1 - eta)*sum(F4mu(i, 0)%F4mu(i, 0));
  mu24bE1[i] = eta*sum(F4mu(i, 0)%F4mu(i, 0));
  mumut4bE0[i] = (1 - eta)*sum(F4mut(i, 0)%F4mu(i, 0));
  mumut4bE1[i] = eta*sum(F4mut(i, 0)%F4mu(i, 0));
  Mu2[i] = sum(F4mu(i, 0)%F4mu(i, 0));
  Eta0sum[i] = pl*(1 - eta);
  Eta1sum[i] = pl*Eta[i];


  // -----------------------------------------------------
  se1l.reset();
  se2l.reset();
  rRinsl.reset();
  rRins2l.reset();

  rRinsgmul.reset();
  rRinsGmutl.reset();
  rDinsGRinsG.reset();
  rDinsgRinsg.reset();
 

}

std::mutex _mtx0;
int paraBlock_GamgamEta::next_GamgamEta(){
  std::lock_guard<std::mutex> lockGuard(_mtx0);
  if(current_idx >= (int)nblocks){
    return -1;
  }
  current_idx++;
  return current_idx-1;
}

void paraBlock_GamgamEta::update_by_thread_GamgamEta(int thread_id){
  while(true){
    int idx = next_GamgamEta();
    if(idx == -1){
      break;
    }
    loop_by_block_gibbs_GamgamEta(idx);
  }
}
