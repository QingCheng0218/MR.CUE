#ifndef GibbsGamgamEta_ptr_hpp
#define GibbsGamgamEta_ptr_hpp

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <thread>
#include <mutex>


using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

class Options_Gibbs3{
public:
  // Constructor definition
  // The complier deciedes which constructor to be called depending on 
  // the number of argument present with the object
  Options_Gibbs3(uword nblocks){
    this -> agm = 0;
    this -> bgm = 0;
    this -> atau1 = 0;
    this -> btau1 = 0;
    this -> atau2 = 0;
    this -> btau2 = 0;
    this -> a = 2;
    this -> b = nblocks;
    this -> CluDes = "PropMajor";
    this -> maxIter = 4000;
    this -> thin = 10;
    this -> burnin = 1000;
  }
  
  Options_Gibbs3(double agm, double bgm, double atau1, double btau1,  double atau2, double btau2, 
                 double a, double b, string CluDes, uword maxIter, uword thin, uword burnin){
    
    this -> agm = agm;
    this -> bgm = bgm;
    this -> atau1 = atau1;
    this -> btau1 = btau1;
    this -> atau2 = atau2;
    this -> btau2 = btau2;
    this -> a = a;
    this -> b = b;
    this -> CluDes = CluDes;
    this -> maxIter = maxIter;
    this -> thin = thin;
    this -> burnin = burnin;
    
  }
  double agm;
  double bgm;
  double atau1;
  double btau1;
  double atau2;
  double btau2;
  double a;
  double b;
  string CluDes;
  uword maxIter;
  uword thin;
  uword burnin;
  uword nblocks;
  
};
struct ObjGibbs3{
  ivec Eta;
  vec Etarate;
  vec EtaIterRate;
  vec Beta0res;
  vec Beta1res;
  vec Sgga2Res;
  vec Tau02Res;
  vec Tau12Res;
  vec WRes;
};

class paraBlock_GamgamEta{
public:
  int current_idx=0;
  int n_thread = 1;
  umat block_inf;
  
  uword nblocks;
  
  int conspar;
  ivec Eta;
  vec mu24bE0, mu24bE1, mumut4bE0, mumut4bE1, Mu2, Eta0sum, Eta1sum;
  double rho, beta0, beta1, logW, sgga2, tau12, tau02xi2;
  
  field<vec> F4se1, F4se2, F4mu, F4mut;
  field<vec> F4rRinsGmut, F4rRinsgmu, F4rGinvsG2, F4rginvsg2, F4rginvse12, F4rGinvse12;
  field<vec> F4rDinsGRinsG, F4rDinsgRinsg;
  field<mat> F4rRins, F4rRins2;
 
  
  
  paraBlock_GamgamEta(uword &nblocks, field<vec> &F4se1, field<vec> &F4se2, field<vec> &F4mu, field<vec> &F4mut,
                      field<vec> &F4rRinsGmut, field<vec> &F4rRinsgmu, field<vec> &F4rGinvsG2, field<vec> &F4rginvsg2,
                      field<vec> &F4rginvse12, field<vec> &F4rGinvse12, field<vec> &F4rDinsGRinsG, field<vec> &F4rDinsgRinsg,
                      field<mat> &F4rRins, field<mat> &F4rRins2,
                      arma::vec &mu24bE0, arma::vec &mu24bE1, arma::vec &mumut4bE0, arma::vec &mumut4bE1,
                      arma::vec &Mu2, arma::ivec Eta, arma::vec &Eta0sum, arma::vec &Eta1sum,
                      double &rho, double &beta0, double &beta1, double &logW, double &sgga2,
                      double &tau12, double &tau02xi2){
    
    this -> nblocks = nblocks;
    this -> F4se1 = F4se1;
    this -> F4se2 = F4se2;
    this -> F4mu = F4mu;
    this -> F4mut = F4mut;
    this -> F4rRinsGmut = F4rRinsGmut;
    this -> F4rRinsgmu = F4rRinsgmu;
    this -> F4rGinvsG2 = F4rGinvsG2;
    this -> F4rginvsg2 = F4rginvsg2;
    this -> F4rginvse12 = F4rginvse12;
    this -> F4rGinvse12 = F4rGinvse12;
    this -> F4rDinsGRinsG = F4rDinsGRinsG;
    this -> F4rDinsgRinsg = F4rDinsgRinsg;
    this -> F4rRins = F4rRins;
    this -> F4rRins2 = F4rRins2;
    this -> mu24bE0 = mu24bE0;
    this -> mu24bE1 = mu24bE1;
    this -> mumut4bE0 = mumut4bE0;
    this -> mumut4bE1 = mumut4bE1;
    this -> Mu2 = Mu2;
    this -> Eta = Eta;
    this -> Eta0sum = Eta0sum;
    this -> Eta1sum = Eta1sum;
    this -> rho = rho;
    this -> beta0 = beta0;
    this -> beta1 = beta1;
    this -> logW = logW;
    this -> sgga2 = sgga2;
    this -> tau12 = tau12;
    this -> tau02xi2 = tau02xi2;
  }
  
  
  int  next_GamgamEta();
  void loop_by_block_gibbs_GamgamEta(int i);
  void update_by_thread_GamgamEta(int thread_id);
  
};

#endif /* GibbsGamgamEta_ptr_hpp */
