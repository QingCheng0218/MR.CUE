#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <thread>
#include <mutex>
#include <ctime>
#include <cmath>
#include "ReadGeneFile.hpp"
#include "CalCorr.hpp"
#include "data_loader.hpp"
#include "GibbsGamgamEta_ptr.hpp"
#include "truncatedNormal.hpp"
#include "function.hpp"

using namespace Rcpp;
using namespace arma;
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::depends(RcppArmadillo)]]

ObjGibbs3 M3indepPXGibbsObj(arma::vec &gammah, arma::vec &Gammah, arma::vec &se1, arma::vec &se2, double &rho,
                            Options_Gibbs3* opts){
  
  // ----------------------------------------------------------------------
  // check number of input arguments
  double agm = opts -> agm;
  double bgm = opts -> bgm;
  double atau0 = opts -> atau1;
  double btau0 = opts -> btau1;
  double atau1 = opts -> atau2;
  double btau1 = opts -> btau2;
  double a = opts -> a;
  double b = opts -> b;
  uword maxIter = opts -> maxIter;
  uword thin = opts -> thin;
  uword burnin = opts -> burnin;
  // ----------------------------------------------------------------------
  // initial values
  int p = Gammah.n_elem;
  double sgga2 = 0.01; double tau02 = 0.01; double tau12 = 0.01;
  double beta0 = 0.01; double beta1 = 0.01; double w = 0.1;
  double xi2 = 0.01;
  
  
  double tau02xi2 = tau02*xi2;
  int numsave = maxIter / thin;
  vec Beta0res = ones(numsave, 1);
  vec Beta1res = ones(numsave, 1);
  vec Sgga2Res = ones(numsave, 1);
  vec Tau02Res = ones(numsave, 1);
  vec Tau12Res = ones(numsave, 1);
  vec Etarate = zeros(p, 1);
  
  vec EtaIterRate = zeros(p, 1);
  vec WRes = ones(numsave, 1);
  imat EtaAll = ones<imat>(p, numsave);
  
  vec mu = 0.01*ones(p, 1);
  vec mut = 0.01*ones(p, 1);
  
  
  ivec Eta = zeros<ivec>(p, 1);
  // ----------------------------------------------------------------------
  double inrho2 = 1./(1 - rho*rho);  // inrho2: 1/(1 - rho^2);
  double rinrho2 = rho / (1 - rho*rho);  // rinrho2: rho/(1 - rho^2);
  vec invse12 = 1. / (se1%se2);
  vec rinrhoInv12 = rinrho2*invse12; // rinrhoInv12: rho/[(1-rho^2)(se1*se2)];
  
  vec sG2 = se2%se2;
  vec sg2 = se1%se1;
  vec invsG2 = 1. / sG2;
  vec invsg2 = 1. / sg2;
  vec GinvsG2 = Gammah / sG2;
  vec ginvsg2 = gammah / sg2;
  
  int l = 0;
  for(int iter = 0; iter < (int_fast32_t)(maxIter + burnin); iter ++){
    double invsgga2 = 1. / sgga2;
    double invtau02xi2 = 1. / (tau02xi2);
    double invtau12 = 1. / tau12;
    
    // ----------------------- //
    // Parameters for Gamma
    // ----------------------- //
    vec v21t, v20t = zeros(p, 1);
    v20t = 1. / (inrho2*invsG2 + invtau02xi2);
    v21t = 1. / (inrho2*invsG2 + invtau12);
    vec muinvT0 = mu*invtau02xi2;
    vec muinvT1 = mu*invtau12;
    vec rhogmgh = rinrhoInv12%(mu - gammah);
    
    for(int j = 0; j < p; j++){
      if(Eta[j]==1){
        double mut1 = (inrho2*GinvsG2[j] + rhogmgh[j] + beta1*muinvT1[j])*v21t[j];
        mut[j] = mut1 + randn()*sqrt(v21t[j]); // mut[j] = mut1;
      }else{
        double mut0 = (inrho2*GinvsG2[j] + rhogmgh[j] + beta0*muinvT0[j])*v20t[j];
        mut[j] = mut0 + randn()*sqrt(v20t[j]); // mut[j] = mut0;
      }
    }
    // ----------------------- //
    // Parameters for gamma
    // ----------------------- //
    
    vec v21, v20 = zeros(p, 1);
    
    v20 = 1. / (inrho2*invsg2 + beta0*beta0*invtau02xi2 + invsgga2);
    v21 = 1. / (inrho2*invsg2 + beta1*beta1*invtau12 + invsgga2);
    vec mutinvT0 = mut*invtau02xi2;
    vec mutinvT1 = mut*invtau12;
    vec rhoGmGh = rinrhoInv12 % (mut - Gammah);
    for(int j = 0; j < p; j++){
      if(Eta[j]==1){
        double mu1 = (inrho2*ginvsg2[j] + rhoGmGh[j] + beta1*mutinvT1[j])*v21[j];
        mu[j] = mu1 + randn()*sqrt(v21[j]); // mu[j] = mu1;
      }else{
        double mu0 = (inrho2*ginvsg2[j] + rhoGmGh[j] + beta0*mutinvT0[j])*v20[j];
        mu[j] = mu0 + randn()*sqrt(v20[j]); // mu[j] = mu0;
      }
      
    }
    // ----------------------- //
    // Update beta0, beta1;
    // ----------------------- //
    
    vec muEta0, muEta1;
    muEta1 = mu % Eta;
    muEta0 = mu % (1 - Eta);
    // ----------------------- //
    // Update beta0, beta1
    // ----------------------- //
    double sig2b0, mub0;
    if(sum(Eta)==p){
      beta0 = 0;
    }else{
      sig2b0 = 1. / (sum(muEta0 % muEta0)*invtau02xi2);
      mub0 = (sum(muEta0 % mut)*invtau02xi2)*sig2b0;
      beta0 = mub0 + randn()*sqrt(sig2b0); // beta0 = mub0;
    }
    double sig2b1, mub1;
    if(sum(Eta)==0){
      beta1 = 0;
    }else{
      sig2b1 = 1. / (sum(muEta1 % muEta1)*invtau12);
      mub1 = (sum(muEta1 % mut)*invtau12)*sig2b1;
      beta1 = mub1 + randn()*sqrt(sig2b1);
    }
    
    
    double tagm, tbgm, tatau0, tbtau0, tatau1, tbtau1, taxi2, tbxi2;
    
    // ----------------------- //
    // Update sgga2
    // ----------------------- //
    
    tagm = agm + p / 2;
    tbgm = as_scalar(mu.t()*mu) / 2 + bgm;
    
    sgga2 =  1 / randg<double>(distr_param(tagm, 1/tbgm)); // sgga2 = tbgm;
    
    // ----------------------- //
    // Update tau02
    // ----------------------- //
    vec err02 = (1 - Eta)%(mut - beta0*mu)%(mut - beta0*mu);
    vec err12 = Eta%(mut - beta1*mu)%(mut - beta1*mu);
    
    
    tatau0 = atau0 + 0.5*sum(1 - Eta);
    tbtau0 = btau0 + 0.5*sum(err02)/xi2;
    vec tmp0 = (mut - beta0*mu)%(mut - beta0*mu);
    vec tmp1 = (mut - beta1*mu)%(mut - beta1*mu);
    
    
    if(sum(Eta)!=p && tbtau0!=0 && tatau0!=0){
      tau02 =  1 / randg<double>(distr_param(tatau0, 1/tbtau0));  // tau02 = tbtau0;
    }
    
    // ----------------------- //
    // Update tau12
    // ----------------------- //
    
    tatau1 = atau1 + 0.5*sum(Eta);
    tbtau1 = btau1 + 0.5*sum(err12);
    
    
    if(tbtau1<0 || tatau1 <0){
      cout << "error!!" << "tbtau1:" << tbtau1 << "tatau1:" << tatau1<< endl;
    }
    
    if(sum(Eta)!=0 && tbtau1!=0 && tatau1!=0){
      tau12 =  1 / randg<double>(distr_param(tatau1, 1/tbtau1)); // tau12 = tbtau1;
    }
    
    
    // ----------------------- //
    // Update xi2
    // ----------------------- //
    taxi2 = 0.5*sum(1 - Eta);
    tbxi2 = 0.5*sum(err02)/tau02;
    if(taxi2<0 || tbxi2 <0){
      cout << "error!!" << "taxi2:" << taxi2 << "tbxi2:" << tbxi2<< endl;
    }
    if(taxi2!=0 && tbxi2!=0){
      xi2 =  1 / randg<double>(distr_param(taxi2, 1/tbxi2));  
    }
    
    tau02xi2 = tau02*xi2;
    if(tau02xi2 < 1e-7){
      tau02xi2 = 1e-7;
    }
    
    // update omega
    // ------------------ //
    double wpa1, wpa2;
    wpa1 = a + sum(Eta);
    wpa2 = b + p - sum(Eta);
    w = R::rbeta(wpa1, wpa2); // w = 0.1;
    
    
    
    // ----------------------- //
    // Update Eta
    // ----------------------- //
    // double tau0 = sqrt(tau02);
    double tau1 = sqrt(tau12);
    double tau02xi = sqrt(tau02xi2);
    vec bmu0 = beta0*mu;
    vec bmu1 = beta1*mu;
    double pr1, pr0, prob;
    vec temp = zeros(p, 1);
    for(int j = 0; j < p; j++){
      pr1 = w*normal_pdf(mut[j], bmu1[j], tau1);
      pr0 = (1 - w)*normal_pdf(mut[j], bmu0[j], tau02xi);
      prob = pr1 / (pr0 + pr1);
      
      Eta[j] = R::rbinom(1, prob);
      temp[j] = prob;
    }
    
    if(iter >= (int)burnin){
      if((iter - burnin) % thin ==0){
        Etarate = (Etarate*l + temp)/(l + 1);
        EtaAll.col(l) = Eta;
        WRes[l] = w;
        Sgga2Res[l] = sgga2;
        Tau02Res[l] = tau02xi2;
        Tau12Res[l] = tau12;
        Beta0res[l] = beta0;
        Beta1res[l] = beta1;
        l += 1;
      }
    }
  }
  
  for(int i = 0; i< (int_fast32_t)EtaAll.n_rows; i++){
    EtaIterRate[i] = (double)sum(EtaAll.row(i))/(double)EtaAll.n_cols;
  }
 
  // The identification of beta1
  if(mean(Tau02Res) > mean(Tau12Res)){
    vec breplace = Beta0res;
    Beta0res = Beta1res;
    Beta1res = breplace;
    breplace = Tau02Res;
    Tau02Res = Tau12Res;
    Tau12Res = Tau02Res;
    Etarate = 1 - Etarate;
    EtaIterRate = 1 - EtaIterRate;
    WRes = 1 - WRes;
  }
  
  
  ObjGibbs3 obj;
  
  obj.Eta = Eta;
  obj.EtaIterRate = EtaIterRate;
  obj.WRes = WRes;
  obj.Etarate = Etarate;
  obj.Beta0res = Beta0res;
  obj.Beta1res = Beta1res;
  obj.Sgga2Res = Sgga2Res;
  obj.Tau02Res = Tau02Res;
  obj.Tau12Res = Tau12Res;
  
  return obj;
  
}


//[[Rcpp::export]]
List MRCUEIndep(arma::vec &gammah, arma::vec &Gammah, arma::vec &se1, arma::vec &se2, double &rho,
                SEXP opts = R_NilValue)
{
  Options_Gibbs3* lp_opt = NULL;
  uword nblocks = gammah.n_elem;
  
  if (!Rf_isNull(opts)){
    Rcpp::List opt(opts);
    lp_opt = new Options_Gibbs3(opt["agm"], opt["bgm"], opt["atau1"], opt["btau1"], 
                                opt["atau2"], opt["btau2"], opt["a"], opt["b"], 
                                    opt["maxIter"], opt["thin"], opt["burnin"]);
  }
  if (Rf_isNull(opts)){
    lp_opt = new Options_Gibbs3(nblocks);
  }
  
  ObjGibbs3 obj = M3indepPXGibbsObj(gammah, Gammah, se1, se2, rho, lp_opt);
  
  double bhat = mean(obj.Beta0res);
  double se = stddev(obj.Beta0res);
  double pvalue = 2*(R::pnorm(abs(bhat / se), 0, 1, 0, 0));
  
  List output = List::create(
    Rcpp::Named("beta.hat") = bhat,
    Rcpp::Named("beta.se") = se,
    Rcpp::Named("beta.p.value") = pvalue,
    Rcpp::Named("Eta") = Rcpp::wrap(obj.Eta),
    Rcpp::Named("WRes") = Rcpp::wrap(obj.WRes),
    Rcpp::Named("EtaIterRate") = Rcpp::wrap(obj.EtaIterRate),
    Rcpp::Named("Etarate") = Rcpp::wrap(obj.Etarate),
    Rcpp::Named("Beta1res") = Rcpp::wrap(obj.Beta0res), // change notation.
    Rcpp::Named("Beta2res") = Rcpp::wrap(obj.Beta1res), // change notation.
    Rcpp::Named("Sgga2Res") = Rcpp::wrap(obj.Sgga2Res),
    Rcpp::Named("Tau12Res") = Rcpp::wrap(obj.Tau02Res), // change notation.
    Rcpp::Named("Tau22Res") = Rcpp::wrap(obj.Tau12Res)  // change notation.
    
  );
  return output;
}


ObjGibbs3 gibbs3group_blockpar(arma::field<vec> F4gammah, arma::field<vec> F4Gammah,
                                  arma::field<vec> F4se1, arma::field<vec> F4se2,
                                  arma::field<mat> F4Rblock, double rho,  int coreNum, 
                                  Options_Gibbs3* opts){
  // ----------------------------------------------------------------------
  // check number of input arguments
  double agm = opts -> agm;
  double bgm = opts -> bgm;
  double atau0 = opts -> atau1;
  double btau0 = opts -> btau1;
  double atau1 = opts -> atau2;
  double btau1 = opts -> btau2;
  double a = opts -> a;
  double b = opts -> b;
  uword maxIter = opts -> maxIter;
  uword thin = opts -> thin;
  uword burnin = opts -> burnin;
  
  // ----------------------------------------------------------------------
  // initial values
  uword nblocks = F4Rblock.size();
  ivec NB = zeros<ivec>(nblocks, 1);
  for (int nn = 0; nn < (int)(nblocks); nn = nn+1){
    NB[nn] = F4gammah(nn, 0).size();
  }
  int p = sum(NB);
  
  double sgga2 = 0.01;
  double beta0 = 0.01; double beta1 = 0.01; double w = 0.1;
  
  double tau12 = 0.01;
  double tau02 = 0.01;
  double xi2 = 0.01;
  
  double tau02xi2 = tau02*xi2;
  
  int numsave = maxIter / thin;
  vec Beta0res = ones(numsave, 1);
  vec Beta1res = ones(numsave, 1);
  vec Sgga2Res = ones(numsave, 1);
  vec Tau02Res = ones(numsave, 1);
  vec Tau12Res = ones(numsave, 1);
  
  
  vec mu = 0.01*ones(p, 1);
  vec mut = 0.01*ones(p, 1);
  
  
  double inrho2 = 1./(1 - rho*rho);
  double rinrho2 = rho / (1 - rho*rho);
  
  
  
  ivec Eta = zeros<ivec>(nblocks, 1);
  vec Etarate = zeros(nblocks, 1);
  
  
  vec EtaIterRate = zeros(nblocks, 1);
  vec WRes = ones(numsave, 1);
  imat EtaAll = ones<imat>(nblocks, numsave);
  
  
  field<vec> F4mu(nblocks, 1), F4mut(nblocks, 1), F4rRinsGmut(nblocks, 1), F4rRinsgmu(nblocks, 1);
  field<vec> F4rGinvsG2(nblocks, 1), F4rginvsg2(nblocks, 1), F4rginvse12(nblocks, 1), F4rGinvse12(nblocks, 1);
  field<mat> F4rinsgRinsG(nblocks, 1), F4rRins(nblocks, 1), F4rRins2(nblocks, 1);
  field<vec> F4rDinsGRinsG(nblocks, 1), F4rDinsgRinsg(nblocks, 1);
  
  
  
  for (int nn = 0; nn < (int)(nblocks); nn = nn+1){
    
    F4mu(nn, 0) = 0.01*ones(NB[nn], 1);
    F4mut(nn, 0) = 0.01*ones(NB[nn], 1);
    vec se1_block = F4se1(nn, 0);
    vec se2_block = F4se2(nn, 0);
    vec sg2_block = pow(se1_block, 2);
    vec sG2_block = pow(se2_block, 2);
    vec bh1_block = F4gammah(nn, 0);
    vec bh2_block = F4Gammah(nn, 0);
    
    mat R_block =  symmatu(F4Rblock(nn, 0));
    F4rGinvsG2(nn, 0) = inrho2*bh2_block / sG2_block;
    F4rginvsg2(nn, 0) = inrho2*bh1_block / sg2_block;
    
    F4rginvse12(nn, 0) = rinrho2*bh1_block / (se1_block % se2_block);
    F4rGinvse12(nn, 0) = rinrho2*bh2_block / (se1_block % se2_block);
    F4rDinsGRinsG(nn, 0) = inrho2*diagvec(R_block) / sG2_block;
    F4rDinsgRinsg(nn, 0) = inrho2*diagvec(R_block) / sg2_block;
    F4rRins(nn, 0) = inrho2*R_block*diagmat(1 / se1_block);
    F4rRins2(nn, 0) = inrho2*R_block*diagmat(1 / se2_block);
    F4rRinsGmut(nn, 0) = inrho2*R_block*diagmat(1 / se2_block)*F4mut(nn, 0);
    F4rRinsgmu(nn, 0) = inrho2*R_block*diagmat(1 / se1_block)*F4mu(nn, 0);
  }
  
  
  int num = 0;
  for(int iter = 0; iter < (int_fast32_t)(maxIter + burnin); iter ++){
    
    double invtau02xi2 = 1. / (tau02xi2);
    double invtau12 = 1. / tau12;
    double logW = log(w / (1 - w));
    vec mu24bE0 = zeros(nblocks, 1);
    vec mu24bE1 = zeros(nblocks, 1);
    vec mumut4bE0 = zeros(nblocks, 1);
    vec mumut4bE1 = zeros(nblocks, 1);
    vec Mu2 = zeros(nblocks, 1);
    
    vec Eta0sum = zeros(nblocks, 1);
    vec Eta1sum = zeros(nblocks, 1);
    
    // ------------------------------------------------------------------
    // set parallel computation for Gamma, gamma, eta;
    paraBlock_GamgamEta parobj_GamgamEta(nblocks, F4se1, F4se2, F4mu, F4mut,
                                         F4rRinsGmut, F4rRinsgmu, F4rGinvsG2, F4rginvsg2,
                                         F4rginvse12, F4rGinvse12, F4rDinsGRinsG, F4rDinsgRinsg,
                                         F4rRins,F4rRins2, mu24bE0, mu24bE1, mumut4bE0, mumut4bE1,
                                         Mu2, Eta, Eta0sum, Eta1sum,
                                         rho, beta0, beta1, logW, sgga2, tau12, tau02xi2);
    
    const int n_thread = coreNum;
    // const int n_thread = 1;
    std::vector<std::thread> threads(n_thread);
    for(int i_thread = 0; i_thread < n_thread; i_thread++){
      threads[i_thread] = std::thread(&paraBlock_GamgamEta::update_by_thread_GamgamEta, &parobj_GamgamEta, i_thread);
    }
    
    
    for(int i = 0; i < n_thread; i++){
      threads[i].join();
    }
    
    // save the parallel result
    mu24bE0 = parobj_GamgamEta.mu24bE0;
    mu24bE1 = parobj_GamgamEta.mu24bE1;
    mumut4bE0 = parobj_GamgamEta.mumut4bE0;
    mumut4bE1 = parobj_GamgamEta.mumut4bE1;
    Mu2 = parobj_GamgamEta.Mu2;
    
    Eta0sum = parobj_GamgamEta.Eta0sum;
    Eta1sum = parobj_GamgamEta.Eta1sum;
    Eta = parobj_GamgamEta.Eta;
    sgga2 = parobj_GamgamEta.sgga2;
    F4mu = parobj_GamgamEta.F4mu;
    F4mut = parobj_GamgamEta.F4mut;
    F4rRinsGmut = parobj_GamgamEta.F4rRinsGmut;
    F4rRinsgmu = parobj_GamgamEta.F4rRinsgmu;
    
    // ------------------------------------------------------------------
    
    // ------------------------------
    // Update beta0 beta1;
    // ------------------------------
    double sig2b0, sig2b1, mub0, mub1;
    sig2b0 = 1. / (invtau02xi2 * sum(mu24bE0));
    sig2b1 = 1. / (invtau12 * sum(mu24bE1));
    mub0 = sig2b0 * (invtau02xi2 * sum(mumut4bE0));
    mub1 = sig2b1 * (invtau12 * sum(mumut4bE1));
    
    if(sum(Eta)==(int_fast32_t)nblocks){
      beta0 = 0;
    }else{
      beta0 = mub0 + randn()*sqrt(sig2b0); // beta0 = mub0;
    }
    
    if(sum(Eta)==0){
      beta1 = 0;
    }else{
      beta1 = mub1 + randn()*sqrt(sig2b1); // beta1 = mub1;
    }
    
    // ------------------------------
    // Update sgga2, tau02, tau12, xi2;
    // ------------------------------
    double tagm, tbgm, tatau0, tbtau0, tatau1, tbtau1, taxi2, tbxi2;
    tagm = agm + p / 2;
    tbgm = sum(Mu2) / 2 + bgm;
    sgga2 =  1 / randg<double>(distr_param(tagm, 1/tbgm)); // sgga2 = tbgm;
    
    
    vec tmp4tau02 = zeros(nblocks, 1);
    vec tmp4tau12 = zeros(nblocks, 1);
    
    for(int l = 0; l < (int_fast32_t)nblocks; l++){
      tmp4tau02[l] = (1 - Eta[l])*sum((F4mut(l, 0) - beta0*F4mu(l, 0))%(F4mut(l, 0) - beta0*F4mu(l, 0)));
      tmp4tau12[l] = Eta[l]*sum((F4mut(l, 0) - beta1*F4mu(l, 0))%(F4mut(l, 0) - beta1*F4mu(l, 0)));
    }
    
    tatau0 = atau0 + 0.5*sum(Eta0sum);
    tbtau0 = btau0 + 0.5*sum(tmp4tau02) / xi2;
    
    if(sum(Eta)!=(int_fast32_t)nblocks && tbtau0!=0 && tatau0!=0){
      tau02 =  1 / randg<double>(distr_param(tatau0, 1/tbtau0)); // tau02 = tbtau0;
    }
    
    tatau1 = atau1 + 0.5*sum(Eta1sum);
    tbtau1 = btau1 + 0.5*sum(tmp4tau12);
    
    if(sum(Eta)!=0 && tbtau1!=0 && tatau1!=0){
      tau12 =  1 / randg<double>(distr_param(tatau1, 1/tbtau1)); // tau12 = tbtau1;
    }
    
    // ------------------------------
    // Update xi2;
    // ------------------------------
    taxi2 = 0.5*sum(Eta0sum);
    tbxi2 = 0.5*sum(tmp4tau02) / tau02;
    if(taxi2!=0 && tbxi2!=0){
      xi2 =  1 / randg<double>(distr_param(taxi2, 1/tbxi2)); // xi2 = tbxi2;
    }
    
    
    tau02xi2 = tau02*xi2;
    if(tau02xi2 < 1e-7){
      tau02xi2 = 1e-7;
    }
    // ------------------------------
    // update omega
    // ------------------------------
    
    double wpa1, wpa2;
    wpa1 = a + sum(Eta);
    wpa2 = b + Eta.n_elem - sum(Eta);
    w = R::rbeta(wpa1, wpa2); // w = 0.1;
    
    if(iter >= (int)burnin){
      if((iter - burnin) % thin ==0){
        // Etarate = (Etarate*num + temp)/(num + 1);
        EtaAll.col(num) = Eta;
        WRes[num] =w;
        Sgga2Res[num] = sgga2;
        Tau02Res[num] = tau02xi2;
        Tau12Res[num] = tau12;
        Beta0res[num] = beta0;
        Beta1res[num] = beta1;
        num += 1;
      }
    }
  }
  
  for(int i = 0; i < (int_fast32_t)EtaAll.n_rows; i++){
    EtaIterRate[i] = (double)sum(EtaAll.row(i))/(double)EtaAll.n_cols;
  }
  

  // The identification of beta1
  if(mean(Tau02Res) > mean(Tau12Res)){
    vec breplace = Beta0res;
    Beta0res = Beta1res;
    Beta1res = breplace;
    breplace = Tau02Res;
    Tau02Res = Tau12Res;
    Tau12Res = Tau02Res;
    Etarate = 1 - Etarate;
    EtaIterRate = 1 - EtaIterRate;
    WRes = 1 - WRes;
  }

  
  ObjGibbs3 obj;
  obj.Eta = Eta;
  obj.WRes = WRes;
  obj.Etarate = Etarate;
  obj.EtaIterRate = EtaIterRate;
  obj.Beta0res = Beta0res;
  obj.Beta1res = Beta1res;
  obj.Sgga2Res = Sgga2Res;
  obj.Tau02Res = Tau02Res;
  obj.Tau12Res = Tau12Res;
  
  return obj;
  
  
  
}



// [[Rcpp::export]]
Rcpp::List MRCUESim(arma::vec &gammah, arma::vec &Gammah, arma::vec &se1, arma::vec &se2,
                    double rho, arma::mat R, arma::umat block_inf, int coreNum, SEXP opts = R_NilValue){
  
  
  
  
  uword nblocks = block_inf.n_rows;
  field<mat> F4Rblock(nblocks, 1);
  field<vec> F4gammah(nblocks, 1);
  field<vec> F4Gammah(nblocks, 1);
  field<vec> F4se1(nblocks, 1);
  field<vec> F4se2(nblocks, 1);
  
  for(int i=0; i< (int)(nblocks); i++){
    F4Rblock(i, 0) = R.submat(block_inf(i, 0), block_inf(i, 0), block_inf(i, 1), block_inf(i, 1));
    F4gammah(i, 0) = gammah.subvec(block_inf(i, 0), block_inf(i, 1));
    F4Gammah(i, 0) = Gammah.subvec(block_inf(i, 0), block_inf(i, 1));
    F4se1(i, 0) = se1.subvec(block_inf(i, 0), block_inf(i, 1));
    F4se2(i, 0) = se2.subvec(block_inf(i, 0), block_inf(i, 1));
  }
  
  
  Options_Gibbs3* lp_opt = NULL;
  if (!Rf_isNull(opts)){
    Rcpp::List opt(opts);
    lp_opt = new Options_Gibbs3(opt["agm"], opt["bgm"], opt["atau1"], opt["btau1"],
                                opt["atau2"], opt["btau2"], opt["a"], opt["b"],
                                    opt["maxIter"], opt["thin"], opt["burnin"]);
  }
  if (Rf_isNull(opts)){
    lp_opt = new Options_Gibbs3(nblocks);
  }
  
  ObjGibbs3 obj = gibbs3group_blockpar(F4gammah, F4Gammah, F4se1, F4se2, F4Rblock, rho, coreNum, 
                                       lp_opt);
  double bhat = mean(obj.Beta0res);
  double se = stddev(obj.Beta0res);
  double pvalue = 2*(R::pnorm(abs(bhat / se), 0, 1, 0, 0));
  
  List output = List::create(
    Rcpp::Named("beta.hat") = bhat,
    Rcpp::Named("beta.se") = se,
    Rcpp::Named("beta.p.value") = pvalue,
    Rcpp::Named("Eta") = Rcpp::wrap(obj.Eta),
    Rcpp::Named("WRes") = Rcpp::wrap(obj.WRes),
    Rcpp::Named("EtaIterRate") = Rcpp::wrap(obj.EtaIterRate),
    Rcpp::Named("Etarate") = Rcpp::wrap(obj.Etarate),
    Rcpp::Named("Beta1res") = Rcpp::wrap(obj.Beta0res),
    Rcpp::Named("Beta2res") = Rcpp::wrap(obj.Beta1res),
    Rcpp::Named("Sgga2Res") = Rcpp::wrap(obj.Sgga2Res),
    Rcpp::Named("Tau12Res") = Rcpp::wrap(obj.Tau02Res),
    Rcpp::Named("Tau22Res") = Rcpp::wrap(obj.Tau12Res)
    
  );
  return output;
}

// [[Rcpp::export]]
Rcpp::List MRCUE(arma::field<vec> F4gammah, arma::field<vec> F4Gammah,
                     arma::field<vec> F4se1, arma::field<vec> F4se2,
                     arma::field<mat> F4Rblock, double rho, int coreNum, SEXP opts = R_NilValue){  
  
  uword nblocks = F4Rblock.size();
  Options_Gibbs3* lp_opt = NULL;
  if (!Rf_isNull(opts)){
    Rcpp::List opt(opts);
    lp_opt = new Options_Gibbs3(opt["agm"], opt["bgm"], opt["atau1"], opt["btau1"],
                                opt["atau2"], opt["btau2"], opt["a"], opt["b"],
                                    opt["maxIter"], opt["thin"], opt["burnin"]);
  }
  if (Rf_isNull(opts)){
    lp_opt = new Options_Gibbs3(nblocks);
  }
  
  ObjGibbs3 obj = gibbs3group_blockpar(F4gammah, F4Gammah, F4se1, F4se2, F4Rblock, rho, coreNum, 
                                          lp_opt);
  
  double bhat = mean(obj.Beta0res);
  double se = stddev(obj.Beta0res);
  double pvalue = 2*(R::pnorm(abs(bhat / se), 0, 1, 0, 0));
  
  List output = List::create(
    Rcpp::Named("beta.hat") = bhat,
    Rcpp::Named("beta.se") = se,
    Rcpp::Named("beta.p.value") = pvalue,
    Rcpp::Named("nblocks") = Rcpp::wrap(nblocks),
    // Rcpp::Named("Indpid") = Rcpp::wrap(Indpid),
    Rcpp::Named("Eta") = Rcpp::wrap(obj.Eta),
    Rcpp::Named("WRes") = Rcpp::wrap(obj.WRes),
    Rcpp::Named("EtaIterRate") = Rcpp::wrap(obj.EtaIterRate),
    Rcpp::Named("Etarate") = Rcpp::wrap(obj.Etarate),
    Rcpp::Named("Beta1res") = Rcpp::wrap(obj.Beta0res),
    Rcpp::Named("Beta2res") = Rcpp::wrap(obj.Beta1res),
    Rcpp::Named("Sgga2Res") = Rcpp::wrap(obj.Sgga2Res),
    Rcpp::Named("Tau12Res") = Rcpp::wrap(obj.Tau02Res),
    Rcpp::Named("Tau22Res") = Rcpp::wrap(obj.Tau12Res)
    
  );
  return output;
}