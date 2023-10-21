#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "time.h"
#include <iostream>
#include "ReadGeneFile.hpp"

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

#define MAX_LEN 20
//-----------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List ReadSNPinfo(std::string stringname, IntegerVector A1, IntegerVector A2, CharacterVector rsname,
                       IntegerVector chr, IntegerVector bp, NumericVector morgan, int N)
{
  FILE *stream;
  int ch, b;
  double mor;
  char s[MAX_LEN + 1], efa, nefa;
  
  stream = fopen(stringname.c_str(), "r");
  clock_t t1 = clock();
  //int i = 0;
  /* Put in various data. */
  for (int i = 0; i < N; i++){
    if (i % 100000 == 0 && i != 0){
      cout << i << "-th SNP" << ",";
      cout << "Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
    }
    
    fscanf(stream, "%i %s %lf %i %c %c", &ch, &s[0], &mor, &b, &efa, &nefa);
    
    chr(i) = ch;
    rsname(i) = s;
    morgan(i) = mor;
    bp(i) = b;
    A1(i) = (int)efa;
    A2(i) = (int)nefa;
    
  }
  
  List output = List::create(Rcpp::Named("chr") = chr,
                             Rcpp::Named("rsname") = rsname,
                             Rcpp::Named("morgan") = morgan,
                             Rcpp::Named("bp") = bp,
                             Rcpp::Named("A1") = A1,
                             Rcpp::Named("A2") = A2);
  return output;
}


// [[Rcpp::export]]
void Read_summarystat(std::string stringname, IntegerVector SA1, IntegerVector SA2, CharacterVector rsname,
                      NumericVector betah, NumericVector s2, NumericVector pvalue, IntegerVector chr, IntegerVector bp, int N){
  
  FILE *stream;
  stream = fopen(stringname.c_str(), "r");
  
  char s[MAX_LEN + 1], efa, nefa;
  double bhat, se2, pv;
  int b, ch;
  
  clock_t t1 = clock();
  
  for (int i = 0; i < (N+8); i++){
    if (i % 200000 == 0 && i != 0){
      cout << i << "-th SNP" << ",";
      cout << "Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
    }
    
    fscanf(stream, "%s %i  %i %c %c %lf %lf %lf", &s[0], &ch, &b, &efa, &nefa, &bhat, &se2, &pv);
    // fscanf(stream, "%s %c %c %lf %lf %lf", &s[0], &efa, &nefa, &bhat, &se2, &samsize);
    if(i > 7){
      //ignore the colnames
      rsname(i - 8) = s;
      //cout << s << efa << nefa << bhat << se2 << samsize << ";" << endl;
      SA1(i - 8) = (int)efa;
      SA2(i - 8) = (int)nefa;
      betah(i - 8) = bhat;
      s2(i - 8) = se2;
      chr(i - 8) = ch;
      bp(i - 8) = b;
      pvalue(i - 8) = pv;
    }
    
    
  }
  // List output = List::create(Rcpp::Named("SA1") = SA1,
  //                            Rcpp::Named("SA2") = SA2,
  //                            Rcpp::Named("betah") = betah,
  //                            Rcpp::Named("s2") = s2,
  //                            Rcpp::Named("chr") = chr,
  //                            Rcpp::Named("bp") = bp,
  //                            Rcpp::Named("pvalue") = pvalue,
  //                            Rcpp::Named("rsname") = rsname);
  // return output;
}
// [[Rcpp::export]]
CharacterVector select(CharacterVector vec_, NumericVector idx_){
  arma::vec vv  = Rcpp::as<arma::vec>(vec_);
  arma::uvec ii = Rcpp::as<arma::uvec>(idx_);
  arma::vec aa  = vv.elem(ii);
  Rcpp::CharacterVector bb(as<Rcpp::NumericMatrix>(wrap(aa)));
  return Rcpp::wrap( bb );
}


// [[Rcpp::export]]
Rcpp::List matchsnp(std::string stringname1, std::string stringname2,
                    std::string stringname3, bool matchExp){
  // summary data 1: stringname1; summary data 2:stringname2; Panel SNP info:stringname3;
  
  //----------------------------------------------------------------------------------------------//
  string bimfile = stringname3;
  bimfile += ".bim";
  
  int N = getLineNum(bimfile);
  cout << "Number of SNPs (panel):" << N << endl;
  
  
  IntegerVector A31(N), A32(N);
  CharacterVector rsname3(N);
  IntegerVector chr3(N), bp3(N);
  NumericVector morgan(N);
  
  //read from SNPinfo file (pass as pointer)
  cout << endl;
  cout << "Start loading SNP info:" << endl;
  clock_t t1 = clock();
  ReadSNPinfo(bimfile, A31, A32, rsname3, chr3, bp3, morgan, N);
  cout << "Finish loading SNP info fro panel data in " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec." << endl;
  
  //----------------------------------------------------------------------------------------------//
  // summary stat SNP info for trait 1, betah, s2, sample_size
  int Ns1 = getLineNum(stringname1) - 1;// ignore the colnames line
  
  cout << "Number of SNPs (exposure):" << Ns1 << endl;
  
  IntegerVector SA11(Ns1), SA12(Ns1), chr1(Ns1), bp1(Ns1);
  CharacterVector rsname1(Ns1);
  NumericVector betah1(Ns1), s12(Ns1), PV1(Ns1), sample_size(Ns1);
  
  cout << endl;
  printf("Start loading exposure summary stat: \n");
  t1 = clock();
  Read_summarystat(stringname1, SA11, SA12, rsname1, betah1, s12, PV1, chr1, bp1, Ns1);
  cout << "Finish loading exposure summary stat in " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec." << endl;
  //----------------------------------------------------------------------------------------------//
  // summary stat SNP info for trait 2, betah, s2, sample_size
  int Ns2 = getLineNum(stringname2) - 1;// ignore the colnames line
  cout << "Number of SNPs (outcome):" << Ns2 << endl;
  
  IntegerVector SA21(Ns2), SA22(Ns2), chr2(Ns2), bp2(Ns2);
  CharacterVector rsname2(Ns2);
  NumericVector betah2(Ns2), s22(Ns2), PV2(Ns2);
  
  cout << endl;
  printf("Start loading  outcome summary stat 2: \n");
  t1 = clock();
  Read_summarystat(stringname2, SA21, SA22, rsname2, betah2, s22, PV2, chr2, bp2, Ns2);
  
  cout << "Finish loading outcome summary stat 2 in " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec." << endl;
  //----------------------------------------------------------------------------------------------//
  //mathcing panel SNP with summary stat and correct direction of minor allele
  CharacterVector rs_tmp = intersect(rsname1, rsname2);
  CharacterVector rs_inter = intersect(rs_tmp, rsname3);
  //----------------------------------------------------------------------------------------------//
  //Panel Data
  IntegerVector idxin3 = match(rsname3, rs_inter);
  
  CharacterVector rsname_4use = rsname3[Rcpp::is_na(idxin3) == false]; // rsname in both panel and ss in the order of panel.
  IntegerVector bp3_4use = bp3[Rcpp::is_na(idxin3) == false]; // bp3 for block diagonal
  IntegerVector chr3_4use = chr3[Rcpp::is_na(idxin3) == false]; // chr3 for block diagonal
  // cout << "size of intersect: " << rsname_4use.size() << endl;
  
  cout << endl;
  //----------------------------------------------------------------------------------------------//
  printf("Start matching SNPs with summary stat for exposure. \n");
  //match snps (rsname_4use; rsname1: trait)
  IntegerVector idxin1 = match(rsname_4use, rsname1);  //index for SNPs in ss
  // cout << "size of summarystat 1:" << idxin1.size() << endl;
  uvec idx = as<uvec>(idxin1) -1;
  fvec tmp = as<fvec>(betah1);
  fvec bhat1 = tmp.elem(idx);
  tmp = as<fvec>(s12);
  fvec se12 = tmp.elem(idx);
  
  tmp = as<fvec>(sample_size);
  fvec samsize = tmp.elem(idx);
  uvec tmp2 = as<uvec>(SA11);
  uvec SA11_ = tmp2.elem(idx);
  tmp2 = as<uvec>(SA12);
  uvec SA12_ = tmp2.elem(idx);
  //----------------------------------------------------------------------------------------------//
  printf("Start matching SNPs with summary stat for outcome. \n");
  //match snps (rsname_4use; rsname1: trait)
  IntegerVector idxin2 = match(rsname_4use, rsname2);  //index for SNPs in ss
  
  idx = as<uvec>(idxin2) -1;
  tmp = as<fvec>(betah2);
  fvec bhat2 = tmp.elem(idx);
  tmp = as<fvec>(s22);
  fvec se22 = tmp.elem(idx);
  
  // tmp = as<fvec>(sample_size);
  // fvec samsize = tmp.elem(idx);
  tmp2 = as<uvec>(SA21);
  uvec SA21_ = tmp2.elem(idx);
  tmp2 = as<uvec>(SA22);
  uvec SA22_ = tmp2.elem(idx);
  //----------------------------------------------------------------------------------------------//
  //match snps (rsname_4use; rsname3: panel)
  printf("Start matching SNPs with panel. \n");
  cout << endl;
  idxin3 = match(rsname_4use, rsname3); //index for SNPs in panel SNPs
  //----------------------------------------------------------------------------------------------//
  //compare direction
  // cout << "Size of matched SNPs: " << rsname_4use.size() << endl; //compare keepIndx in R: Height_1000Q_match.R
  idx = as<uvec>(idxin3) -1;
  //mat R = tmp1.rows(idx);
  tmp2 = as<uvec>(A31);
  uvec A31_ = tmp2.elem(idx);
  tmp2 = as<uvec>(A32);
  uvec A32_ = tmp2.elem(idx);
  //----------------------------------------------------------------------------------------------//
  //ascii: (A:65; C:67; G:71; T:84) (a:97;c:99,g:103;t:116)
  //replace lower letter to upper letter in SA11_, SA12_, SA11_ and SA12_;
  idx = find(SA11_ > 85);
  SA11_.elem(idx) = SA11_.elem(idx) - 32;
  idx = find(SA12_ > 85);
  SA12_.elem(idx) = SA12_.elem(idx) - 32;
  idx = find(SA21_ > 85);
  SA21_.elem(idx) = SA21_.elem(idx) - 32;
  idx = find(SA22_ > 85);
  SA22_.elem(idx) = SA22_.elem(idx) - 32;
  idx = find(A31_ > 85);
  A31_.elem(idx) = A31_.elem(idx) - 32;
  idx = find(A32_ > 85);
  A32_.elem(idx) = A32_.elem(idx) - 32;
  // cout << "check error 2" << endl;
  //compare A31_ SA11_, A32_ SA12_
  //A31_: replace T with A,replace G with C
  idx = find(A31_ == 84);
  uvec idxrepl1(idx.n_elem);
  idxrepl1.fill(65);
  A31_.elem(idx) = idxrepl1;
  
  idx = find(A31_ == 71);
  uvec idxrepl2(idx.n_elem);
  idxrepl2.fill(67);
  A31_.elem(idx) = idxrepl2;
  
  //SA11_: replace T with A,replace G with C
  idx = find(SA11_ == 84);
  uvec idxrepl3(idx.n_elem);
  idxrepl3.fill(65);
  SA11_.elem(idx) = idxrepl3;
  
  idx = find(SA11_ == 71);
  uvec idxrepl4(idx.n_elem);
  idxrepl4.fill(67);
  SA11_.elem(idx) = idxrepl4;
  
  //A32_: replace T with A,replace G with C
  idx = find(A32_ == 84);
  uvec idxrepl5(idx.n_elem);
  idxrepl5.fill(65);
  A32_.elem(idx) = idxrepl5;
  
  idx = find(A32_ == 71);
  uvec idxrepl6(idx.n_elem);
  idxrepl6.fill(67);
  A32_.elem(idx) = idxrepl6;
  
  //SA12_: replace T with A,replace G with C
  idx = find(SA12_ == 84);
  uvec idxrepl7(idx.n_elem);
  idxrepl7.fill(65);
  SA12_.elem(idx) = idxrepl7;
  
  idx = find(SA12_ == 71);
  uvec idxrepl8(idx.n_elem);
  idxrepl8.fill(67);
  SA12_.elem(idx) = idxrepl8;
  
  //SA21_: replace T with A,replace G with C
  idx = find(SA21_ == 84);
  uvec idxrepl9(idx.n_elem);
  idxrepl9.fill(65);
  SA21_.elem(idx) = idxrepl9;
  
  idx = find(SA21_ == 71);
  uvec idxrepl10(idx.n_elem);
  idxrepl10.fill(67);
  SA21_.elem(idx) = idxrepl10;
  
  //SA22_: replace T with A,replace G with C
  idx = find(SA22_ == 84);
  uvec idxrepl11(idx.n_elem);
  idxrepl11.fill(65);
  SA22_.elem(idx) = idxrepl11;
  
  idx = find(SA22_ == 71);
  uvec idxrepl12(idx.n_elem);
  idxrepl12.fill(67);
  SA22_.elem(idx) = idxrepl12;
  // cout << "SA22_:" << SA22_.size() << endl;
  
  
  // #############################################################
  //remove index which are ambiguous;
  uvec idx1, idx2;
  idx1 = find((SA11_ + SA12_) == (A31_ + A32_));
  idx2 = find((SA11_ + SA12_) == (SA21_ + SA22_));
  // idx1 = find((A31_ + A32_) == (SA11_ + SA12_));
  // idx2 = find((A31_ + A32_) == (SA21_ + SA22_));
  idx = intersect(idx1, idx2);
  uvec tmp1 = as<uvec>(idxin3);
  uvec idxin4 = tmp1.elem(idx);// the final index.
  uvec A31_r = A31_.elem(idx), A32_r = A32_.elem(idx);
  uvec SA11_r = SA11_.elem(idx), SA12_r = SA12_.elem(idx);
  uvec SA21_r = SA21_.elem(idx), SA22_r = SA22_.elem(idx);
  
  IntegerVector idxtmp = wrap(idx);
  CharacterVector rsname_final = rsname_4use[idxtmp];
  IntegerVector bp_final = bp3_4use[idxtmp];
  IntegerVector chr_final = chr3_4use[idxtmp];
  //------------------------------------------------//
  fvec bh1, s12_, bh2, s22_;
  // the panel data comparing with exposure data.
  uvec idx4panel;
  
  if(matchExp)
  {
    // match with exposure data.
    // the exposure data
    bh1 = bhat1.elem(idx);
    s12_ = se12.elem(idx);
    //------------------------------------------------//
    // the outcome data comparing with exposure data.
    bh2 = bhat2.elem(idx);
    s22_ = se22.elem(idx);
    
    fvec ind2(SA11_r.n_elem);
    idx2 = find(SA11_r == SA21_r);
    ind2.ones();
    ind2 = -ind2;
    ind2.elem(idx2).ones();
    bh2 = bh2 % ind2;
    //------------------------------------------------//
    // return the index of matchdata whose minor allel does not equal to the exposure data.
    idx4panel = find(SA11_r != A31_r);
  }
  else
  {
    bh1 = bhat1.elem(idx);
    s12_ = se12.elem(idx);
    
    fvec ind1(A31_r.n_elem);
    idx1 = find(A31_r == SA11_r);
    ind1.ones();
    ind1 = -ind1;
    ind1.elem(idx1).ones();
    cout << "direction (compare with trait_1000Q_match.R): " << sum(ind1) << endl; //compare sum(ind) in R: Height_1000Q_match.R
    
    bh1 = bh1 % ind1;
    
    //------------------------------------------------------------------------------------//
    bh2 = bhat2.elem(idx);
    s22_ = se22.elem(idx);
    
    fvec ind2(A31_r.n_elem);
    idx2 = find(A31_r == SA21_r);
    ind2.ones();
    ind2 = -ind2;
    ind2.elem(idx2).ones();
    cout << "direction (compare with trait_1000Q_match.R): " << sum(ind2) << endl; //compare sum(ind) in R: Height_1000Q_match.R
    cout << "Size of matched SNPs (remove ambiguous SNPs): " << A31_r.n_elem << endl;
    bh2 = bh2 % ind2;
  }
  // #############################################################
  
  
  
  List output = List::create(Rcpp::Named("bh1") = bh1,
                             Rcpp::Named("bh2") = bh2,
                             Rcpp::Named("s12") = s12_,
                             Rcpp::Named("s22") = s22_,
                             Rcpp::Named("rsname") = rsname_final,
                             Rcpp::Named("bp") = bp_final,
                             Rcpp::Named("chr") = chr_final,
                             Rcpp::Named("idxin") = idxin4,
                             Rcpp::Named("idx4panel") = idx4panel);
  return output;
  
}


// [[Rcpp::export]]
Rcpp::List matchscreen(std::string screenname, std::string stringname1,
                       std::string stringname2, std::string stringname3,
                       double pva_cutoff, bool matchExp = false){
  
  // screening dataset:screenname;
  // summary data 1(exposure): stringname1;
  // summary data 2(outcome): stringname2;
  // Panel SNP info:stringname3;
  //----------------------------------------------------------------------------------------------//
  // summary stat SNP info for trait 1, betah, s2, sample_size
  int Nscreen = getLineNum(screenname) - 1;// ignore the colnames line
  
  IntegerVector SA01(Nscreen), SA02(Nscreen), chr0(Nscreen), bp0(Nscreen);
  CharacterVector rsname0(Nscreen);
  NumericVector betah0(Nscreen), s02(Nscreen), PV0(Nscreen), sample_size(Nscreen);
  
  Read_summarystat(screenname, SA01, SA02, rsname0, betah0, s02, PV0, chr0, bp0, Nscreen);
  LogicalVector idx_cutoff = PV0 < pva_cutoff;
  rsname0 = rsname0[idx_cutoff  == 1];
  //----------------------------------------------------------------------------------------------//
  // summary stat SNP info for trait 1, betah, s2, sample_size
  int Ns1 = getLineNum(stringname1) - 1;// ignore the colnames line
  
  IntegerVector SA11(Ns1), SA12(Ns1), chr1(Ns1), bp1(Ns1);
  CharacterVector rsname1(Ns1);
  NumericVector betah1(Ns1), s12(Ns1), PV1(Ns1);
  
  Read_summarystat(stringname1, SA11, SA12, rsname1, betah1, s12, PV1, chr1, bp1, Ns1);
  cout << "rsname1.size()" <<rsname1.size()<<endl;
  
  CharacterVector scre_rs1 = intersect(rsname0, rsname1);
  IntegerVector scre_idx1 = match(rsname1, scre_rs1);
  rsname1 = rsname1[Rcpp::is_na(scre_idx1) == false];
  SA11 = SA11[Rcpp::is_na(scre_idx1) == false];
  SA12 = SA12[Rcpp::is_na(scre_idx1) == false];
  betah1 = betah1[Rcpp::is_na(scre_idx1) == false];
  s12 = s12[Rcpp::is_na(scre_idx1) == false];
  PV1 = PV1[Rcpp::is_na(scre_idx1) == false];
  chr1 = chr1[Rcpp::is_na(scre_idx1) == false];
  bp1 = bp1[Rcpp::is_na(scre_idx1) == false];
  cout << "rsname1.size()" <<rsname1.size()<<endl;
  //  //----------------------------------------------------------------------------------------------//
  // summary stat SNP info for trait 2, betah, s2, sample_size
  int Ns2 = getLineNum(stringname2) - 1; // ignore the colnames line
  
  IntegerVector SA21(Ns2), SA22(Ns2), chr2(Ns2), bp2(Ns2);
  CharacterVector rsname2(Ns2);
  NumericVector betah2(Ns2), s22(Ns2), PV2(Ns2);
  
  
  Read_summarystat(stringname2, SA21, SA22, rsname2, betah2, s22, PV2, chr2, bp2, Ns2);
  CharacterVector scre_rs2 = intersect(rsname0, rsname2);
  IntegerVector scre_idx2 = match(rsname2, scre_rs2);
  rsname2 = rsname2[Rcpp::is_na(scre_idx2) == false];
  SA21 = SA21[Rcpp::is_na(scre_idx2) == false];
  SA22 = SA22[Rcpp::is_na(scre_idx2) == false];
  betah2 = betah2[Rcpp::is_na(scre_idx2) == false];
  s22 = s22[Rcpp::is_na(scre_idx2) == false];
  PV2 = PV2[Rcpp::is_na(scre_idx2) == false];
  chr2 = chr2[Rcpp::is_na(scre_idx2) == false];
  bp2 = bp2[Rcpp::is_na(scre_idx2) == false];
  
  // //----------------------------------------------------------------------------------------------//
  string bimfile = stringname3;
  bimfile += ".bim";
  
  int N = getLineNum(bimfile);
  cout << "Number of SNPs (panel):" << N << endl;
  
  
  IntegerVector A31(N), A32(N);
  CharacterVector rsname3(N);
  IntegerVector chr3(N), bp3(N);
  NumericVector morgan(N);
  
  //read from SNPinfo file (pass as pointer)
  cout << endl;
  cout << "Start loading SNP info:" << endl;
  clock_t t1 = clock();
  ReadSNPinfo(bimfile, A31, A32, rsname3, chr3, bp3, morgan, N);
  cout << "Finish loading SNP info for panel data in " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec." << endl;
  
  //----------------------------------------------------------------------------------------------//
  //mathcing panel SNP with summary stat and correct direction of minor allele
  CharacterVector rs_tmp = intersect(rsname1, rsname2);
  CharacterVector rs_inter = intersect(rs_tmp, rsname3);
  //----------------------------------------------------------------------------------------------//
  // cout << "check error 0" <<endl;
  //Panel Data
  IntegerVector idxin3 = match(rsname3, rs_inter);
  
  CharacterVector rsname_4use = rsname3[Rcpp::is_na(idxin3) == false]; // rsname in both panel and ss in the order of panel.
  IntegerVector bp3_4use = bp3[Rcpp::is_na(idxin3) == false]; // bp3 for block diagonal
  IntegerVector chr3_4use = chr3[Rcpp::is_na(idxin3) == false]; // chr3 for block diagonal
  // cout << "size of intersect: " << rsname_4use.size() << endl;
  
  // cout << "check error 1" <<endl;
  //----------------------------------------------------------------------------------------------//
  printf("Start matching SNPs with summary stat for trait 1. \n");
  //match snps (rsname_4use; rsname1: trait)
  IntegerVector idxin1 = match(rsname_4use, rsname1);  //index for SNPs in ss
  // cout << "size of summarystat 1:" << idxin1.size() << endl;
  uvec idx = as<uvec>(idxin1) -1;
  fvec tmp = as<fvec>(betah1);
  fvec bhat1 = tmp.elem(idx);
  tmp = as<fvec>(s12);
  fvec se12 = tmp.elem(idx);
  
  tmp = as<fvec>(sample_size);
  fvec samsize = tmp.elem(idx);
  uvec tmp2 = as<uvec>(SA11);
  uvec SA11_ = tmp2.elem(idx);
  tmp2 = as<uvec>(SA12);
  uvec SA12_ = tmp2.elem(idx);
  //----------------------------------------------------------------------------------------------//
  printf("Start matching SNPs with summary stat for trait 2. \n");
  //match snps (rsname_4use; rsname1: trait)
  IntegerVector idxin2 = match(rsname_4use, rsname2);  //index for SNPs in ss
  cout << "size of summarystat 2:" << idxin2.size() << endl;
  idx = as<uvec>(idxin2) -1;
  tmp = as<fvec>(betah2);
  fvec bhat2 = tmp.elem(idx);
  tmp = as<fvec>(s22);
  fvec se22 = tmp.elem(idx);
  
  // tmp = as<fvec>(sample_size);
  // fvec samsize = tmp.elem(idx);
  tmp2 = as<uvec>(SA21);
  uvec SA21_ = tmp2.elem(idx);
  tmp2 = as<uvec>(SA22);
  uvec SA22_ = tmp2.elem(idx);
  //----------------------------------------------------------------------------------------------//
  //match snps (rsname_4use; rsname3: panel)
  printf("Start matching SNPs with panel. \n");
  idxin3 = match(rsname_4use, rsname3); //index for SNPs in panel SNPs
  //----------------------------------------------------------------------------------------------//
  //compare direction
  // cout << "Size of matched SNPs: " << rsname_4use.size() << endl; //compare keepIndx in R: Height_1000Q_match.R
  idx = as<uvec>(idxin3) -1;
  //mat R = tmp1.rows(idx);
  tmp2 = as<uvec>(A31);
  uvec A31_ = tmp2.elem(idx);
  tmp2 = as<uvec>(A32);
  uvec A32_ = tmp2.elem(idx);
  //----------------------------------------------------------------------------------------------//
  //ascii: (A:65; C:67; G:71; T:84) (a:97;c:99,g:103;t:116)
  //replace lower letter to upper letter in SA11_, SA12_, SA11_ and SA12_;
  idx = find(SA11_ > 85);
  SA11_.elem(idx) = SA11_.elem(idx) - 32;
  idx = find(SA12_ > 85);
  SA12_.elem(idx) = SA12_.elem(idx) - 32;
  idx = find(SA21_ > 85);
  SA21_.elem(idx) = SA21_.elem(idx) - 32;
  idx = find(SA22_ > 85);
  SA22_.elem(idx) = SA22_.elem(idx) - 32;
  idx = find(A31_ > 85);
  A31_.elem(idx) = A31_.elem(idx) - 32;
  idx = find(A32_ > 85);
  A32_.elem(idx) = A32_.elem(idx) - 32;
  // cout << "check error 2" << endl;
  //compare A31_ SA11_, A32_ SA12_
  //A31_: replace T with A,replace G with C
  idx = find(A31_ == 84);
  uvec idxrepl1(idx.n_elem);
  idxrepl1.fill(65);
  A31_.elem(idx) = idxrepl1;
  
  idx = find(A31_ == 71);
  uvec idxrepl2(idx.n_elem);
  idxrepl2.fill(67);
  A31_.elem(idx) = idxrepl2;
  
  //SA11_: replace T with A,replace G with C
  idx = find(SA11_ == 84);
  uvec idxrepl3(idx.n_elem);
  idxrepl3.fill(65);
  SA11_.elem(idx) = idxrepl3;
  
  idx = find(SA11_ == 71);
  uvec idxrepl4(idx.n_elem);
  idxrepl4.fill(67);
  SA11_.elem(idx) = idxrepl4;
  
  //A32_: replace T with A,replace G with C
  idx = find(A32_ == 84);
  uvec idxrepl5(idx.n_elem);
  idxrepl5.fill(65);
  A32_.elem(idx) = idxrepl5;
  
  idx = find(A32_ == 71);
  uvec idxrepl6(idx.n_elem);
  idxrepl6.fill(67);
  A32_.elem(idx) = idxrepl6;
  
  //SA12_: replace T with A,replace G with C
  idx = find(SA12_ == 84);
  uvec idxrepl7(idx.n_elem);
  idxrepl7.fill(65);
  SA12_.elem(idx) = idxrepl7;
  
  idx = find(SA12_ == 71);
  uvec idxrepl8(idx.n_elem);
  idxrepl8.fill(67);
  SA12_.elem(idx) = idxrepl8;
  
  //SA21_: replace T with A,replace G with C
  idx = find(SA21_ == 84);
  uvec idxrepl9(idx.n_elem);
  idxrepl9.fill(65);
  SA21_.elem(idx) = idxrepl9;
  
  idx = find(SA21_ == 71);
  uvec idxrepl10(idx.n_elem);
  idxrepl10.fill(67);
  SA21_.elem(idx) = idxrepl10;
  
  //SA22_: replace T with A,replace G with C
  idx = find(SA22_ == 84);
  uvec idxrepl11(idx.n_elem);
  idxrepl11.fill(65);
  SA22_.elem(idx) = idxrepl11;
  
  idx = find(SA22_ == 71);
  uvec idxrepl12(idx.n_elem);
  idxrepl12.fill(67);
  SA22_.elem(idx) = idxrepl12;
  cout << "SA22_:" << SA22_.size() << endl;
  // #############################################################
  //remove index which are ambiguous;
  uvec idx1, idx2;
  idx1 = find((SA11_ + SA12_) == (A31_ + A32_));
  idx2 = find((SA11_ + SA12_) == (SA21_ + SA22_));
  // idx1 = find((A31_ + A32_) == (SA11_ + SA12_));
  // idx2 = find((A31_ + A32_) == (SA21_ + SA22_));
  idx = intersect(idx1, idx2);
  uvec tmp1 = as<uvec>(idxin3);
  uvec idxin4 = tmp1.elem(idx);// the final index.
  uvec A31_r = A31_.elem(idx), A32_r = A32_.elem(idx);
  uvec SA11_r = SA11_.elem(idx), SA12_r = SA12_.elem(idx);
  uvec SA21_r = SA21_.elem(idx), SA22_r = SA22_.elem(idx);
  
  IntegerVector idxtmp = wrap(idx);
  CharacterVector rsname_final = rsname_4use[idxtmp];
  IntegerVector bp_final = bp3_4use[idxtmp];
  IntegerVector chr_final = chr3_4use[idxtmp];
  //------------------------------------------------//
  fvec bh1, s12_, bh2, s22_;
  // the panel data comparing with exposure data.
  uvec idx4panel;
  if(matchExp)
  {
    // match with exposure data.
    // the exposure data
    bh1 = bhat1.elem(idx);
    s12_ = se12.elem(idx);
    //------------------------------------------------//
    // the outcome data comparing with exposure data.
    bh2 = bhat2.elem(idx);
    s22_ = se22.elem(idx);

    fvec ind2(SA11_r.n_elem);
    idx2 = find(SA11_r == SA21_r);
    ind2.ones();
    ind2 = -ind2;
    ind2.elem(idx2).ones();
    bh2 = bh2 % ind2;
    //------------------------------------------------//
    // return the index of matchdata whose minor allel does not equal to the exposure data.
    idx4panel = find(SA11_r != A31_r);
  }
  else
  {
    bh1 = bhat1.elem(idx);
    s12_ = se12.elem(idx);

    fvec ind1(A31_r.n_elem);
    idx1 = find(A31_r == SA11_r);
    ind1.ones();
    ind1 = -ind1;
    ind1.elem(idx1).ones();
    cout << "direction (compare with trait_1000Q_match.R): " << sum(ind1) << endl; //compare sum(ind) in R: Height_1000Q_match.R
    
    bh1 = bh1 % ind1;
    
    //------------------------------------------------------------------------------------//
    bh2 = bhat2.elem(idx);
    s22_ = se22.elem(idx);

    fvec ind2(A31_r.n_elem);
    idx2 = find(A31_r == SA21_r);
    ind2.ones();
    ind2 = -ind2;
    ind2.elem(idx2).ones();
    cout << "direction (compare with trait_1000Q_match.R): " << sum(ind2) << endl; //compare sum(ind) in R: Height_1000Q_match.R
    cout << "Size of matched SNPs (remove ambiguous SNPs): " << A31_r.n_elem << endl;
    bh2 = bh2 % ind2;
  }
  

  List output = List::create(Rcpp::Named("bh1") = bh1,
                             Rcpp::Named("bh2") = bh2,
                             Rcpp::Named("s12") = s12_,
                             Rcpp::Named("s22") = s22_,
                             Rcpp::Named("rsname") = rsname_final,
                             Rcpp::Named("bp") = bp_final,
                             Rcpp::Named("chr") = chr_final,
                             Rcpp::Named("idxin") = idxin4,
                             Rcpp::Named("idx4panel") = idx4panel);
  
  return output;
  
}
