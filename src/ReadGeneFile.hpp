#ifndef ReadGeneFile_hpp
#define ReadGeneFile_hpp

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <bitset>
#include <math.h>

#include "time.h"
#include <iostream>

using namespace Rcpp;
using namespace arma;
using namespace std;

int getLineNum(std::string filename);

void getFourGentype(int* geno, std::bitset<8> bits);

void vec_sub(arma::vec& v0, double m_value, double squaresum);

void readPlink(std::string stringname, int N, int P, unsigned* X);

// Rcpp::List ReadSNPinfo(std::string stringname, IntegerVector A1, IntegerVector A2, CharacterVector rsname,
//                        IntegerVector chr, IntegerVector bp, NumericVector morgan, int N);
// Rcpp::List Read_summarystat(std::string stringname, IntegerVector SA1, IntegerVector SA2, CharacterVector rsname,
//                             NumericVector betah, NumericVector s2, NumericVector pvalue, IntegerVector chr, IntegerVector bp, int N)；
// CharacterVector select(CharacterVector vec_, NumericVector idx_);
// Rcpp::List matchscreen(std::string screenname, std::string stringname1,
//                        std::string stringname2, std::string stringname3,
//                        double pva_cutoff);
// Rcpp::List matchsnp(std::string stringname1, std::string stringname2, std::string stringname3);
// arma::Col<int> get_interval(arma::Col<int> index, uword i);

arma::Col<int> get_interval(arma::Col<int> index, uword i);
#endif /* ReadGeneFile_hpp */