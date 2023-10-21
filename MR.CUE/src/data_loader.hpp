#ifndef data_loader_hpp
#define data_loader_hpp

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "time.h"
#include <iostream>
#include "ReadGeneFile.hpp"

using namespace Rcpp;
using namespace arma;
using namespace std;

Rcpp::List matchsnp(std::string stringname1, std::string stringname2, std::string stringname3, bool matchExp = false);
Rcpp::List matchscreen(std::string screenname, std::string stringname1, std::string stringname2, std::string stringname3, double pva_cutoff, bool matchExp = false);
#endif /* data_loader_hpp */