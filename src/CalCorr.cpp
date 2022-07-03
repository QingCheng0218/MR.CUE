#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdio.h>
#include <bitset>
#include <math.h>
#include <boost/algorithm/string.hpp>
#include <boost/range.hpp>
#include <boost/range/algorithm.hpp>
#include "time.h"
#include <iostream>
#include <thread>
#include <vector>
#include "ReadGeneFile.hpp"



using namespace Rcpp;
using namespace arma;
using namespace std;
using namespace boost;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

// --------------------------------------------------------
// New change: related 'cal_blocks' function (Apr 11, 2021);
// New change: remove the package pdsce. (Aug 24, 2021);
// --------------------------------------------------------
class DataBlock{
public:
  int chrom;
  int block_no;
  int start;
  int end;
  DataBlock(int chrom, int block_no);
  DataBlock();
  DataBlock(const DataBlock& block);
  fmat corr;
  double heritability_beta;
  double heritability;
};

DataBlock::DataBlock(int chrom, int block_no){
  this->chrom = chrom;
  this->block_no = block_no;
  this->heritability_beta = 0;
}

DataBlock::DataBlock(){

}

DataBlock::DataBlock(const DataBlock& block){
  this->chrom = block.chrom;
  this->block_no = block.block_no;
  this->start = block.start;
  this->end = block.end;
}


// [[Rcpp::export]]
vector<umat> load_block_file(string block_file){
  ivec chrom_vec;

  vector<umat> blocks_mat;
  vector<int *> blocks;
  std::ifstream ifs(block_file.c_str());

  std::string line;
  int chromsome;
  string snpname;
  vector <string> fields;
  int i = 0;
  std::getline(ifs, line);

  while (std::getline(ifs, line)) // read one line from ifs
  {
    std::istringstream iss(line); // access line as a stream
    i++;
  }

  ifs.close();
  ifs.open(block_file.c_str());
  std::getline(ifs, line);
  chrom_vec.resize(i);
  umat block_mat(i, 2);
  i = 0;
  while (std::getline(ifs, line)) // read one line from ifs
  {
    std::istringstream iss(line); // access line as a stream
    boost::algorithm::split(fields, line, boost::is_any_of("\t"));
    chromsome = atoi(fields[0].substr(3).c_str());
    chrom_vec[i] = chromsome - 1;
    block_mat(i, 0) = atof(fields[1].c_str());
    block_mat(i, 1) = atof(fields[2].c_str());
    i++;
  }

  for (int k = 0; k < 22; k++){
    blocks_mat.push_back(block_mat.rows(find(chrom_vec == k)));
  }
  ifs.close();
  return blocks_mat;
}

void cal_blocks(vector<DataBlock> &datablocks, arma::ivec bp, arma::ivec chr, std::string block_file){
  vector<arma::fmat> corr_blocks;
  vector<umat>  blocks = load_block_file(block_file.c_str());

  Col<int> chrom_list;
  chrom_list = arma::unique(chr);
  ivec chr_index(chrom_list.size());

  double sum2 = 0;
  for(int i = 1; i<=(int)(chrom_list.size()); i++){
    int chr_i = chrom_list[i - 1];
    chr_index[i - 1] = sum(chr==chr_i);
    // chr_index[i - 1] = sum(chr==i);
    sum2 += chr_index[i - 1];
  }

  ivec tmpsum = cumsum(chr_index);
  uword n_chrom = chr_index.size();
  vector<int> block_vector;
  int last_block_no = -1;

  for (uword i = 0; i < n_chrom; i++) {
    ivec u_i = get_interval(tmpsum, i);
    int chri_idx = chrom_list[i] - 1;
    umat block_i = blocks[chri_idx];
    // umat block_i = blocks[i];
    uword n_block = block_i.n_rows;
    int last_chrom;

    for (int kk = 0; kk < (int)(u_i.size()); kk++){
      bool in = false;
      int block_no = -1;
      int k = u_i[kk];
      for (int j = 0; j < (int)(n_block); j++){
        in = (bp[k] <= (int)(block_i(j, 1))) && (bp[k] >= (int)(block_i(j, 0)));
        if (in){
          block_no = j;
          break;
        }
      }
      if (block_no != -1){
        if (block_no < last_block_no){
          // int last_chrom = datablocks[datablocks.size() - 1].chrom;
          if (last_chrom == (int)(i)){
            cout << "Error Order" << endl;
          }
        }
        if (block_no != last_block_no && last_block_no != -1){
          DataBlock block((int)i, last_block_no);
          block.start = block_vector[0];
          block.end = block_vector[block_vector.size() - 1];
          datablocks.push_back(block);
          block_vector.clear();
        }
        
        // In case that the last location of last chrom equals to the current location,
        // they may combine the blocks, we wanna start a new block (different chroms should not be assigned in the same block).
        if(block_no == last_block_no&&kk==0){
          // cout << "if change blocks!!" << endl;
          DataBlock block((int)i, last_block_no);
          block.start = block_vector[0];
          block.end = block_vector[block_vector.size() - 1];
          datablocks.push_back(block);
          block_vector.clear();
        }
        
        block_vector.push_back(k);
        last_block_no = block_no;
      }
      last_chrom = i;
    }
  }

  if (block_vector.size() > 0){
    DataBlock block((int)(n_chrom - 1), last_block_no);
    block.start = block_vector[0];
    block.end = block_vector[block_vector.size() - 1];
    datablocks.push_back(block);
  }

  // Col<int> chrom_vec;
  // chrom_vec.resize(datablocks.size());
  // chrom_vec.zeros();
  // for (int i = 0; i < (int)(datablocks.size()); i++){
  //   int chrom_idx = datablocks[i].chrom;
  //   uvec idx_u = find(chrom_list == chrom_idx + 1);
  //   uword idx = idx_u[0];
  //   chrom_vec[i] = (int)idx;
  // }

}
// [[Rcpp::export]]
List test_blocks(arma::ivec bp, arma::ivec chr, std::string block_file){
  vector<DataBlock> datablocks;
  vector<arma::fmat> corr_blocks;
  vector<umat>  blocks = load_block_file(block_file.c_str());

  Col<int> chrom_list;

  chrom_list = arma::unique(chr);
  ivec chr_index(chrom_list.size());

  double sum2 = 0;
  for(int i = 1; i<=(int)(chrom_list.size()); i++){
    int chr_i = chrom_list[i - 1];
    chr_index[i - 1] = sum(chr==chr_i);
    // chr_index[i - 1] = sum(chr==i);
    sum2 += chr_index[i - 1];
  }
  // cout<< "chr_index:"<<chr_index.t()<<endl;
  ivec tmpsum = cumsum(chr_index);
  uword n_chrom = chr_index.size();
  vector<int> block_vector;
  int last_block_no = -1;
  // cout << "check error 0 "<< endl;
  // cout << "n_chrom: "<< n_chrom << endl;
  // cout <<"tmpsum:"<< tmpsum<<endl;


  for (uword i = 0; i < n_chrom; i++) {
    ivec u_i = get_interval(tmpsum, i);
    // cout << "check error" <<endl;
    // cout <<"u_i:"<<u_i.t()<<endl;
    int chri_idx = chrom_list[i] - 1;
    umat block_i = blocks[chri_idx];
    uword n_block = block_i.n_rows;
    int last_chrom;

    for (int kk = 0; kk < (int)(u_i.size()); kk++){
      bool in = false;
      int block_no = -1;
      int k = u_i[kk];
      for (int j = 0; j < (int)(n_block); j++){
        in = (bp[k] <= (int)(block_i(j, 1))) && (bp[k] >= (int)(block_i(j, 0)));
        if (in){
          block_no = j;
          break;
        }
      }
      if (block_no != -1){
        if (block_no < last_block_no){
          if (last_chrom == chri_idx){
            cout << "Error Order" << endl;
          }
        }
        if (block_no != last_block_no && last_block_no != -1){
          DataBlock block((int)i, last_block_no);
          block.start = block_vector[0];
          block.end = block_vector[block_vector.size() - 1];
          datablocks.push_back(block);
          block_vector.clear();
        }
        block_vector.push_back(k);
        last_block_no = block_no;
      }
      last_chrom = i;

    }
  }

  // cout << "check error 1"<< endl;
  if (block_vector.size() > 0){
    DataBlock block((int)(n_chrom - 1), last_block_no);
    block.start = block_vector[0];
    block.end = block_vector[block_vector.size() - 1];
    datablocks.push_back(block);
  }


  Col<int> chrom_vec;
  chrom_vec.resize(datablocks.size());
  chrom_vec.zeros();


  // for (int i = 0; i < (int)(datablocks.size()); i++){
  //   int chrom_idx = datablocks[i].chrom;
  //   uvec idx_u = find(chrom_list == chrom_idx + 1);
  //   // cout <<"i:" << i << "idx_u: " << idx_u<<endl;
  //   cout<<"chrom_idx :"<<chrom_idx <<endl;
  //   uword idx = idx_u[0];
  //   chrom_vec[i] = (int)idx;
  // }


  // -------------------------------------------------------
  // cout << "check error2 "<< endl;
  uword nblocks = datablocks.size();
  umat block_inf = zeros<umat>(nblocks, 2);
  for(int ii = 0; ii<(int)(nblocks); ii++){
    block_inf(ii, 0) = datablocks[ii].start;
    block_inf(ii, 1) = datablocks[ii].end;
  }
  List output = List::create(
    Rcpp::Named("block_inf") = block_inf
  );
  return output;

}



// [[Rcpp::export]]
arma::ivec std_setdiff(arma::ivec& x, arma::ivec& y) {
  
  // std::vector<int> a = arma::conv_to< std::vector<int> >::from(arma::sort(x));
  // std::vector<int> b = arma::conv_to< std::vector<int> >::from(arma::sort(y));
  std::vector<int> out;
  
  std::set_difference(x.begin(), x.end(), y.begin(), y.end(),
                      std::inserter(out, out.end()));
  
  // return arma::conv_to<arma::ivec>::from(out);
  return out;
}
// [[Rcpp::export]]
ivec LDclump(arma::mat &R, double ld_r2_thresh){
  int nsnp = R.n_rows;
  
  R.diag() = zeros(nsnp);
  R = trimatl(R);
  
  ivec a = regspace<ivec>(0, 1, nsnp - 1);
  ivec idx1 = zeros<ivec>(nsnp*nsnp, 1);
  ivec idx2 = zeros<ivec>(nsnp*nsnp, 1);
  for(int i=0; i< nsnp; i++){
    uvec ii = regspace<uvec>(i*nsnp, 1, ((i+1)*nsnp - 1));
    idx1.elem(ii) = a;
    idx2.elem(ii) = i*ones<ivec>(nsnp, 1);
  }
  
  uvec q1 = find(R%R > ld_r2_thresh);
  ivec id1, id2;
  
  ivec IDsave;
  
  if(q1.n_elem){
    id1 = idx2.elem(q1);
    id2 = idx1.elem(q1);
    
    idx1.reset();
    idx2.reset();
    
    ivec slct_indx = join_cols(id1, id2);
    ivec indxbuf = arma::unique(slct_indx);
    
    // count how many SNPs in high LD
    int nproc = indxbuf.n_elem;
    ivec n_slct_snp = zeros<ivec>(nproc, 1);
    
    for( int i = 0; i < nproc; i = i + 1 ){
      n_slct_snp[i] = sum(slct_indx==indxbuf[i]);
    }
    
    // decide the index to remove
    nproc = id1.n_elem;
    int n1, n2;
    for(int i = 0; i < nproc; i = i + 1){
      n1 = n_slct_snp[as_scalar(find(indxbuf==id1[i]))];
      n2 = n_slct_snp[as_scalar(find(indxbuf==id2[i]))];
      if(n1 < n2){
        int t;
        t = id1[i];
        id1[i] = id2[i];
        id2[i] = t;
      }
      
    }
    
    id1 = arma::unique(id1);
    // return the index to save(C indicator start from 0)
    IDsave = std_setdiff(a, id1);
  }else{
    IDsave = a;
    
  }
  
  return IDsave;
  
}

class paraBlock_CorR{

public:
  int current_idx=0;
  uword Ngene_active;
  int n_thread = 1;
  double ld_r2_thresh;
  vector<DataBlock> datablocks;
  arma::Mat<unsigned>* X;
  // Chroms* chroms;
  Col<int> chrom_vec;

  uword nblocks;
  field<mat> F4Rblock;
  field<ivec> F4index;
  arma::ivec bp;
  arma::ivec chr;
  arma::uvec avbIndex;
  std::string block_file;
  std::string stringname3;
  double lam;
  arma::uvec Nb;
  arma::uvec Nidex;


  paraBlock_CorR(uword &nblocks, vector<DataBlock> &datablocks, arma::Mat<unsigned>* X, field<mat>& F4Rblock, 
                 field<ivec>& F4index, double &ld_r2_thresh, arma::uvec &Nb, arma::uvec &Nidex, double &lam){
    this -> nblocks = nblocks;
    this -> datablocks = datablocks;
    this -> X = X;
    this -> F4Rblock = F4Rblock;
    this -> F4index= F4index;
    this -> ld_r2_thresh = ld_r2_thresh;
    this -> Nb = Nb;
    this -> Nidex = Nidex;
    this -> lam = lam;
  }

  arma::fmat calCorr();
  arma::mat cal_blockcor(Mat<unsigned>& X);

  int  next_gibbs();
  void update_by_thread_gibbs(int thread_id);
  mat loop_by_block_CarCorr(int i);
  void update_by_thread_CarCorr(int thread_id);


};

arma::mat paraBlock_CorR::cal_blockcor(Mat<unsigned>& X){
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
std::mutex _mtx1;
int paraBlock_CorR::next_gibbs(){
  std::lock_guard<std::mutex> lockGuard(_mtx1);
  if(current_idx >= (int)nblocks){
    return -1;
  }
  current_idx++;
  return current_idx-1;
}

mat paraBlock_CorR::loop_by_block_CarCorr(int i){
  Mat<unsigned> subX = X->cols(datablocks[i].start, datablocks[i].end);
  // arma::mat corr0 = cal_blockcor(subX);
  arma::mat corr = cal_blockcor(subX);
  // arma::mat corr;
  int p1 = corr.n_rows;
  // arma::mat LAM = zeros(p1, p1);
  // LAM.fill(lam);
  // 
  // if(lam > 0.5)
  // {
  //   corr = corr0;
  //   arma::mat eyeI(p1, p1);
  //   eyeI.eye();
  //   corr = corr0;
  //   corr *= lam;
  //   corr += (1 - lam)*eyeI;
  // }
  // else
  // {
  //   pdsoftObj out = pdsoft(corr0, LAM);
  //   corr = out.theta;
  // }
  
  arma::mat eyeI(p1, p1);
  eyeI.eye();
  // corr = corr0;
  corr *= lam;
  corr += (1 - lam)*eyeI;

  F4Rblock(i, 0) = corr;
  ivec id;
  id = LDclump(corr, ld_r2_thresh);
  Nidex[i] = id.n_elem;
  
  if(id.n_elem!=0){
    F4index(i, 0) = id;
  }
  
  Nb(i) = datablocks[i].end - datablocks[i].start + 1;

  return(corr);
}

void paraBlock_CorR::update_by_thread_CarCorr(int thread_id){
  while(true){
    int idx = next_gibbs();
    if(idx == -1){
      break;
    }
    loop_by_block_CarCorr(idx);
  }
}

// [[Rcpp::export]]
List Cal_blockR(arma::ivec &bp, arma::ivec &chr, arma::uvec &avbIndex, arma::uvec &idx4panel, std::string block_file,
                std::string stringname3,  double ld_r2_thresh, int coreNum, double lam){

  vector<DataBlock> datablocks;
  cout << "Start partition the block:" << endl;
  clock_t t1 = clock();

  cal_blocks(datablocks, bp, chr, block_file);
  
  cout << "Start partition the block: in " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec." << endl;
  cout << endl;
  
  cout << "Start loading the pannel data:" << endl;
  t1 = clock();
  cout<<"check error 00"<<endl;

  string famfile = stringname3;
  famfile +=".fam";
  string bimfile = stringname3;
  bimfile += ".bim";
  cout<<"check error 01"<<endl;
  int N = getLineNum(famfile);
  int P = getLineNum(bimfile);
  long long Nlong = (long long)N;
  long long Plong = (long long)P;
  cout<<"check error 02"<<"N:"<<Nlong << "P:" <<Plong << endl;
  unsigned* X0 = new unsigned[Nlong * Plong];
  cout<<"check error 03"<<"N:"<<N << "P:" <<P<<endl;
  readPlink(stringname3, N, P, X0);
  cout<<"check error 04"<<endl;
  arma::Mat<unsigned>* Xdata =  new arma::Mat<unsigned>(X0, N, P, false, false);
  //
  uword M = Xdata -> n_cols;
  arma::Mat<unsigned>* X = Xdata;
  arma::Mat<unsigned> subX0;
  if(avbIndex.size() != 0 && avbIndex.size() < M){
    subX0 = Xdata->cols(avbIndex);
    X = &subX0;
  }
  
  // revise the genotype data for matching with exposure.
  if(idx4panel.n_elem > 0){
    // for(int i = 0; i < idx4panel.n_elem; i++){
    //   X->cols(i) = 2 - X->cols(idx4panel.elem(i));
    // }
    cout <<"Revise the genotype data for matching with exposure." << endl;
    
    X->cols(idx4panel) = 2 - X->cols(idx4panel);
  }

  uword nblocks = datablocks.size();
  uvec Nb;
  Nb.zeros(nblocks, 1);

  // block_inf: save the start and end information of datablocks
  umat block_inf = zeros<umat>(nblocks, 2);
  for(int ii = 0; ii<(int)(nblocks); ii++){
    block_inf(ii, 0) = datablocks[ii].start;
    block_inf(ii, 1) = datablocks[ii].end;
  }

  
  field<mat> F4Rblock(nblocks, 1);
  field<ivec> F4index(nblocks, 1);
  uvec Nidex = zeros<uvec>(nblocks, 1);
  // umat block_inf = zeros<umat>(nblocks, 2);
  // uvec Nb = zeros<uvec>(nblocks, 1);
  
  cout << "Start caculating the correlation between SNPs :" << endl;
  t1 = clock();
  // -----------------------------------------------------------------------
  // parallel for correlation matrix
  // set parallele structure object
  paraBlock_CorR parobj(nblocks, datablocks, X, F4Rblock, F4index,  ld_r2_thresh, Nb, Nidex, lam);
  //set parallel computation
  const int n_thread = coreNum;
  std::vector<std::thread> threads(n_thread);
  for(int i_thread = 0; i_thread < n_thread; i_thread++){
    threads[i_thread] = std::thread(&paraBlock_CorR::update_by_thread_CarCorr, &parobj, i_thread);
  }
  
  for(int i = 0; i < n_thread; i++){
    threads[i].join();
  }
  
  F4Rblock = parobj.F4Rblock;
  F4index = parobj.F4index;
  Nidex = parobj.Nidex;
  Nb = parobj.Nb;
  cout << "Finish caculating the correlation between SNPs in " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec." << endl;
  cout << endl;
  // ---------------------------------------------------------------------------
  // cout << "Start caculating the correlation of overlap:" << endl;
  // t1 = clock();
  ivec Indpid = zeros<ivec>(sum(Nidex), 1);
  uvec Nsumidex = cumsum<uvec>(Nidex);
  uvec Nbsum = cumsum<uvec>(Nb);
  uvec tmpid;
  
  
  for(int ii = 0; ii<(int)(nblocks); ii++){
    if(Nidex[ii]!=0){
      if(ii==0){
        tmpid = arma::regspace<uvec>(0, 1, Nsumidex[ii] - 1);
        Indpid.elem(tmpid) = F4index(ii, 0);
      }else{
        tmpid = arma::regspace<uvec>(Nsumidex[ii - 1], 1, Nsumidex[ii] - 1);
        Indpid.elem(tmpid) = Nbsum[ii - 1] + F4index(ii, 0);
        
      }
    }
  }
  List output = List::create(
    // Rcpp::Named("genotype") = subX0,
    // Rcpp::Named("R") = R,
    Rcpp::Named("Indpid") = Indpid,
    Rcpp::Named("block_inf") = block_inf,
    Rcpp::Named("F4Rblock") = F4Rblock,
    Rcpp::Named("Nb") = parobj.Nb,
    Rcpp::Named("nblocks") = nblocks
  );
  return output;
}



// [[Rcpp::export]]
List Cal_block_Rmatrix(arma::ivec &bp, arma::ivec &chr, arma::uvec &avbIndex, arma::uvec &idx4panel, std::string block_file,
                       std::string stringname3,  double ld_r2_thresh, int coreNum, double lam){
  
  vector<DataBlock> datablocks;
  
  cal_blocks(datablocks, bp, chr, block_file);
  
  string famfile = stringname3;
  famfile +=".fam";
  string bimfile = stringname3;
  bimfile += ".bim";
  int N = getLineNum(famfile);
  int P = getLineNum(bimfile);
  long long Nlong = (long long)N;
  long long Plong = (long long)P;
  unsigned* X0 = new unsigned[Nlong * Plong];
  readPlink(stringname3, N, P, X0);
  arma::Mat<unsigned>* Xdata =  new arma::Mat<unsigned>(X0, N, P, false, false);
  //
  uword M = Xdata -> n_cols;
  arma::Mat<unsigned>* X = Xdata;
  arma::Mat<unsigned> subX0;
  if(avbIndex.size() != 0 && avbIndex.size() < M){
    subX0 = Xdata->cols(avbIndex);
    X = &subX0;
  }
  
  // revise the genotype data for matching with exposure.
  if(idx4panel.n_elem > 0){
    
    X->cols(idx4panel) = 2 - X->cols(idx4panel);
  }
  
  uword nblocks = datablocks.size();
  uvec Nb;
  Nb.zeros(nblocks, 1);
  
  // block_inf: save the start and end information of datablocks
  umat block_inf = zeros<umat>(nblocks, 2);
  for(int ii = 0; ii<(int)(nblocks); ii++){
    block_inf(ii, 0) = datablocks[ii].start;
    block_inf(ii, 1) = datablocks[ii].end;
  }
  
  field<mat> F4Rblock(nblocks, 1);
  field<ivec> F4index(nblocks, 1);
  uvec Nidex = zeros<uvec>(nblocks, 1);
  // -----------------------------------------------------------------------
  // parallel for correlation matrix
  // set parallele structure object
  paraBlock_CorR parobj(nblocks, datablocks, X, F4Rblock, F4index,  ld_r2_thresh, Nb, Nidex, lam);
  //set parallel computation
  const int n_thread = coreNum;
  std::vector<std::thread> threads(n_thread);
  for(int i_thread = 0; i_thread < n_thread; i_thread++){
    threads[i_thread] = std::thread(&paraBlock_CorR::update_by_thread_CarCorr, &parobj, i_thread);
  }
  
  for(int i = 0; i < n_thread; i++){
    threads[i].join();
  }
  
  F4Rblock = parobj.F4Rblock;
  F4index = parobj.F4index;
  Nidex = parobj.Nidex;
  Nb = parobj.Nb;
  
  ivec Indpid = zeros<ivec>(sum(Nidex), 1);
  uvec Nsumidex = cumsum<uvec>(Nidex);
  uvec Nbsum = cumsum<uvec>(Nb);
  uvec tmpid;
  
  
  for(int ii = 0; ii<(int)(nblocks); ii++){
    if(Nidex[ii]!=0){
      if(ii==0){
        tmpid = arma::regspace<uvec>(0, 1, Nsumidex[ii] - 1);
        Indpid.elem(tmpid) = F4index(ii, 0);
      }else{
        tmpid = arma::regspace<uvec>(Nsumidex[ii - 1], 1, Nsumidex[ii] - 1);
        Indpid.elem(tmpid) = Nbsum[ii - 1] + F4index(ii, 0);
        
      }
    }
  }
  // ---------------------------------------------------------------------------
  // resize the F4Block
  mat R;
  R.zeros(avbIndex.size(), avbIndex.size());
  
  for(int i = 0; i< (int)(nblocks); i++){
    R.submat(datablocks[i].start, datablocks[i].start, datablocks[i].end, datablocks[i].end) = F4Rblock(i, 0);
  }
  
  mat IndR  = zeros(Indpid.size(), Indpid.size());
  for(int j = 0; j <(int)(Indpid.size()); j++){
    for(int i = 0; i <(int)(Indpid.size()); i++){
      IndR(i, j) = R(Indpid[i], Indpid[j]);
    }
    // vec indj = R.col(Indpid[j]);
    // 
    // // vec indj0 = indj.elem(Indpid);
    // uvec idx = as<uvec>(Indpid);
    // cout<< indj.elem(idx).t() << endl;
    // IndR.col(j) = indj.elem(Indpid);
  }

  List output = List::create(
    Rcpp::Named("IndR") = IndR,
    Rcpp::Named("Indpid") = Indpid,
    Rcpp::Named("R") = R,
    Rcpp::Named("block_inf") = block_inf,
    Rcpp::Named("F4Rblock") = F4Rblock
    // Rcpp::Named("X") = subX0
    // Rcpp::Named("nblocks") = nblocks
  );
  return output;
}





// [[Rcpp::export]]
List Cal_blockinf(arma::ivec &bp, arma::ivec &chr, std::string block_file){
  
  vector<DataBlock> datablocks;
  // clock_t t1 = clock();
  
  cal_blocks(datablocks, bp, chr, block_file);
  uword nblocks = datablocks.size();
  uvec Nb;
  Nb.zeros(nblocks, 1);
  // block_inf: save the start and end information of datablocks
  umat block_inf = zeros<umat>(nblocks, 2);
  for(int ii = 0; ii<(int)(nblocks); ii++){
    block_inf(ii, 0) = datablocks[ii].start;
    block_inf(ii, 1) = datablocks[ii].end;
  }
  
  List output = List::create(
    Rcpp::Named("block_inf") = block_inf,
    Rcpp::Named("Nb") = Nb,
    Rcpp::Named("nblocks") = nblocks
  );
  return output;
}


// [[Rcpp::export]]
Rcpp::List IndepSummary(arma::ivec &bp, arma::ivec &chr, arma::uvec &avbIndex, std::string &block_file,
                        std::string stringname3,  
                        arma::vec &bh1, arma::vec &bh2, arma::vec &se1, arma::vec &se2,
                        int coreNum,  double lam, const double &ld_r2_thresh){
  
  uvec idx4panel;
  
  cout <<"idx4panel:"<< idx4panel.n_elem<< endl;

  List Rblockres = Cal_blockR(bp, chr, avbIndex, idx4panel, block_file,
                              stringname3, ld_r2_thresh,  coreNum, lam);
  

  
  ivec Indpid = Rblockres["Indpid"];
  int inlen = Indpid.n_elem;
  arma::vec bh1_ind = zeros(inlen, 1); arma::vec bh2_ind = zeros(inlen, 1);
  arma::vec se1_ind = zeros(inlen, 1); arma::vec se2_ind = zeros(inlen, 1);
  for(int j = 0; j < inlen; j ++){
    bh1_ind[j] = bh1[Indpid[j]];
    bh2_ind[j] = bh2[Indpid[j]];
    se1_ind[j] = se1[Indpid[j]];
    se2_ind[j] = se2[Indpid[j]];
  }
  // arma::vec z1; arma::vec z2; arma::vec rhores;
  // z1 = bh1_ind / se1_ind;
  // z2 = bh2_ind / se2_ind;
  
  // double rho =overlaprho(bh1_ind, bh2_ind, se1_ind, se2_ind);
  // arma::vec a = -1*pth*ones(2, 1);
  // arma::vec b = pth*ones(2, 1);
  // rhores = truncEstfun(a, b, z1, z2, 4000, 1000, 10);
  List output = List::create(
    Rcpp::Named("bh1_ind") = Rcpp::wrap(bh1_ind),
    Rcpp::Named("bh2_ind") = Rcpp::wrap(bh2_ind),
    Rcpp::Named("se1_ind") = Rcpp::wrap(se1_ind),
    Rcpp::Named("se2_ind") = Rcpp::wrap(se2_ind)
    
  );
  return output;
  
}


// [[Rcpp::export]]
double testR(double rho, int n){
  double pvalue, tmp;
  tmp = - abs(sqrt(n-2)*(rho/sqrt(1- rho*rho)));
  pvalue = R::pt(tmp, 1.0*(n-2), true, false);
  pvalue = pvalue * 2;
  return pvalue;
}





// [[Rcpp::export]]
imat comb(int p)
{
  //2- combinations from a given set S of n elements
  vector<int> myVector;
  std::string bitmask(2, 1); // K leading 1's
  bitmask.resize(p, 0); // N-K trailing 0's
  // cout << bitmask << endl;
  // print integers and permute bitmask
  do {
    for (int i = 0; i < p; ++i) // [0..N-1] integers
    {
      if (bitmask[i]){
        myVector.push_back(i);
      }
    }
  } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
  
  imat at = zeros<imat>(p*(p - 1) / 2, 2);
  for (int j = 0; j < p*(p - 1) / 2; j++){
    at(j, 0) = as_scalar(myVector[2 * j]);
    at(j, 1) = as_scalar(myVector[2 * j + 1]);
  }
  
  return at;
}

// [[Rcpp::export]]
vec Mat2Vec(mat R){
  // convering upper triangular matrix into vector 
  imat at = comb(R.n_rows);
  int p = at.n_rows;
  vec RV = zeros(p, 1);
  // cout << at << endl;
  int j, k;
  for(int i = 0; i < (int)p; i++){
    j = at(i, 0);
    k = at(i, 1);
    RV[i] = R(j, k);
  }
  
  return RV;
}

// [[Rcpp::export]]
mat Vec2Mat(vec RV, int p1){
  imat at = comb(p1);
  int p = at.n_rows;
  mat R = eye(p1, p1);
  int j, k;
  for(int i = 0; i < (int)p; i++){
    j = at(i, 0);
    k = at(i, 1);
    R(j, k) = RV[i];
    R(k, j) = RV[i];
  }
  
  return R;
}


// [[Rcpp::export]]
List Cal_block_Rvec(arma::ivec &bp, arma::ivec &chr, arma::uvec &avbIndex,  std::string block_file,
                    std::string stringname3,  int coreNum, double lam){
  // int temp = 1;
  
  vector<DataBlock> datablocks;

  cal_blocks(datablocks, bp, chr, block_file);

  string famfile = stringname3;
  famfile +=".fam";
  string bimfile = stringname3;
  bimfile += ".bim";
  int N = getLineNum(famfile);
  int P = getLineNum(bimfile);
  long long Nlong = (long long)N;
  long long Plong = (long long)P;
  unsigned* X0 = new unsigned[Nlong * Plong];
  readPlink(stringname3, N, P, X0);
  arma::Mat<unsigned>* Xdata =  new arma::Mat<unsigned>(X0, N, P, false, false);
  //
  uword M = Xdata -> n_cols;
  arma::Mat<unsigned>* X = Xdata;
  arma::Mat<unsigned> subX0;
  if(avbIndex.size() != 0 && avbIndex.size() < M){
    subX0 = Xdata->cols(avbIndex);
    X = &subX0;
  }



  uword nblocks = datablocks.size();
  uvec Nb;
  Nb.zeros(nblocks, 1);

  // block_inf: save the start and end information of datablocks
  umat block_inf = zeros<umat>(nblocks, 2);
  for(int ii = 0; ii<(int)(nblocks); ii++){
    block_inf(ii, 0) = datablocks[ii].start;
    block_inf(ii, 1) = datablocks[ii].end;
  }

  field<mat> F4Rblock(nblocks, 1);
  field<ivec> F4index(nblocks, 1);
  uvec Nidex = zeros<uvec>(nblocks, 1);
  // -----------------------------------------------------------------------
  // parallel for correlation matrix
  double ld_r2_thresh = 0.0001; // For independent SNPs.
  // set parallele structure object
  paraBlock_CorR parobj(nblocks, datablocks, X, F4Rblock, F4index,  ld_r2_thresh, Nb, Nidex, lam);
  //set parallel computation
  const int n_thread = coreNum;
  std::vector<std::thread> threads(n_thread);
  for(int i_thread = 0; i_thread < n_thread; i_thread++){
    threads[i_thread] = std::thread(&paraBlock_CorR::update_by_thread_CarCorr, &parobj, i_thread);
  }

  for(int i = 0; i < n_thread; i++){
    threads[i].join();
  }

  F4Rblock = parobj.F4Rblock;
  F4index = parobj.F4index;
  Nidex = parobj.Nidex;
  Nb = parobj.Nb;

  // ---------------------------------------------------------------------------
  // resize the F4Block
  field<vec> F4Rvec(nblocks, 1);
  for(int i = 0; i< (int)(nblocks); i++){
    F4Rvec(i, 0) = Mat2Vec(F4Rblock(i, 0));
  }


  List output = List::create(
    // Rcpp::Named("IndR") = IndR,
    // Rcpp::Named("Indpid") = Indpid,
    Rcpp::Named("F4Rblock") = F4Rblock,
    Rcpp::Named("block_inf") = block_inf,
    Rcpp::Named("F4Rvec") = F4Rvec
    // Rcpp::Named("F4Rblock") = F4Rblock,
    // Rcpp::Named("Nb") = parobj.Nb
    // Rcpp::Named("nblocks") = nblocks
  );
  return output;

}