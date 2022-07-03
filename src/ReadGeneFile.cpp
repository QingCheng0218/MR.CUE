#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <thread>
#include <mutex>
#include "ReadGeneFile.hpp"
#include "CalCorr.hpp"
#include "data_loader.hpp"
#include <ctime>
#include "truncatedNormal.hpp"

using namespace Rcpp;
using namespace arma;
using namespace std;


#define MAX_LEN 20
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
int getLineNum(std::string filename){
  
  FILE *pf = fopen(filename.c_str(), "r"); // 鎵撳紑鏂囦欢
  char buf[10000];
  int lineCnt = 0;
  if (!pf) // 鍒ゆ柇鏄惁鎵撳紑鎴愬姛
    return -1;
  while (fgets(buf, 10000, pf)) // fgets寰幆璇诲彇锛岀洿鍒版枃浠舵渶鍚庯紝鎵嶄細杩斿洖NULL
    lineCnt++; // 绱琛屾暟
  fclose(pf);
  return lineCnt;
}

// void getFourGentype(int* geno, std::bitset<8> bits){
//   int idx = 0;
//   for (int j = 0; j < 8; j = j + 2) {
//     if (bits[j] && bits[j + 1]){
//       geno[idx] = 0;
//     }
//     else if (!bits[j] && !bits[j + 1]){
//       geno[idx] = 2;
//     }
//     else if (!bits[j] && bits[j + 1]){
//       geno[idx] = 1;
//     }
//     else if (bits[j] && !bits[j + 1]){
//       geno[idx] = 3;
//     }
//     idx++;
//   }
// }

void getFourGentype(unsigned* geno, std::bitset<8> bits){
  int idx = 0;
  for (int j=0; j < 8; j = j + 2) {
    if(bits[j] && bits[j+1]){
      geno[idx] = 0;
    }else if(!bits[j] && !bits[j+1]){
      geno[idx] = 2;
    }else if(!bits[j] && bits[j+1]){
      geno[idx] = 1;
    }else if(bits[j] && !bits[j+1]){
      geno[idx] = 3;
    }
    idx++;
  }
}



void vec_sub(arma::vec& v0, double m_value, double squaresum){
  v0 -= m_value;
  if (squaresum > 0)
    v0 /= squaresum;
}

// void readPlink(std::string stringname, int N, int P, unsigned* X){
//   
//   // string stringname = dir + dataname;
//   FILE *fp;
//   unsigned char buff[3];
//   string bedfile = stringname + ".bed";
//   fp = fopen(bedfile.c_str(), "rb");
//   if (!fp) return;
//   fread(buff, sizeof(char), 3, fp);
//   
//   std::bitset<8> magic1(buff[0]);
//   std::bitset<8> magic2(buff[1]);
//   std::bitset<8> mode0(buff[2]);
//   
//   if (magic1.to_ulong() != 108 || magic2.to_ulong() != 27){
//     //   cout <<"Error Identifier of plink binary file" << endl;
//   }
//   
//   unsigned long mode = mode0.to_ulong();
//   if (mode == 0){
//     printf("individual-Major Order:improper type of plink file");
//     exit(EXIT_FAILURE);
//   }
//   
//   int n = 0;
//   long charNum = ceil(N*1.0 / 4) * 10000;
//   int leftGenoNum = ceil(N*1.0 / 4)*P;
//   int nblock = ceil(N*1.0 / 4);
//   int nSNP = 0;
//   while (!feof(fp)) {
//     if (leftGenoNum <= 0)
//       break;
//     if (leftGenoNum <= charNum){
//       charNum = leftGenoNum;
//     }
//     char* genotype = new char[charNum];
//     fread(genotype, sizeof(char), charNum, fp);
//     int* geno = new int[4];
//     int nSNPc = int(charNum / nblock); //number of SNPs of this iteration
//     int idx = 0;
//     for (int i = 0; i < nSNPc; i++) {
//       for (int j = 0; j < nblock - 1; j++){
//         std::bitset<8> bits(genotype[idx]);
//         getFourGentype(geno, bits);
//         memcpy(X + nSNP * N + j * 4, geno, 4 * sizeof(int));
//         idx++;
//         leftGenoNum -= 1;
//       }
//       int left = N - (nblock - 1) * 4;
//       std::bitset<8> bits(genotype[idx]);
//       getFourGentype(geno, bits);
//       memcpy(X + nSNP * N + (nblock - 1) * 4, geno, left*sizeof(int));
//       idx++;
//       leftGenoNum -= 1;
//       nSNP++;
//     }
//     delete[] geno;
//     delete[] genotype;
//     n++;
//     //    cout <<n << " processing"<<endl;
//   }
//   
//   
// }

void readPlink(string stringname,int N, int P, unsigned* X){
  
  // string stringname = dir + dataname;
  FILE *fp;
  unsigned char buff[3];
  string bedfile = stringname +".bed";
  fp = fopen(bedfile.c_str(), "rb");
  if (!fp) return;
  fread(buff, sizeof(char), 3, fp);
  
  std::bitset<8> magic1(buff[0]);
  std::bitset<8> magic2(buff[1]);
  std::bitset<8> mode0(buff[2]);
  
  if(magic1.to_ulong() != 108 || magic2.to_ulong() != 27){
    //   cout <<"Error Identifier of plink binary file" << endl;
  }
  
  unsigned long mode =  mode0.to_ulong();
  if(mode == 0){
    printf ("individual-Major Order:improper type of plink file");
    exit (EXIT_FAILURE);
  }
  //     cout << "SNP-Major Order" << endl;
  // }else if(mode == 0){
  //    cout << "individual-Major Order" << endl;
  // }
  // X = new int[N*P];
  //cout << N << "x" << P << endl;
  //cout << sizeof(*X) << endl;
  long n = 0;
  long long charNum = ceil(N*1.0/4)*10000;
  long long leftGenoNum = ceil(N*1.0/4)*P;
  long nblock = ceil(N*1.0/4);
  long nSNP = 0;
  //cout << "nblock: " << nblock << endl;
  //cout << "leftGenoNum: " << leftGenoNum << endl;
  while (!feof(fp)) {
    if(leftGenoNum <= 0)
      break;
    if(leftGenoNum <= charNum){
      charNum  = leftGenoNum;
    }
    char* genotype = new char[charNum];
    fread(genotype, sizeof(char), charNum, fp);
    unsigned* geno = new unsigned[4];
    long nSNPc = long(charNum / nblock); //number of SNPs of this iteration
    
    //cout << n << "-th line: ";
    //cout << "nSNPc: "<< nSNPc << endl;
    // cout << sizeof(int) << endl;
    long long idx = 0;
    for (long i=0; i < nSNPc; i++) {
      
      /*if ( n >= 4 ){
      cout << i << "-th snp" << endl;
      if (i == 0){
      //cout << "break 1 ... " << nSNP << ";" << N << ";" << leftGenoNum << ";" << idx << endl;
      }
    }*/
      for(long j=0; j < nblock - 1; j++){
        long long indx = (long long)(nSNP) * (long long)(N) + (long long)(j*4);
        
        // if ( n == 3 && i == 9999  && j == nblock - 2){
        /*if ( n == 3 && j == nblock - 2 && i > 5000 ){
        // long long indxt = nSNP * N + j*4;
        // cout << "break 2 ... " << endl;
        cout << "break 1 ... "<< n << "-th line, " << i << "-th snp: " << N <<";" << j  << ";" << 
          indx << ";" << leftGenoNum << ";" << idx << endl;
      }*/
        
        std::bitset<8> bits(genotype[idx]);
        getFourGentype(geno,bits);
        memcpy(X + indx, geno, 4*sizeof(unsigned));
        
        
        idx++;
        leftGenoNum -= 1;
  }
      //if ( n == 4  && i == 0 ){
      //	cout << "break 4 ... " << endl;
      //}
      // cout << "break 1 ... " << endl;
      long left = N - (nblock - 1)*4;
      std::bitset<8> bits(genotype[idx]);
      getFourGentype(geno,bits);
      
      long long indx2 = (long long)(nSNP) * (long long)(N) + (long long)(nblock - 1)*4;
      long long indx3 = left*sizeof(unsigned);
      /*if ( n == 3 && i > 5000 ){
      cout << "break 2 ... " << n << "-th line, " << i << "-th snp: " << nSNP << ";" << 
        N << ";" << nblock << ";" << indx2 << ";" << indx3 << endl;
      }*/
      memcpy(X + indx2, geno, indx3);
      idx++;
      leftGenoNum -= 1;
      nSNP ++;
      }
    
    delete[] geno;
    delete[] genotype;
    n++;
    //    cout <<n << " processing"<<endl;
    }
  
  
}
arma::Col<int> get_interval(arma::Col<int> index, uword i){
  Col<int> u_i;
  if(i == 0){
    u_i = linspace< Col<int> >(0, index(0) - 1, index(0));
  }else{
    uword total = index(i) - index(i-1);
    u_i = linspace< Col<int> >(index(i-1), index(i) - 1, total);
  }
  return u_i;
}