/*
 * ============================================================================
 *  
 *        Authors: Christina SHI <hshi@cse.cuhk.edu.hk>
 *                 Kevin Yip <kevinyip@cse.cuhk.edu.hk>
 *
 * ============================================================================
 */

#include "true2falseKmer_DP.h"

double true2falseKmer_DP(string errorFile, size_t K){
  vector<double> baseErrors;
  ifstream fin(errorFile);
  double tmp;
  while(fin>>tmp){
    baseErrors.push_back(tmp);
  }
  size_t seq_len=baseErrors.size();
  
  double trueP;
  vector<double> DP(K+1, 0), new_DP(K+1, 0);
  
  tmp=1;
  for(int x=0; x<K; x++){
    tmp*=(1-baseErrors[x]);
  }
  DP[0]=tmp;
  for(int x=1; x<=K; x++){
    tmp=baseErrors[x-1];
    for(int y=x; y<K; y++){
      tmp*=(1-baseErrors[y]);
    }
    DP[x] = tmp;
  }
  
  trueP=DP[0];
  for(int x=K; x<seq_len; x++){
    new_DP[0]=DP[0]*(1-baseErrors[x]);
    for(int y=1; y<=K; y++){
      new_DP[y-1]+=DP[y]*(1-baseErrors[x]);
    }
    new_DP[K]=baseErrors[x];
    
    trueP += new_DP[0];
    DP = new_DP;
    new_DP.assign(K+1, 0);
  }

  return trueP/(seq_len-K+1-trueP);
}
