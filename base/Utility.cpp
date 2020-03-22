#include "Utility.h"

size_t NGx(vector<size_t> seq_lens, int x, int ref_len){
  std::sort(seq_lens.rbegin(), seq_lens.rend());
  size_t cutoff_len = ref_len * x / 100;
  size_t current_sum, result;
  current_sum = 0;
  for(auto len : seq_lens){
    current_sum += len;
    if(current_sum >= cutoff_len){
      result = len;
      break;
    }
  }
  return result;
}

double median(vector<int>& nums){
  if(nums.size()==0){
    return 0;
  }else if(nums.size()==1){
    return nums[0];
  }
  sort(nums.begin(), nums.end());
  int tmp = nums.size()/2;
  if(nums.size()%2==0){
    return (nums[tmp-1]+nums[tmp])/2.0;  
  }else{
    return nums[tmp];
  }
}

double semi_global_align(string seq_A, string seq_B){
  double match_score = 1;
  double mismatch_score = -1, gap_penalty = -1;

  size_t seq_a_len = seq_A.length();
  size_t seq_b_len = seq_B.length();

  double** dp = new double*[seq_a_len+1];
  for(size_t x = 0; x <= seq_a_len; x++){
    dp[x] = new double[seq_b_len + 1];
    memset(dp[x], 0, sizeof(double)*(seq_b_len + 1));
  }
  
  for(size_t row = seq_a_len; row > 0;){
    row--;
    for(size_t col = seq_b_len; col > 0;){
      col--;
      double score_tmp;
      if(seq_A[row] == seq_B[col]){
        score_tmp = dp[row+1][col+1] + match_score;
      }
      else{
        score_tmp = dp[row+1][col+1] + mismatch_score;
      }
      score_tmp = max(score_tmp, dp[row+1][col] + gap_penalty);
      score_tmp = max(score_tmp, dp[row][col+1] + gap_penalty);
    
      dp[row][col] = score_tmp;
    }
  }

  double result = 0;

  /*
  //output dp table
  for(size_t x = 0; x <= seq_a_len; x++){
    for(size_t y = 0; y <= seq_b_len; y++){
      cout<< dp[x][y]<<" ";
    }
    cout<<endl;
  }
  */

  //result is the maximum value in first row and first column.
  for(size_t x = 0; x <= seq_a_len; x++){
    result = max(result, dp[x][0]);
  }
  for(size_t x = 0; x <= seq_b_len; x++){
    result = max(result, dp[0][x]);
  }
  //delete allocated memory
  for(size_t x = 0; x <= seq_a_len; x++){
    delete []dp[x];
  }
  delete []dp;

  return result;
}

double pearson_corr(matrix_double A, matrix_double B){
  double av_A, av_B;
  av_A = average(A);
  av_B = average(B);

  A -= av_A;
  B -= av_B;
  
  double num_ele = A.size() * A[0].size();
  //sum of square
  double dev_A = inner_product(A, A) / num_ele;
  double dev_B = inner_product(B, B) / num_ele;
  
  double result  = inner_product(A, B) / num_ele / (sqrt(dev_A) * sqrt(dev_B));
  return result;
}

double average(const matrix_double& m){
  double sum = 0;
  double size = 0;;
  if(m.size() == 0){
    return 0;
  }
  for(size_t r = 0; r < m.size(); r++){
    for(size_t c = 0; c < m[r].size(); c++){
      sum += m[r][c];
    }
    size += m[r].size();
  }
  return sum/size;
}
matrix_double& operator-=(matrix_double& m, const double c){
  for(size_t r = 0; r < m.size(); r++){
    for(size_t col = 0; col < m[r].size(); col++){
      m[r][col] -= c;
    }
  }
  return m;
}

//assume A and B are of the same size
double inner_product(const matrix_double& A, const matrix_double& B){
  double result = 0;
  for(size_t r = 0; r < A.size(); r++){
    for(size_t c = 0; c < A[r].size(); c++){
      result += A[r][c] * B[r][c];
    }
  }
  return result;
}

double jaccard_index(sketch_t A, sketch_t B, int sketch_size){
  sketch_t AuB = new u_int[sketch_size];
  memset(AuB, 0, sizeof(u_int) * sketch_size);

  size_t index = 0;
  for(int i = 0, j = 0; i < sketch_size && j < sketch_size && index < sketch_size;){
    if(A[i] == 0 || B[j] == 0){
      if(A[i] != 0){
        while(index < sketch_size && i < sketch_size){
          AuB[index++] = A[i++];
        }
      }else if(B[i] != 0){
        while(index < sketch_size && j < sketch_size){
          AuB[index++] = B[j++];
        }
      }
      break;
    }
    
    if(A[i] < B[j]){
      AuB[index++] = A[i++];
    }else{
      AuB[index++] = B[j++];
    }
  }
  
  int intersec = 0, uni = index;
  for(size_t x = 0, i = 0, j = 0; x < index; x++){
    while(i < sketch_size && A[i] < AuB[x]){i++;}
    while(j < sketch_size && B[j] < AuB[x]){j++;}
    if(i == sketch_size || j == sketch_size){
      break;
    }
    if(A[i] == AuB[x] && B[j] == AuB[x]){
      intersec++;
    }
  }

 /*
  for(size_t i = 0, j = 0; i < sketch_size && j < sketch_size;){
    if(A[i] == 0 || B[j] == 0){
      while(i < sketch_size && A[i] != 0){
        i++;
        uni++;
      }
      while(j < sketch_size && B[j] != 0){
        j++;
        uni++;
      }
      break;
    }
    while(j < sketch_size && B[j] < A[i]){
      uni++;
      j++;
    }
    if(A[i] == B[j]){
      i++; j++;
      uni++;
    }else{
      i++;
      uni++;
    }
  }
  */

  delete [] AuB;

  return (double)intersec / (double)uni;
}

double jaccard_index(set<size_t> s1, set<size_t> s2){
  size_t insec, unin;
  insec=unin=0;
  set<size_t>::iterator it1 = s1.begin();
  set<size_t>::iterator it2 = s2.begin();
  
  while(it1 != s1.end() && it2 != s2.end()){
    while(*it1 < *it2 && it1 != s1.end()){
      it1++;
      //unin++;
    }
    if(*it1 == *it2 && it1 != s1.end()){
      //unin++;
      insec++;
      it1++; it2++;
    }
    while(*it1 > *it2 && it2 != s2.end()){
      it2++;
      //unin++;
    }
  }

  /*
  while(it1 != s1.end()){
    it1++;
    unin++;
  }
  while(it2 != s2.end()){
    it2++;
    unin++;
  }
  */
  unin = min(s1.size(), s2.size());
  return (double)insec/unin;
}
