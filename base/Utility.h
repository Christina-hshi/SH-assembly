/* implementing generic functions
*/

#pragma once

#include "global.h"
#include "Hash.h"
#include "DNA_string.h"

using namespace std;

typedef vector<vector<double>> matrix_double;
typedef vector<double> vec_double;

typedef vector<vector<int>> matrix_int;
typedef vector<int> vec_int;


struct Contig
{
  DNAString seq;
  double median_abundance;//median abundances of k-mers in
  Contig(string s="", double a=0){
    seq = s;
    median_abundance = a;
  }
  Contig(DNAString s, double a=0){
    seq = s;
    median_abundance = a;
  }
  void clear(){
    seq.clear();
  }
  bool is_empty() const{
    return (seq.length() == 0);
  }
};

/*Compute median from a vector of int
*/
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

/** Get reverse complement of DNA sequence
*/
inline string RC_DNA(const string seq){
  string result(seq);
  
  map<char, char> nucleo_C{{'A', 'T'}, {'a', 't'}, {'T', 'A'}, {'t', 'a'}, {'G', 'C'}, {'g', 'c'}, {'C', 'G'}, {'c', 'g'}, {'N', 'N'}, {'n', 'n'}};
  for(size_t x = 0, y = seq.length()-1; x < seq.length(); x++, y--){
    result[y] = nucleo_C[seq[x]];
    /*
    switch(seq[x]){
    case 'A':
    //case 'a':
      result[y] = 'T';
      break;
    case 'T':
    //case 't':
      result[y] = 'A';
      break;
    case 'G':
    //case 'g':
      result[y] = 'C';
      break;
    case 'C':
    //case 'c':
      result[y] = 'G';
      break;
    default:
      throw invalid_argument("DNA sequence shouldn't contain " + string(1, seq[x]) + "!");
      break;
    }
    */
  }
  return result;
}

inline char RC_DNAbase(const char base){
  switch(base){
    case 'A':
    case 'a':
      return 'T';
      break;
    case 'T':
    case 't':
      return 'A';
      break;
    case 'G':
    case 'g':
      return 'C';
      break;
    case 'C':
    case 'c':
      return 'G';
      break;
    default:
      throw invalid_argument( base + " is not a valid DNA base.");
      break;
  }
}

inline string& to_upper_DNA(string& seq){
  for(size_t x = 0; x < seq.length(); x++){
    /*
    switch(seq[x]){
      case 'a':
        seq[x] = 'A';
        break;
      case 'c':
        seq[x] = 'C';
        break;
      case 'g':
        seq[x] = 'G';
        break;
      case 't':
        seq[x] = 'T';
        break;
      case 'n':
        seq[x] = 'N';
        break;
      case 'A':
      case 'C':
      case 'G':
      case 'T':
      case 'N':
        break;
      default:
        throw std::invalid_argument("DNA sequence shouldn't contain " + string(1, seq[x]) + "!");
        break;
    }
    */
    if(seq[x]>0x60){
      seq[x] -= 32;
    }
  }
  return seq;
}

inline string& to_lower_DNA(string& seq){
  for(size_t x = 0; x < seq.length(); x++){
    switch(seq[x]){
      case 'A':
        seq[x] = 'a';
        break;
      case 'C':
        seq[x] = 'c';
        break;
      case 'G':
        seq[x] = 'g';
        break;
      case 'T':
        seq[x] = 't';
        break;
      case 'N':
        seq[x] = 'n';
        break;
      case 'a':
      case 'c':
      case 'g':
      case 't':
      case 'n':
        break;
      default:
        throw std::invalid_argument("DNA sequence shouldn't contain " + string(1, seq[x]) + "!");
        break;
    }
  }
  return seq;
}

template<std::size_t N>
struct Bit_comparer{
  bool operator()(const bitset<N>& a, const bitset<N>& b) const{
    for(size_t i = 0; i < N; i++){
      if(a[i]&&!b[i]){ return false;}
      if(!a[i]&&b[i]){ return true;}
    }
    return false;
  }
};

template<typename T>
set<T> intersection(set<T> a, set<T> b){
  set<T> result;
  if(a.size() < b.size()){
    for(auto ele:a){
      if(b.find(ele) != b.end()){
        result.insert(ele);
      }
    }
  }else{
    for(auto ele:b){
      if(a.find(ele) != a.end()){
        result.insert(ele);
      }
    }
  }
  return result;
}

//encode A C G T by 00 01 10 11
//K <= 16
inline vector<u_int> extract_kmer16(const string& seq, const int K){
  vector<u_int> kmers, kmers_RC;
  if(K>16){
    throw std::invalid_argument("K:"+to_string(K)+"is not supported, valid K <= 16.");
  }
 
  size_t seq_len = seq.length();
  u_int kmer_idx = 0; u_int kmer_idx_RC = 0;
  /*find the first k-mer without N
   */
  size_t start_pos = 0;
  while(start_pos < seq_len){
    size_t tmp = seq.find('N', start_pos);
    if(tmp == string::npos){
      if(seq_len-start_pos < K){
        return kmers;
      }else{
        break;
      }
    }else if(tmp - start_pos >= K){
      break;
    }else{
      start_pos = tmp+1;
    }
  }
  for(size_t x = start_pos; x < (start_pos+K-1) && x < seq_len; x++){
    kmer_idx <<= 2;
    kmer_idx_RC >>= 2;
    switch(seq[x]){
      case 'A':
        kmer_idx_RC |= 3221225472;
        break;
      case 'C':
        kmer_idx |= 1;
        kmer_idx_RC |= 2147483648;
        break;
      case 'G':
        kmer_idx |= 2;
        kmer_idx_RC |= 1073741824;
        break;
      case 'T':
        kmer_idx |= 3;
        break;
      default:
        break;
    }
  }

  u_int last_kmer;
  u_int new_start_pos;
  int repeat_time = 0;
  uint64_t kmer_mask = pow(4, K)-1;;
  
  for(int x = start_pos+K-1; x < seq_len; x++){
    kmer_idx <<= 2;
    kmer_idx_RC >>= 2;
    switch(seq[x]){
      case 'A':
        kmer_idx_RC |= 3221225472;
        break;
      case 'C':
        kmer_idx |= 1;
        kmer_idx_RC |= 2147483648;
        break;
      case 'G':
        kmer_idx |= 2;
        kmer_idx_RC |= 1073741824;
        break;
      case 'T':
        kmer_idx |= 3;
        break;
      case 'N':
        new_start_pos = x+1;
        while(new_start_pos<seq_len){
          size_t tmp = seq.find('N', new_start_pos);
          if(tmp == string::npos){
            if(seq_len-new_start_pos < K){
              repeat_time = seq_len -x;
              x = seq_len + 1;//end loop
              break;
            }else{
              repeat_time = new_start_pos + K - 1 - x;
              x = new_start_pos;
              break;
            }
          }else if(tmp-new_start_pos < K){
            new_start_pos = tmp+1;
          }else{
            x = new_start_pos;
            repeat_time = new_start_pos + K - 1 - x;
            break;
          }
        }
        if(new_start_pos == seq_len){
          repeat_time = 1;
          x=seq_len+1;
        }
        last_kmer = kmers.back();
        kmers.insert(kmers.end(), repeat_time, last_kmer);
        last_kmer = kmers_RC.back();
        kmers_RC.insert(kmers_RC.end(), repeat_time, last_kmer);
       
        //process the k bit of k-mer
        kmer_idx = 0;
        kmer_idx_RC = 0;
        for(int t = x+K; x<t && x<seq_len; x++){
          kmer_idx <<= 2;
          kmer_idx_RC >>= 2;
          switch(seq[x]){
            case 'A':
              kmer_idx_RC |= 3221225472;
              break;
            case 'C':
              kmer_idx |= 1;
              kmer_idx_RC |= 2147483648;
              break;
            case 'G':
              kmer_idx |= 2;
              kmer_idx_RC |= 1073741824;
              break;
            case 'T':
              kmer_idx |= 3;
              break;
            default:
              break;
          }
        }
        break;
      default:
        break;
    }
    if(x>seq_len){
      break;
    }
    kmers.push_back(kmer_idx & kmer_mask);
    kmers_RC.push_back(kmer_idx_RC >> 2*(16-K));
  }
  kmers.insert(kmers.end(), kmers_RC.crbegin(), kmers_RC.crend());
  return kmers;
}

//encode A C G T by 00 01 10 11
//K <= 32
inline vector<uint64_t> extract_kmer32(const string& seq, const int K){
  vector<uint64_t> kmers, kmers_RC;
  if(K>32){
    throw std::invalid_argument("K:"+to_string(K)+"is not supported, valid K <= 32.");
  }

  size_t seq_len = seq.length();
  uint64_t kmer_idx = 0; uint64_t kmer_idx_RC = 0;
  /*find the first k-mer without N
   */
  size_t start_pos = 0;
  while(start_pos < seq_len){
    size_t tmp = seq.find('N', start_pos);
    if(tmp == string::npos){
      if(seq_len-start_pos < K){
        return kmers;
      }else{
        break;
      }
    }else if(tmp - start_pos >= K){
      break;
    }else{
      start_pos = tmp+1;
    }
  }
  for(size_t x = start_pos; x < (start_pos+K-1) && x < seq_len; x++){
    kmer_idx <<= 2;
    kmer_idx_RC >>= 2;
    switch(seq[x]){
      case 'A':
        kmer_idx_RC |= 0xC000000000000000;
        break;
      case 'C':
        kmer_idx |= 1;
        kmer_idx_RC |= 0x8000000000000000;
        break;
      case 'G':
        kmer_idx |= 2;
        kmer_idx_RC |= 0x4000000000000000;
        break;
      case 'T':
        kmer_idx |= 3;
        break;
      default:
        break;
    }
  }

  uint64_t last_kmer;
  u_int new_start_pos;
  int repeat_time = 0;
  uint64_t kmer_mask;
  if(K==32){
    kmer_mask = 0xFFFFFFFFFFFFFFFF;
  }else{
    kmer_mask = pow(4, K)-1;
  }
  for(int x = start_pos+K-1; x < seq_len; x++){
    kmer_idx <<= 2;
    kmer_idx_RC >>= 2;
    switch(seq[x]){
      case 'A':
        kmer_idx_RC |= 0xC000000000000000;
        break;
      case 'C':
        kmer_idx |= 1;
        kmer_idx_RC |= 0x8000000000000000;
        break;
      case 'G':
        kmer_idx |= 2;
        kmer_idx_RC |= 0x4000000000000000;
        break;
      case 'T':
        kmer_idx |= 3;
        break;
      case 'N':
        new_start_pos = x+1;
        while(new_start_pos<seq_len){
          size_t tmp = seq.find('N', new_start_pos);
          if(tmp == string::npos){
            if(seq_len-new_start_pos < K){
              repeat_time = seq_len -x;
              x = seq_len + 1;//end loop
              break;
            }else{
              repeat_time = new_start_pos + K - 1 - x;
              x = new_start_pos;
              break;
            }
          }else if(tmp-new_start_pos < K){
            new_start_pos = tmp+1;
          }else{
            x = new_start_pos;
            repeat_time = new_start_pos + K - 1 - x;
            break;
          }
        }
        if(new_start_pos == seq_len){
          repeat_time = 1;
          x=seq_len+1;
        }
        last_kmer = kmers.back();
        kmers.insert(kmers.end(), repeat_time, last_kmer);
        last_kmer = kmers_RC.back();
        kmers_RC.insert(kmers_RC.end(), repeat_time, last_kmer);
       
        //process the k bit of k-mer
        kmer_idx = 0;
        kmer_idx_RC = 0;
        for(int t = x+K; x<t && x<seq_len; x++){
          kmer_idx <<= 2;
          kmer_idx_RC >>= 2;
          switch(seq[x]){
            case 'A':
              kmer_idx_RC |= 0xC000000000000000;
              break;
            case 'C':
              kmer_idx |= 1;
              kmer_idx_RC |= 0x8000000000000000;
              break;
            case 'G':
              kmer_idx |= 2;
              kmer_idx_RC |= 0x4000000000000000;
              break;
            case 'T':
              kmer_idx |= 3;
              break;
            default:
              break;
          }
        }
        break;
      default:
        break;
    }
    if(x>seq_len){
      break;
    }
    kmers.push_back(kmer_idx & kmer_mask);
    kmers_RC.push_back(kmer_idx_RC >> 2*(32-K));
  }
  kmers.insert(kmers.end(), kmers_RC.crbegin(), kmers_RC.crend());
  return kmers;
}

inline vector<uint32_t> extract_kmer_hash32(const string& seq, const int K){
  vector<uint32_t> kmers, kmers_RC;
  int seq_len = seq.length();
  if(seq_len < K){
    return kmers;
  }
  string RC_seq = RC_DNA(seq);
  size_t start_pos = 0;
  while(start_pos <= seq_len-K){
    kmers.push_back(MurmurHash2(seq.substr(start_pos, K)));
    kmers_RC.push_back(MurmurHash2(RC_seq.substr(start_pos, K)));
    start_pos++;
  }
  kmers.insert(kmers.end(), kmers_RC.begin(), kmers_RC.end());
  return kmers;
}

inline array<uint64_t, 4> kmers_before(uint64_t kmer, const int K){
  if(K>32){
    throw std::invalid_argument("K: "+to_string(K)+"  is larger than expected(<=64)");
  }
  array<uint64_t, 4> kmers;
  kmer >>= 2;
  for(uint64_t x= 0; x< 4; x++){
    kmers[x] = kmer|(x<<((K-1)*2));
  }
  return kmers;
}

inline array<string, 4> kmers_before(string kmer){
  array<string, 4> kmers;
  kmer=kmer.substr(0, kmer.length()-1);
  kmers[0] = 'A'+kmer;
  kmers[1] = 'T'+kmer;
  kmers[2] = 'G'+kmer;
  kmers[3] = 'C'+kmer;
  return kmers;
}

inline array<uint64_t, 4> kmers_after(uint64_t kmer, const int K){
  if(K>32){
    throw std::invalid_argument("K: "+to_string(K)+"  is larger than expected(<=64)");
  }
  array<uint64_t, 4> kmers;
  uint64_t kmer_mask;
  if(K==32){
    kmer_mask = 0xFFFFFFFFFFFFFFFF;
  }else{
    kmer_mask = (1LL<<(K*2)) - 1;
  }
  kmer <<= 2;
  for(uint64_t x= 0; x< 4; x++){
    kmers[x] = (kmer|x)&kmer_mask;
  }
  return kmers;
}

inline array<string, 4> kmers_after(string kmer){
  array<string, 4> kmers;
  kmer=kmer.substr(1);
  kmers[0] = kmer+'A';
  kmers[1] = kmer+'T';
  kmers[2] = kmer+'G';
  kmers[3] = kmer+'C';
  return kmers;
}

inline uint64_t RC(uint64_t kmer, int K){
  uint64_t RC_kmer = 0;
  for(int x = 0; x<K; x++){
    switch(kmer&3){
    case 0:
      RC_kmer |= 3;
      break;
    case 1:
      RC_kmer |= 2;
      break;
    case 2:
      RC_kmer |= 1;
      break;
    case 3:
      break;
    }
    if(x!=K-1){
      RC_kmer <<= 2;
      kmer >>= 2;
    }
  }
  return RC_kmer;
}

inline string toDNASeq(uint64_t kmer, int K){
  string seq(K, ' ');
  for(int x = K-1; x>=0; x--){
    switch(kmer&3){
    case 0:
      seq[x] = 'A';
      break;
    case 1:
      seq[x] = 'C';
      break;
    case 2:
      seq[x] = 'G';
      break;
    case 3:
      seq[x] = 'T';
      break;
    }
    kmer >>= 2;
  }
  return seq;
}

/*
inline void display_mallinfo(){
  struct  mallinfo mi;

  mi = mallinfo();

  printf("Total non-mmapped bytes (arena):       %d\n", mi.arena);
  printf("# of free chunks (ordblks):            %d\n", mi.ordblks);
  printf("# of free fastbin blocks (smblks):     %d\n", mi.smblks);
  printf("# of mapped regions (hblks):           %d\n", mi.hblks);
  printf("Bytes in mapped regions (hblkhd):      %d\n", mi.hblkhd);
  printf("Max. total allocated space (usmblks):  %d\n", mi.usmblks);
  printf("Free bytes held in fastbins (fsmblks): %d\n", mi.fsmblks);
  printf("Total allocated space (uordblks):      %d\n", mi.uordblks);
  printf("Total free space (fordblks):           %d\n", mi.fordblks);
  printf("Topmost releasable block (keepcost):   %d\n", mi.keepcost);
}
*/
inline void display_mallinfo(){}

template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {
  // initialize original index locations
  vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

template <typename InputIterator>
inline InputIterator advance(InputIterator& curr, const InputIterator& end, int n){
  size_t remaining(std::distance(curr, end));
  if (remaining < n){
    curr=end;
  }else{
    std::advance(curr, n);
  }
}

template <typename t_name>
inline t_name sum(vector<t_name> nums){
  t_name result = 0;
  for(auto ele:nums){
    result += ele;
  }
  return result;
}

template <typename InputIterator>
inline void shuffle(const InputIterator& start, const InputIterator& end){
  size_t len=std::distance(start, end);
  int swap_time = 3;
  size_t idx1, idx2;
  for(int x = 0; x< swap_time; x++){
    idx1 = rand()%len;
    idx2 = rand()%len;
    if(idx1 == idx2){
      continue;
    }else{
      auto it1 = std::next(start, idx1);
      auto it2 = std::next(start, idx2);
      auto tmp = *it1;
      *it1 = *it2;
      *it2 = tmp;
    }
  }
}

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
inline const std::string DisplayCurrentDateTime() {
  time_t     now = time(0);
  struct tm  tstruct;
  char       buf[80];
  tstruct = *localtime(&now);
  strftime(buf, sizeof(buf), "%Y-%m-%d %X", &tstruct);
  cout<<buf<<endl;
  return buf;
}

inline int shared_kmer_num(vector<u_int> s1, vector<u_int> s2){
  int count = 0;
  for(auto& ele2 : s2){
    for(auto& ele1 : s1){
      if(ele1 == ele2){
        count++;
        break;
      }
    }
  }
  return count;
}

inline int shared_kmer_num(vector<u_int> v1, vector<u_int> s2, int t){
  int count = 0;
  int continuous = 0;
  set<u_int> s1(v1.begin(), v1.end());
  for(auto& ele2 : s2){
    bool flag=false;
    
    if(s1.find(ele2) != s1.end()){
      continuous++;
      flag = true;
    }

    if(!flag){
      if(continuous >= t){
        count += continuous;
      }
      continuous = 0;
    }
  }
  if(continuous >= t){
    count += continuous;
  }
  return count;
}

inline double shared_kmer_ratio2(vector<u_int> v1, vector<u_int> s2, int t){
  double result = 0;
  int continuous = 0;
  
  set<u_int> s1(v1.begin(), v1.end());

  int delim = s2.size();
  for(auto& ele2 : s2){
    bool flag=false;
    
    if(s1.find(ele2) != s1.end()){
      continuous++;
      flag = true;
    }
    
    if(!flag){
      if(continuous >= t){
        result += continuous*continuous;
      }
      continuous = 0;
    }
  }
  if(continuous >= t){
    result += continuous*continuous;
  }
  return sqrt(result/(delim*delim));
}

inline void display_stack_limit(){
  struct rlimit limit;
  getrlimit (RLIMIT_STACK, &limit);
  printf ("\nStack Limit = %ld and %ld max\n", limit.rlim_cur, limit.rlim_max);
}
/*
 * Definition of a tandom repeat
 *  same pieces of DNA repeats several times
 *  1. splits seq into small pieces(kmer) k 
 *  2. count the number of occurance of each kmer
 *  3. sum up the occurance of kmers repeats at least t times
 *  4. if the sum is large than h, then it is a tandom repeat
 */
inline bool isTandomRepeat(string seq){
  int k = 4;
  int t = 10;
  int h = 100;
  vector<u_int> kmers = extract_kmer16(seq, k);
  
  map<u_int, int> kmer_freq;
  for(auto kmer : kmers){
    if(kmer_freq.find(kmer) == kmer_freq.end()){
      kmer_freq[kmer] = 1;
    }else{
      kmer_freq[kmer]++;
    }
  }
  int sum = 0;
  for(auto ele : kmer_freq){
    if(ele.second >= t){
      sum += ele.second;
    }
  }
  return sum>=h?true:false;
}

inline int shared_kmer_num(set<u_int> s1, set<u_int> s2){
  int count = 0;
  for(auto& ele : s2){
    if(s1.find(ele) != s1.end()){
      count++;
    }
  }
  return count;
}

/**Semi-global alignment
* Default parameters are used
* mismatch, gap(opening&extending) penalty: -1
* match score: 1
*/
double semi_global_align(string seq_A, string seq_B);

template<typename T>
double mean(vector<T> vec){
  double sum = 0;
  for(auto ele : vec){
    sum += ele;
  }
  return sum/vec.size();
}

template<typename T>
double variance(vector<T> vec){
  double ave = mean(vec);
  double var = 0;
  for(auto ele : vec){
    var += (ave-ele)*(ave-ele);
  }
  return var;
}

inline double jaccard_coe(u_int a, u_int b){
  u_int join = a & b;
  u_int combine = a | b;

#if defined(COMPILER_GNU) || defined(COMPILER_CLANG) || defined(COMPILE_INTEL)
  return (double)__builtin_popcount(join) / (double)__builtin_popcount(combine);
#else
  cout<<"no built in popcount function"<<endl;
  return 0;
#endif
  
}

/*
  jaccard_index(A, B) = S(AUB)nS(A)nS(B) / S(AUB), where A and B are two sets or arrays, with elements in ascending order. S(AUB) is the smallest top n values for AUB, n is the size of A and B.
  Assume no zero hash value
*/
double jaccard_index(sketch_t A, sketch_t B, int sketch_size);


/*
  jaccard index between two sets
  jaccard_index(set1, set2) = |set1 n set2|/|set1 U set2|
*/
double jaccard_index(set<size_t> s1, set<size_t> s2);

/**take average
*/
double average(const matrix_double& m);
double average(const vec_double& v);

/**overload - operator
*/
inline vec_int& operator+=(vec_int& A, const vec_int& B){
  for(size_t x = 0; x < A.size(); x++){
    A[x] += B[x];
  }
  return A;
}

matrix_double& operator-=(matrix_double& m, const double c);
//matrix_double& operator-=(matrix_double& A, const matrix_double& B);
//matrix_double operator-(const matrix_double& m, const double c);

//vec_double& operator-=(vec_double& v, double c);
//vec_double& operator-=(vec_double& v1, const vec_double& v2);
//vec_double operator-(const vec_double& v, const double c);

/** inner product
*/
double inner_product(const matrix_double& A, const matrix_double& B);

/**calculate pearson correlation
* appliable to both matrix and vector
*/
double pearson_corr(matrix_double A, matrix_double B);
//double pearson_corr(vec_double A, vec_double B);

/*
inline vector<int*>& operator=(vector<int*>& A, const vector<int*>& B){
  int n = 5;
  A.clear();
  for(auto& ele : B){
    int tmp[n];
    std::copy(ele, ele + n, tmp);
    A.push_back(tmp);
  }
  return A;
}
*/

template <typename T>
T convert_to (const std::string &str)
{
  std::istringstream ss(str);
  T num;
  ss >> num;
  return num;
}

inline void contig_summary(vector<string> contigs, int ref_len=4800000){
  size_t contig_num = 0;
  contig_num += contigs.size();
  
  int total_len = 0;
  int NG50, N50; NG50 = N50 = 0;
  vector<int> contig_lens(contig_num);

  int tmp_idx = 0;
  for(auto contig:contigs){
    contig_lens[tmp_idx] = contig.length();
    total_len += contig_lens[tmp_idx];
    tmp_idx++;
  }

  std::sort(contig_lens.rbegin(), contig_lens.rend());
  int sum_tmp = 0;
  for(auto len:contig_lens){
    sum_tmp += len;
    if(sum_tmp>=total_len/2){
      N50 = len;
      break;
    }
  }
  
  sum_tmp = 0;
  for(auto len:contig_lens){
    sum_tmp += len;
    if(sum_tmp>=ref_len/2){
      NG50 = len;
      break;
    }
  }
  cout<<endl<<"Contig statistics: "<<endl;
  cout<<contig_num<<" contigs, "<<total_len<<" bp in total."<<endl;
  cout<<"Min contig len "<<contig_lens.back()<<", max contig len "<<contig_lens.front()<<"."<<endl;
  cout<<"Contig N50 "<<N50<<", contig NG50 "<<NG50<<endl;
  return ;
}

inline void contig_summary(vector<Contig> contigs, int ref_len=4800000){
  size_t contig_num = 0;
  contig_num += contigs.size();
  
  int total_len = 0;
  int NG50, N50; NG50 = N50 = 0;
  vector<int> contig_lens(contig_num);

  int tmp_idx = 0;
  for(auto contig:contigs){
    contig_lens[tmp_idx] = contig.seq.length();
    total_len += contig_lens[tmp_idx];
    tmp_idx++;
  }

  std::sort(contig_lens.rbegin(), contig_lens.rend());
  int sum_tmp = 0;
  for(auto len:contig_lens){
    sum_tmp += len;
    if(sum_tmp>=total_len/2){
      N50 = len;
      break;
    }
  }
  
  sum_tmp = 0;
  for(auto len:contig_lens){
    sum_tmp += len;
    if(sum_tmp>=ref_len/2){
      NG50 = len;
      break;
    }
  }
  cout<<endl<<"Contig statistics: "<<endl;
  cout<<contig_num<<" contigs, "<<total_len<<" bp in total."<<endl;
  cout<<"Min contig len "<<contig_lens.back()<<", max contig len "<<contig_lens.front()<<"."<<endl;
  cout<<"Contig N50 "<<N50<<", contig NG50 "<<NG50<<endl;
  return ;
}
