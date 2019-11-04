#pragma once

#include <iostream>
#include <algorithm>
#include <cstring>
#include <vector>
#include <set>
#include <unordered_set>
#include <bitset>
#include <cassert>
#include <fstream>

#include <boost/thread/thread.hpp>
#include <boost/lockfree/queue.hpp>
#include <boost/lockfree/spsc_queue.hpp>
#include <boost/atomic.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/poisson.hpp>

#include <time.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/mman.h>

#include "clipp.h"
#include "gqf.h"
#include "base/multithread_io.h"
#include "base/nthash.hpp"
//#include "base/cdflib.hpp"
#include "chunk.h"
//#include "kmer.h"
#include "reader.h"
//#include "hashutil.h"
#include <mutex>

#define BITMASK(nbits) ((nbits) == 64 ? 0xffffffffffffffff : (1ULL << (nbits)) \
												- 1ULL)
#define QBITS_LOCAL_QF 16
#define SPARE_EMPTY_LOCAL_QFS 16

//#define VARIANCE(mean) (1.05*mean-0.371)
//#define CDF(mean, var, x) (0.5*erfc(-(x-mean)/sqrt(2*var)))
//cumulative density function
//#define CDF(mean, x) (0.5*erfc(-(x-mean)/sqrt(2*VARIANCE(mean))))

//probability density function for normal
/*
double PDF(double x, double m)
{
  double s = sqrt(VARIANCE(m));
  static const double inv_sqrt_2pi = 0.3989422804014327;
  double a = (x - m) / s;
  return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
}
*/

using namespace std;
using namespace kmercounting;

enum CQF_MODE{EMPTY, FILEMAP, MEMORY};
enum FILE_TYPE{FASTA, FASTQ};
enum FILE_MODE{TEXT, GZIP, BZIP2};  

bool getFileReader(FILE_MODE fmode, const char* seq_file, reader* file_reader);
/*
double CDF2mean(double cdf_desired){
  //a linear search
  double mean=7;
  double cdf,
  CDF(mean, 1);
  if(cdf < cdf_desired){
    mean--;
    cdf = CDF(mean, 1);
    while(cdf < cdf_desired && mean > 0){
      mean--;
      cdf = CDF(mean, 1); 
    }
    mean += 0.5;
  }else{
    mean++;
    cdf = CDF(mean, 1);
    while(cdf > cdf_desired){
      mean++;
      cdf = CDF(mean, 1);
    }
    mean -= 0.5;
  }
  return mean;
}
*/
inline double pdfpoi(double x, double mean){
  boost::math::poisson_distribution<> p(mean);
  return boost::math::pdf(p, x);
}

inline double pdfpoi_positive(double x, double mean){
  boost::math::poisson_distribution<> p(mean);
  return (boost::math::pdf(p, x)-boost::math::pdf(p, 0))/(1-boost::math::pdf(p, 0));
}

inline double cdfpoi(double x, double mean){
  //double cdf, ccdf;
  //cumpoi(&x, &mean, &cdf, &ccdf);
  boost::math::poisson_distribution<> p(mean);
  return boost::math::cdf(p, x);
}

inline double cdfpoi_positive(double x, double mean){
  boost::math::poisson_distribution<> p(mean);
  return (boost::math::cdf(p, x) - boost::math::cdf(p, 0))/(1-boost::math::cdf(p,0));
}

inline double cdfpoi_positive(const boost::math::poisson_distribution<>& p, double x){
  double tmp = boost::math::cdf(p, 0);
  return (boost::math::cdf(p, x) - tmp)/(1-tmp);
}

int mean_CDF2deNoise(double mean, double cdf_desired){
  int start, end, mid;
  start = 0, end = mean+1;
  boost::math::poisson_distribution<> p(mean);
  double cdf0 = boost::math::cdf(p, 0);
  auto cdf_positive = [&](double x){
    return (boost::math::cdf(p, x)-cdf0)/(1-cdf0); 
  };
  while(cdf_positive(end) < cdf_desired){
    end *= 2;
  }
  //binary search
  while(start<=end){
    if(start==end){
      return start;
    }else if((start+1)==end){
      double tmp1, tmp2;
      tmp1 = cdf_positive(start);
      tmp2 = cdf_positive(end);
      if(tmp2 <= cdf_desired){
        return end;
      }else if(tmp1 <= cdf_desired){
        return start;
      }else{
        return max(start-1, 0);
      }
    }
    mid = (start+end)/2;
    double cdf = cdf_positive(mid);
    if(cdf < cdf_desired){
      start = mid+1;  
    }else if(cdf > cdf_desired){
      end = mid-1;
    }else{
      return start;
    }
  }
  //cerr<<"Problem to decide number of times deNoise being called with mean:"<<mean<<" cdf_desired:"<<cdf_desired<<endl;
  return start;
}

bool nextSmaller(vector<int>& bits){
  int len = bits.size();
  if(bits[len-1]){
    int last_empty=len-2;
    while(last_empty >= 0 && bits[last_empty]){
      last_empty--;
    }
    if(last_empty<0){
      return false;
    }
    int first_nonempty_before_last_empty = last_empty-1;
    while(first_nonempty_before_last_empty >=0 && !bits[first_nonempty_before_last_empty]){
      first_nonempty_before_last_empty--;
    }
    if(first_nonempty_before_last_empty < 0){
      return false;
    }
    bits[first_nonempty_before_last_empty]=0;
    bits[first_nonempty_before_last_empty+1]=1;
    if(first_nonempty_before_last_empty == last_empty-1){
      return true;
    }
    int tmp_len = len-last_empty-1;
    int x = first_nonempty_before_last_empty+2;
    for(int y=0; y < tmp_len; x++, y++){
      bits[x] = 1;
    }
    for(; x < len; x++){
      bits[x] = 0;
    }
    return true;
  }else{
    int last_nonempty=len-2;
    while(last_nonempty >=0 && !(bits[last_nonempty])){
      last_nonempty--;
    }
    if(last_nonempty<0){
    return false;
    }
    bits[last_nonempty]=0;
    bits[last_nonempty+1]=1;
    return true;
  }
}

void PrintProbOfFalseRemoval(double mean, const vector<uint64_t>& blocks){
  size_t block_num = blocks.size();
  double pro_blocks[block_num];
  uint64_t kmer_total=0;
  
  cerr<<"mean: "<<mean<<endl;
  cerr<<block_num<<" blocks:"<<endl<<"\t";
  for(auto ele:blocks)
    cerr<<ele<<"\t";
  cerr<<endl;

  for(auto ele:blocks){
    kmer_total += ele;
  }
  for(int x = 0; x < block_num; x++){
    pro_blocks[x] = (double)blocks[x]/kmer_total;
  }
  double permut = 1;
  double tmp, tmp1, p_falseRemoval;
  p_falseRemoval = 0;
  cerr<<"*Probability of removing true k-mers with frequency from [1, "<<block_num<<"]:"<<endl
    <<std::right<<std::setw(12)<<"frequency"
    <<"  p"<<endl;
  for(int x = 1; x <= block_num; x++){
    permut *= x;
    vector<int> bits(block_num, 0);
    for(int y = 0; y < x; y++){
      bits[y]=1;
    }
    tmp1 = 0;
    do{
      tmp = 1;
      for(int y = 0; y <block_num; y++){
        if(bits[y]){
          tmp *= pro_blocks[y];
        }
      }
      tmp1 += tmp;
    }while(nextSmaller(bits));
    tmp1 *= permut;
    p_falseRemoval += tmp1*pdfpoi(x, mean);
    cerr<<std::right<<std::setw(12)<<x<<"  "
      <<tmp1<<endl;
  }
  cerr<<"probability of false removal: "<<endl
    <<"\t"<<p_falseRemoval<<"(removed true k-mers with frequency <="<<block_num<<") among all true k-mers"<<endl;
  double p_falseRemoval_L1 = p_falseRemoval - pdfpoi(1, mean);
  cerr<<"\t"<<p_falseRemoval_L1<<"(removed true k-mers with frequency >1 and <="<<block_num<<") among all true k-mers"<<endl;
  tmp=0;
  for(int x = 2; x<= blocks.size(); x++){
    tmp += pdfpoi(x, mean);
  }
  double p_falseRemoval_L1_within = p_falseRemoval_L1/tmp;
  cerr<<"\t"<<p_falseRemoval_L1_within<<" among true k-mers with frequency >1 and <="<<block_num<<endl; 
  return ;
}

double ProbOfFalseRemoval(double mean_true_kmer, const vector<uint64_t>& blocks){
  size_t block_num = blocks.size();
  double pro_blocks[block_num];
  uint64_t kmer_total=0;
  for(auto ele:blocks){
    kmer_total += ele;
  }
  for(int x = 0; x < block_num; x++){
    pro_blocks[x] = (double)blocks[x]/kmer_total;
  }
  double permut = 1;
  double tmp, tmp1, prob = 0;
  for(int x = 1; x <= block_num; x++){
    permut *= x;
    vector<int> bits(block_num, 0);
    for(int y = 0; y < x; y++){
      bits[y]=1;
    }
    tmp1 = 0;
    do{
      tmp = 1;
      for(int y = 0; y <block_num; y++){
        if(bits[y]){
          tmp *= pro_blocks[y];
        }
      }
      tmp1 += tmp;
    }while(nextSmaller(bits));
    prob += permut*tmp1*pdfpoi(x, mean_true_kmer);
  }
  return prob;
}

enum RunMode{ExtractKmer, DeNoise, Idle, Unknown};
struct CQF_runtime_mt{
  struct Work_Queue{
    atomic<uint64_t> current_index;
    uint64_t max_index;
    mutex mut;
  };
  atomic<uint64_t> ndistinct_elts;
  atomic<uint64_t> nelts;
  //atomic<int> workDone;
  atomic<RunMode> mode;
  uint8_t workerNum;
  //atomic<uint8_t> workerDone;
  atomic<uint8_t> workerReady;
  uint8_t min_C;
  uint32_t min_DeNoise_len;
  uint64_t ndistinct_true_elts;
  Work_Queue work_queue;
  double p_true_kmer_singleton;//when the probability that a true k-mer occurs once less than "p_true_kmer_singleton", then do the DeNoise
  uint64_t ndistinct_elts_for_DeNoise;//DeNosing when there are at least number of distinct elts
  atomic<uint32_t> num_deNoise;
  const bool end_deNoise;
  //atomic<uint8_t> DeNoise_runs;//record how many times the DeNoise procedure execute
  vector<uint64_t> numberOfKmer_runs;//record number of k-mers added after the last DeNoise procedure

  CQF_runtime_mt(uint8_t t, uint32_t min_DN_len, uint64_t ndistinct_true, uint64_t ndistinct_for_DeNoise, uint32_t n_deNoise, bool end_deNois=true, double p_true=0.01, uint8_t c=2): end_deNoise(end_deNois){
    ndistinct_elts = 0;
    nelts = 0;
    mode = RunMode::ExtractKmer;
    workerNum = t;
    workerReady = t;
    min_C = c;
    min_DeNoise_len = min_DN_len;
    ndistinct_true_elts = ndistinct_true;
    p_true_kmer_singleton = p_true;
    ndistinct_elts_for_DeNoise = ndistinct_for_DeNoise;
    //DeNoise_runs=0;
    numberOfKmer_runs.push_back(0);
    num_deNoise = n_deNoise;
  }

  bool needDeNoise(){
    if(ndistinct_elts < ndistinct_elts_for_DeNoise){
      return false;
    }else{
      return true;
    }
  }
};

struct flush_object{
	QF *local_qf;
	QF *main_qf;
	uint32_t count {0};
	uint32_t ksize {28};
  CQF_runtime_mt* runtime;
};

struct file_pointer{
	std::unique_ptr<reader> freader{nullptr};
	char* part{nullptr};
	char* part_buffer{nullptr};
	FILE_MODE fmode{FILE_MODE::TEXT};
	uint64_t size{0};
	uint64_t part_filled{0};
};

struct seqFile_batch{
  vector<string> file_names;
  boost::lockfree::queue<file_pointer*, boost::lockfree::fixed_sized<true> > ip_files;
  boost::atomic<int> num_files;
  FILE_TYPE ftype;  
  FILE_MODE fmode;
  uint32_t OVERHEAD_SIZE;

  seqFile_batch(const vector<string> file_names, FILE_TYPE ft, FILE_MODE fm): ip_files(file_names.size()){
    ftype = ft;
    fmode = fm;
    num_files = 0;
    OVERHEAD_SIZE = 65535;
    this->file_names = file_names;
    
    for(auto fname : file_names){
      auto* fr = new reader;
      if(getFileReader(fmode, fname.c_str(), fr)){
        file_pointer* fp = new file_pointer;
        fp->fmode = fm;
        fp->freader.reset(fr);
        fp->part_buffer = new char[OVERHEAD_SIZE];
        ip_files.push(fp);
        num_files++;
      }else{
        delete fr;
      }
    }
  }

  ~seqFile_batch(){
    file_pointer* fp;
    while(num_files){
      if(ip_files.pop(fp)){
        if(fp->fmode == FILE_MODE::TEXT){
          fclose(fp->freader->in);
        }else if(fp->fmode == FILE_MODE::GZIP){
          gzclose(fp->freader->in_gzip);
        }else if(fp->fmode == FILE_MODE::BZIP2){
          if(fp->freader->in){
            BZ2_bzReadClose(&(fp->freader->bzerror), fp->freader->in_bzip2);
            fclose(fp->freader->in);
          }
        }
        delete[] fp->part_buffer;
        delete fp;
        num_files--;
      }
    }
  }
};

class CQF_mt{
public:
  QF* qf;
  CQF_MODE cqf_mode;
  QFi qfi;
  uint64_t qb;
  uint64_t hb;

  uint64_t default_value;

  uint16_t t;//number of threads

  CQF_mt(){}
  CQF_mt(uint64_t qb, uint64_t hb, uint16_t t, uint32_t seed=2038074761){
    this->qb = qb;
    this->hb = hb;
    this->t  = t;
    uint64_t nslots;
    if(qb>64){
      throw std::invalid_argument("Only support with less than 64 bits.");
    }else{
      nslots = qb==64?std::numeric_limits<uint64_t>::max():(1ULL << qb);
    }
    qf = new QF();
    if(qb==hb){
      qf_init(qf, nslots, this->hb+8, 0, true, "", seed);
      default_value = 0;//0xFF;
    }else{
      qf_init(qf, nslots, this->hb, 0, true, "", seed);
      default_value = 0; 
    }
  }

  //Basic functions
  void reset(){
    qf_reset(qf);
  }
  
  bool insert(uint64_t key, uint64_t count){
    return qf_insert(qf, key, default_value, count, true, true);
  }
  void remove(uint64_t key, uint64_t count){
    qf_remove(qf, key, default_value, count);
  }
  void remove_all(uint64_t key){
    qf_delete_key_value(qf, key, default_value);
  }
  
  uint64_t count(uint64_t key){
    return qf_count_key_value(qf, key, default_value);
  }
 
  bool set_iterator_2pos(uint64_t pos){
    return qf_iterator(qf, &qfi, pos);
  }

  //return true if the end has not been reached
  bool get(uint64_t& key, uint64_t& count){
    uint64_t tmp;
    return qfi_get(&qfi, &key, &tmp, &count)?false:true;
  }
  //return whether the next entry is found or not
  bool next(){
    return qfi_next(&qfi)?false:true;
  }

  void mmap2file(const string filename){
    cqf_mode = FILEMAP;
    qf_read(qf, filename.c_str());
  }

  void munmap2file(){
    return qf_destroy(qf, false);
  }

#ifdef GRAPH_TRAVERSE
  bool is_traveled(uint64_t key){
    uint64_t count;
    return qf_count_key_value_is_traveled(qf, key, default_value, &count);
  }
  void set_traveled(uint64_t key){
    uint64_t count;
    qf_count_key_value_set_traveled(qf, key, default_value, &count);
  }
  bool next_untraveled(){
    return qfi_next_untraveled(&qfi)?false:true;
  }
  //return whether it is traveled before setting
  bool count_key_value_set_traveled(uint64_t key, uint64_t& count){
    return qf_count_key_value_set_traveled(qf, key, default_value, &count);
  }
  bool count_key_value_is_traveled(uint64_t key, uint64_t& count){
    return qf_count_key_value_is_traveled(qf, key, default_value, &count);
  }
#endif

  void load(const string filename){
    qf = new QF();
    qf_deserialize(qf, filename.c_str());
  }

  void save(const string filename){
    qf_serialize(qf, filename.c_str()); 
  }
  //basic function end

  //functions for multi-thread
  void build_KmerSpectrum(const vector<string>& fnames, const FILE_TYPE ftype, const FILE_MODE fmode, int ksize, uint64_t n_distinct_true_elts, uint64_t n_distinct_elts_for_DeNoise, uint32_t n_deNoise, bool end_deNoise, double p_true_kmer_singleton);

  //funcitons for multi-thread end
  void print_metadata(){
    cerr<<"#metadata"<<endl
      <<"size: "<<qf->metadata->size<<endl
      <<"seed: "<<qf->metadata->seed<<endl
      <<"nslots: "<<qf->metadata->nslots<<endl
      <<"xnslots: "<<qf->metadata->xnslots<<endl
      <<"key_bits: "<<qf->metadata->key_bits<<endl
      <<"value_bits: "<<qf->metadata->value_bits<<endl
      <<"key_remainder_bits: "<<qf->metadata->key_remainder_bits<<endl
      <<"bits_per_slots: "<<qf->metadata->bits_per_slot<<endl
      //<<"range: "<<qf->metadata->range<<endl
      <<"nelts: "<<qf->metadata->nelts<<endl
      <<"ndistinct_elts: "<<qf->metadata->ndistinct_elts<<endl
      <<"noccupied_slots: "<<qf->metadata->noccupied_slots<<endl
      <<"num_locks: "<<qf->metadata->num_locks<<endl;
  }

  ~CQF_mt(){
    if(cqf_mode==CQF_MODE::EMPTY){
      return ;
    }else if(cqf_mode==CQF_MODE::FILEMAP){
      //munmap2file();
      qf_destroy(qf, false);
    }else{
      qf_destroy(qf, true);
      //qf_free(qf);
    }
    delete qf;
    qf=NULL;
    //qfi_free(qfi);
    //qfi=NULL;
  }
};

/* check if it's the end of the file. */
inline bool is_eof(reader &file_reader, FILE_MODE fmode){
	if(fmode == FILE_MODE::TEXT)
		return feof(file_reader.in) != 0;
	else if(fmode == FILE_MODE::GZIP)
		return gzeof(file_reader.in_gzip) != 0;
	else if(fmode == FILE_MODE::BZIP2)
		return file_reader.bzerror == BZ_STREAM_END;

	return true;
}

/* move the pointer to the end of the next newline. */
bool skip_next_eol(char *part, int64_t &pos, int64_t max_pos){
	int64_t i;
	for(i = pos; i < max_pos-2; ++i)
		if((part[i] == '\n' || part[i] == '\r') && !(part[i+1] == '\n' ||
																								 part[i+1] == '\r'))
			break;

	if(i >= max_pos-2)
		return false;
	pos = i+1;

	return true;
}  

/* dump the contents of a local QF into the main QF */
static void dump_local_qf_to_main(flush_object *obj){
	uint64_t new_elts, total_elts; new_elts = total_elts = 0;
  bool isNew;
  QFi local_cfi;

	if (qf_iterator(obj->local_qf, &local_cfi, 0)) {
		do {
			uint64_t key = 0, value = 0, count = 0;
			qfi_get(&local_cfi, &key, &value, &count);
			qf_insert_advance(obj->main_qf, key, 0, count, true, true, isNew);
      if(isNew){
        new_elts++;
      }
      total_elts += count;
		} while (!qfi_next(&local_cfi));
		qf_reset(obj->local_qf);
	}
  obj->runtime->ndistinct_elts += new_elts;
  obj->runtime->nelts += total_elts;
}

/* convert a chunk of the fastq file into kmers */
void reads_to_kmers(chunk &c, flush_object *obj){
	uint64_t new_elts, total_elts; new_elts=total_elts=0;
  bool isNew;
  auto fs = c.get_reads();
	auto fe = c.get_reads();
	auto end = fs + c.get_size();
	while (fs && fs!=end) {
		fs = static_cast<char*>(memchr(fs, '\n', end-fs)); // ignore the first line
		fs++; // increment the pointer

		fe = static_cast<char*>(memchr(fs, '\n', end-fs)); // read the read
		string read(fs, fe-fs);
		/*cerr << read << endl;*/

    uint64_t hash;
    uint64_t hash_RC;
start_read:
		if (read.length() < obj->ksize) // start with the next read if length is smaller than K
			goto next_read;
		{
      NTPC64(read.c_str(), obj->ksize, hash, hash_RC);
      /*
			 * first try and insert in the main QF.
			 * If lock can't be accuired in the first attempt then
			 * insert the item in the local QF.
			 */
			if(hash < hash_RC){
      if (!qf_insert_advance(obj->main_qf, hash%obj->main_qf->metadata->range, 0, 1,
										 true, false, isNew)) {
				qf_insert(obj->local_qf, hash%obj->local_qf->metadata->range, 0, 1, false, false);
				obj->count++;
				// check of the load factor of the local QF is more than 50%
				if (obj->count > 1ULL<<(QBITS_LOCAL_QF-1)) {
					dump_local_qf_to_main(obj);
					obj->count = 0;
				}
			}else{
        if(isNew){
          new_elts++;
        }else{
          total_elts++;
        }
      }
      }else{
      if (!qf_insert_advance(obj->main_qf, hash_RC%obj->main_qf->metadata->range, 0, 1,
										 true, false, isNew)) {
				qf_insert(obj->local_qf, hash_RC%obj->local_qf->metadata->range, 0, 1,
									false, false);
				obj->count++;
				// check of the load factor of the local QF is more than 50%
				if (obj->count > 1ULL<<(QBITS_LOCAL_QF-1)) {
					dump_local_qf_to_main(obj);
					obj->count = 0;
				}
			}else{
        if(isNew){
          new_elts++;
        }else{
          total_elts++;
        }
      }
      }
			for(uint32_t i=obj->ksize; i<read.length(); i++) { //next kmers
				if(read[i]=='N'){
          read = read.substr(i+1, read.length());
          goto start_read;
        }
 				/*
				 * first try and insert in the main QF.
				 * If lock can't be accuired in the first attempt then
				 * insert the item in the local QF.
				 */
				NTPC64(read[i-obj->ksize], read[i], obj->ksize, hash, hash_RC);
        if(hash < hash_RC){
        if (!qf_insert_advance(obj->main_qf, hash%obj->main_qf->metadata->range, 0, 1, true,
											 false, isNew)) {
					qf_insert(obj->local_qf, hash%obj->local_qf->metadata->range, 0, 1, false,
										false);
					obj->count++;
					// check of the load factor of the local QF is more than 50%
					if (obj->count > 1ULL<<(QBITS_LOCAL_QF-1)) {
						dump_local_qf_to_main(obj);
						obj->count = 0;
					}
				}else{
          if(isNew){
            new_elts++;
          }else{
            total_elts++;
          }
        }
        }else{
        if (!qf_insert_advance(obj->main_qf, hash_RC%obj->main_qf->metadata->range, 0, 1, true, false, isNew)) {
					qf_insert(obj->local_qf, hash_RC%obj->local_qf->metadata->range, 0, 1, false, false);
					obj->count++;
					// check of the load factor of the local QF is more than 50%
					if (obj->count > 1ULL<<(QBITS_LOCAL_QF-1)) {
						dump_local_qf_to_main(obj);
						obj->count = 0;
					}
				}else{
          if(isNew){
            new_elts++;
          }else{
            total_elts++;
          }
        }
        }
			}
		}

next_read:
		fs = ++fe;		// increment the pointer
		fs = static_cast<char*>(memchr(fs, '\n', end-fs)); // ignore one line
		fs++; // increment the pointer
		fs = static_cast<char*>(memchr(fs, '\n', end-fs)); // ignore one more line
		fs++; // increment the pointer
	}
	free(c.get_reads());
  obj->runtime->ndistinct_elts += new_elts;
  obj->runtime->nelts += (new_elts+total_elts);
}


/* read a part of the fastq file. */
static bool fastq_read_parts(FILE_MODE fmode, file_pointer *fp){
	char *& _part = (fp->part);
	uint64_t& _size = fp->size;
	char*& part_buffer = (fp->part_buffer);
	uint64_t& part_filled = fp->part_filled;
	reader& file_reader = *(fp->freader.get());

	uint32_t OVERHEAD_SIZE = 65535;
	uint64_t part_size = 1ULL << 23;
	char *part = (char *)malloc((part_size + OVERHEAD_SIZE)*sizeof(char));
	memcpy(part, part_buffer, part_filled);

	if(is_eof(file_reader, fmode))
		return false;

	uint64_t readed = 0;

	if (fmode == FILE_MODE::TEXT)
		readed = fread(part+part_filled, 1, part_size, file_reader.in);
	else if (fmode == FILE_MODE::GZIP)
		readed = gzread(file_reader.in_gzip, part+part_filled, (int) part_size);
	else if (fmode == FILE_MODE::BZIP2)
		readed = BZ2_bzRead(&file_reader.bzerror, file_reader.in_bzip2,
												part+part_filled, (int) part_size);
	else 
		readed = 0;

	int64_t total_filled = part_filled + readed;
	int64_t i;
	if(part_filled >= OVERHEAD_SIZE)
	{
		cerr << "Error: Wrong input file!\n";
		exit(EXIT_FAILURE);
	}
	if(is_eof(file_reader, fmode))
	{
		_part = part;
		_size = total_filled;
		part = NULL;
		return true;
	}
	// Looking for a FASTQ record at the end of the area
	{
		int64_t line_start[9];
		int32_t j;
		i = total_filled - OVERHEAD_SIZE / 2;
		for(j = 0; j < 9; ++j)
		{
			if(!skip_next_eol(part, i, total_filled))
				break;
			line_start[j] = i;
		}
		_part = part;
		if(j < 9)
			_size = 0;
		else
		{
			int k;
			for(k = 0; k < 4; ++k)
			{
				if(part[line_start[k]+0] == '@' && part[line_start[k+2]+0] == '+')
				{
					if(part[line_start[k+2]+1] == '\n' || part[line_start[k+2]+1] == '\r')
						break;
					if(line_start[k+1]-line_start[k] == line_start[k+3]-line_start[k+2] &&
						 memcmp(part+line_start[k]+1, part+line_start[k+2]+1,
										line_start[k+3]-line_start[k+2]-1) == 0)
						break;
				}
			}
			if(k == 4)
				_size = 0;
			else
				_size = line_start[k];
		}
	}

	copy(_part+_size, _part+total_filled, part_buffer);
	part_filled = total_filled - _size;

	return true;
}

/* read a part of the fastq file, parse it, convert the reads to kmers, and
 * insert them in the CQF
 */
static bool fastq_to_uint64kmers_prod(flush_object* obj, seqFile_batch* seqFiles){
	file_pointer* fp;
  CQF_runtime_mt* runtime = obj->runtime;
  while(true){
    runtime->workerReady--;
    if(runtime->mode == RunMode::ExtractKmer){
      while(seqFiles->num_files){
        if (seqFiles->ip_files.pop(fp)) {
          if (fastq_read_parts(fp->fmode, fp)) {
            seqFiles->ip_files.push(fp);
            chunk c(fp->part, fp->size);
            reads_to_kmers(c, obj);
            if (obj->count) {
              dump_local_qf_to_main(obj);
              obj->count = 0;
            }             
            if(runtime->num_deNoise && runtime->needDeNoise()){//runtime->ndistinct_elts > runtime->ndistinct_true_elts*2){
              break;
            }
          } else {
            /* close the file */
            if (fp->fmode == FILE_MODE::TEXT)
              fclose(fp->freader->in);
            else if (fp->fmode == FILE_MODE::GZIP)
              gzclose(fp->freader->in_gzip);
            else if (fp->fmode == FILE_MODE::BZIP2)
              if (fp->freader->in) {
                BZ2_bzReadClose(&(fp->freader->bzerror), fp->freader->in_bzip2);
                fclose(fp->freader->in);
              }
            delete[] fp->part_buffer;
            delete fp;
            seqFiles->num_files--;
          }
        }
      }
      //if(obj->count){
      //  dump_local_qf_to_main(obj);
      //  obj->count = 0;
      //}
      runtime->work_queue.mut.lock();
      if(runtime->workerReady == runtime->workerNum-1){
        if( (runtime->num_deNoise && runtime->needDeNoise()) || (runtime->end_deNoise && !seqFiles->num_files) ){
          if(runtime->num_deNoise){
            runtime->num_deNoise--;
          }
          //(!seqFiles->num_files || runtime->needDeNoise()){
          //runtime->ndistinct_elts > runtime->ndistinct_true_elts*2){
          //runtime->DeNoise_runs++;
          runtime->numberOfKmer_runs.back() = runtime->nelts - runtime->numberOfKmer_runs.back();
          runtime->work_queue.current_index = find_first_nonempty_slot(obj->main_qf, 0);
          cerr<<currentDateTime()<<endl;
          cerr<<"Ready for DeNoise: ndistinct_elts/total_elts."<<runtime->ndistinct_elts<<"/"<<runtime->nelts<<" ndistinct_true_elts."<<runtime->ndistinct_true_elts<<endl;
          runtime->mode = RunMode::DeNoise;
        }else if(!seqFiles->num_files){
          runtime->mode = RunMode::Idle;
        }
        runtime->workerReady++;
        runtime->work_queue.mut.unlock();
        continue;
      }else{
        runtime->workerReady++;
        runtime->work_queue.mut.unlock();
        while(runtime->mode == RunMode::ExtractKmer){
          //runtime->workerReady < runtime->workerNum){
          boost::this_thread::sleep_for(boost::chrono::milliseconds(200));   
        }
        continue;
      }
    }else if(runtime->mode == RunMode::DeNoise){
      uint64_t start, end;
      runtime->work_queue.mut.lock();
      runtime->work_queue.mut.unlock();
      while(runtime->work_queue.current_index < runtime->work_queue.max_index){
        runtime->work_queue.mut.lock();
        start = runtime->work_queue.current_index;
        if(start < runtime->work_queue.max_index){
          end = (start+runtime->min_DeNoise_len)>(obj->main_qf->metadata->nslots)?(obj->main_qf->metadata->nslots):(start+runtime->min_DeNoise_len);
          end = find_first_empty_slot(obj->main_qf, end)-1;
          runtime->work_queue.current_index = find_first_nonempty_slot(obj->main_qf, end+1);
          runtime->work_queue.mut.unlock();
          qf_clean_singleton_with_lock(obj->main_qf, start, end, runtime);
        }else{
          runtime->work_queue.mut.unlock();
          break;
        }
      }
      runtime->work_queue.mut.lock();
      if(runtime->workerReady == runtime->workerNum -1){
        if(seqFiles->num_files){
          runtime->mode = RunMode::ExtractKmer;
          runtime->numberOfKmer_runs.push_back(runtime->nelts);
        }else{
          runtime->mode = RunMode::Idle;
        }
        //runtime->mode = seqFiles->num_files?RunMode::ExtractKmer : RunMode::Idle;
        runtime->workerReady++;
        runtime->work_queue.mut.unlock();
        cerr<<"Finished DeNoise: ndistinct_elts/total_elts."<<runtime->ndistinct_elts<<"/"<<runtime->nelts<<" ndistinct_true_elts."<<runtime->ndistinct_true_elts<<endl;
        cerr<<currentDateTime()<<endl;
        continue;
      }else{
        runtime->workerReady++;
        runtime->work_queue.mut.unlock();
        while(runtime->mode == RunMode::DeNoise){
          //(runtime->workerReady < runtime->workerNum){
          boost::this_thread::sleep_for(boost::chrono::milliseconds(200));   
        }
        continue;
      }
    }else if(runtime->mode == RunMode::Idle){
      free(obj);
      return true;
    }else{
      free(obj);
      return false;
    }
  }
}

bool getFileReader(FILE_MODE fmode, const char* seq_file, reader* file_reader){
  uint64_t gzip_buffer_size = 1ULL << 26;
  uint64_t bzip2_buffer_size = 1ULL << 26;

  if(fmode == FILE_MODE::TEXT){
    if ((file_reader->in = fopen(seq_file, "rb")) == NULL)
      return false;
  }else if(fmode == FILE_MODE::GZIP) {
    if ((file_reader->in_gzip = gzopen(seq_file, "rb")) == NULL)
      return false;
    gzbuffer(file_reader->in_gzip, gzip_buffer_size);
  }else if(fmode == FILE_MODE::BZIP2) {
    file_reader->in = fopen(seq_file, "rb");
    if (!file_reader->in)
      return false;
    setvbuf(file_reader->in, NULL, _IOFBF, bzip2_buffer_size);
    if((file_reader->in_bzip2 = BZ2_bzReadOpen(&file_reader->bzerror,
            file_reader->in, 0, 0, NULL,
            0)) == NULL) {
      fclose(file_reader->in);
      return false;
    }
  }
  return true;
}

void CQF_mt::build_KmerSpectrum(const vector<string>& fnames, const FILE_TYPE ftype, const FILE_MODE fmode, int ksize, uint64_t n_distinct_true_elts, uint64_t n_distinct_elts_for_DeNoise, uint32_t n_deNoise, bool end_deNoise, double p_true_kmer_singleton){
  time_t start_time = time(NULL);
  seqFile_batch* seqFiles = new seqFile_batch(fnames, ftype, fmode);
  QF local_qfs[t];
 
  CQF_runtime_mt* runtime=new CQF_runtime_mt(t, NUM_SLOTS_TO_LOCK<<4, n_distinct_true_elts, n_distinct_elts_for_DeNoise, n_deNoise, end_deNoise, p_true_kmer_singleton);
  runtime->work_queue.max_index = qf->metadata->nslots;
  runtime->work_queue.current_index = 0;

  boost::thread_group prod_threads;
  for(int x = 0; x< t; x++){
    qf_init(&local_qfs[x], (1ULL << QBITS_LOCAL_QF), hb, 0, true, "", qf->metadata->seed);
    flush_object* obj = (flush_object*)malloc(sizeof(flush_object));
    obj->local_qf = &local_qfs[x];
    obj->main_qf = qf;
    obj->ksize = ksize;
    obj->runtime = runtime;
    prod_threads.add_thread(new boost::thread(fastq_to_uint64kmers_prod, obj, seqFiles));
  }
  
  prod_threads.join_all();

  for(int x = 0; x<t; x++){
    qf_destroy(&local_qfs[x], true);
  }
  delete seqFiles;
  cerr<<"Time for building K-mer spectrum without dumping to disk: "<<difftime(time(NULL), start_time)<<" seconds."<<endl;
  qf->metadata->nelts = runtime->nelts;
  qf->metadata->ndistinct_elts = runtime->ndistinct_elts;
  if(!end_deNoise){
    runtime->numberOfKmer_runs.pop_back();
  }
  cerr<<"Estimated probability of true k-mers with freq<="
    <<(end_deNoise?n_deNoise+1:n_deNoise)<<" is: "
    <<cdfpoi((end_deNoise?n_deNoise+1:n_deNoise), (qf->metadata->nelts/(double)qf->metadata->ndistinct_elts))<<endl;
  delete runtime;
}

/* Clean up singletons inside buckets [start_bucket_id, end_bucket_id].
 */
void qf_clean_singleton_with_lock(const QF* qf, uint64_t start_bucket_id, uint64_t end_bucket_id, CQF_runtime_mt* runtime){
	uint64_t start, end, start_block_idx, end_block_idx, removed_elts;
	start_block_idx = start_bucket_id/SLOTS_PER_BLOCK;
	end_block_idx = end_bucket_id/SLOTS_PER_BLOCK;
  removed_elts = 0;
  
  if(start_block_idx==end_block_idx){ 
		qf_spin_lock(&qf->mem->locks[start_bucket_id/NUM_SLOTS_TO_LOCK], true);
		qf_clean_singleton_discrete(qf, start_bucket_id, end_bucket_id, &removed_elts);
    qf_spin_unlock(&qf->mem->locks[start_bucket_id/NUM_SLOTS_TO_LOCK]);
  }else{
    start = start_bucket_id;
    end = find_first_empty_slot(qf, start+1)-1;
		qf_spin_lock(&qf->mem->locks[start_bucket_id/NUM_SLOTS_TO_LOCK], true);
    while((end/SLOTS_PER_BLOCK) == start_block_idx){
      qf_clean_singleton(qf, start, end, &removed_elts);
      start = find_first_nonempty_slot(qf, end+1);
      end = find_first_empty_slot(qf, start+1)-1;
    }
    qf_spin_unlock(&qf->mem->locks[start_bucket_id/NUM_SLOTS_TO_LOCK]);
    if(start/SLOTS_PER_BLOCK == start_block_idx){
      qf_clean_singleton_with_lock_atStart(qf, start, end, &removed_elts);
      start = find_first_nonempty_slot(qf, end+1);
      end = find_first_empty_slot(qf, start+1)-1;
    }
    while(start < end_bucket_id){
      if(end/SLOTS_PER_BLOCK == end_block_idx){
				qf_spin_lock(&qf->mem->locks[end_bucket_id/NUM_SLOTS_TO_LOCK], true);
        qf_clean_singleton_discrete(qf, start, end_bucket_id, &removed_elts);
        qf_spin_unlock(&qf->mem->locks[end_bucket_id/NUM_SLOTS_TO_LOCK]);
				break;
			}else{
				qf_clean_singleton(qf, start, end, &removed_elts);
			}
      start = find_first_nonempty_slot(qf, end+1);
      end = find_first_empty_slot(qf, start+1)-1;
    }
  }
  runtime->nelts -= removed_elts;
  runtime->ndistinct_elts -= removed_elts;
}
