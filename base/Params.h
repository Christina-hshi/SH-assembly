/*
 * ============================================================================
 *  
 *        Authors: Christina SHI <hshi@cse.cuhk.edu.hk>
 *                 Kevin Yip <kevinyip@cse.cuhk.edu.hk>
 *
 * ============================================================================
 */

#pragma once

#include "base/Utility.h"

struct Params{
  //general parameters
  size_t thread_num;

  //k-mer related parameters
  size_t K; //k-mer size
  int kmer_abundance_min; //k-mer with lower abundance won't be used 
  int solid_kmer_abundance_min;
  int solid_kmer_abundance_max;

  //for extending unitigs
  //int dis_deviation_max;//for contigs suggested to be connected by aligning reads (distance: x), they are not connected by the DBG iff we don't see the connection by searching for nodes within (x + max_dis_deviation) bp. 

  //input and output
  string readFileList;
  FILE_TYPE ftype;
  FILE_MODE fmode;
  string cqfFile;
  string output;

  // //read related 
  // size_t readLen_mean;
  // size_t readLen_max;

  // //contig related parameters 
  // int CONTIG_min_cov;//minimum coverage for a K-mer to use
  // int CONTIG_min_SR;//minimum number of supportive reads for two unitigs near hub to be paired. and also for reliable connection between two unitigs. 
  // int CONTIG_min_len;//minimum length of contigs to output.
  // int CONTIG_tip_len;//branches with len no longer than contig_tip_len are removed while there is alternative branch with length larger than contig_tip_len, which is regarded as true path.

  // //For solve repeat
  // int REPEAT_p;//reliable distance with at least p(ratio) supportive connections  
  
  // //scaffold related parameters
  // int SCAF_min_MAPQ;//minimum mapping quality for an alignment to use
  // //int SCAF_min_MAPLEN;//minimum aligned length of an qualified alignment

  // //For output
  // bool OUT_repeat;
  // int OUT_CONTIG_min_len;
  bool OUT_haploid;//whether output in haploid or diploid mode

  //For recording statistics
  //atomic<int> STATS_repeat_num;
  //atomic<int> STATS_repeat_num_solvedBy_read;
  Params(){}
  friend ostream& operator<<( ostream &output, const Params& params);
};

ostream& operator<<( ostream &output, const Params& params){
  output<<"Params:"<<endl
  <<"\tkmer size:                "<<params.K<<endl
  <<"\tkmer min. abundance:      "<<params.kmer_abundance_min<<endl
  <<"\tsolid kmer min. abundance:"<<params.solid_kmer_abundance_min<<endl
  <<"\tsolid kmer max. abundance:"<<params.solid_kmer_abundance_max<<endl
  <<"\tthreads:                  "<<params.thread_num<<endl;
  return output; 
}

enum SeqLib_type{FR, RF, TT};

struct SeqLib{
  SeqLib_type type;
  uint32_t insert_size;
  string fname;
};
