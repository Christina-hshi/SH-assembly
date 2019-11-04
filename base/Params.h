#pragma once

struct Params{
  size_t K; //k-mer size
  size_t readLen_mean;
  size_t readLen_max;
  size_t thread_num;
  
  //contig related parameters 
  int CONTIG_min_cov;//minimum coverage for a K-mer to use
  int CONTIG_min_SR;//minimum number of supportive reads for two unitigs near hub to be paired. and also for reliable connection between two unitigs. 
  int CONTIG_min_len;//minimum length of contigs to output.
  int CONTIG_tip_len;//branches with len no longer than contig_tip_len are removed while there is alternative branch with length larger than contig_tip_len, which is regarded as true path.

  //For solve repeat
  int REPEAT_p;//reliable distance with at least p(ratio) supportive connections  
  
  //scaffold related parameters
  int SCAF_min_MAPQ;//minimum mapping quality for an alignment to use
  //int SCAF_min_MAPLEN;//minimum aligned length of an qualified alignment

  //For output
  bool OUT_repeat;
  int OUT_CONTIG_min_len;
  bool OUT_haploid;//whether output in haploid or diploid mode

  //For recording statistics
  //atomic<int> STATS_repeat_num;
  //atomic<int> STATS_repeat_num_solvedBy_read;
  Params(){}
};

enum SeqLib_type{FR, RF, TT};

struct SeqLib{
  SeqLib_type type;
  uint32_t insert_size;
  string fname;
};
