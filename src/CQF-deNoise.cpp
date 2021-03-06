/*
 * ============================================================================
 *  Filename:  CQF-deNoise.cpp
 *  
 *  Codes by
 *    "CQF-deNoise"
 *        Authors: Christina SHI <hshi@cse.cuhk.edu.hk>
 *                 Kevin Yip <kevinyip@cse.cuhk.edu.hk>
 *        
 * ============================================================================
 */

#define GRAPH_TRAVERSE

#include "cqf/CQF_mt.h"
#include "cqf/true2falseKmer_DP.h"

boost::program_options::variables_map get_opts(int argc, char* argv[]){
  namespace po = boost::program_options;
  po::options_description desc(string(argv[0])+"  <options>\nOptions");
  desc.add_options() 
    ("help,h", "print help messages") 
    (",k", po::value<int>()->required(), "k-mer size") 
    ("trueKmer,n", po::value<uint64_t>()->required(), "number of unique true k-mers")
    (",N", po::value<uint64_t>()->required(), "total number of k-mers")
    ("alpha,e", po::value<double>()->default_value(-1), "average base error rate, when specified, the <errorProfile> is ignored")
    ("errorProfile", po::value<string>()->default_value(""), "error profile in a file, each line with error rate for the corresponding base, e.g. the error rate of the second base is specified in the second line")
    ("fr", po::value<double>()->default_value(0), "tolerable rate of true k-mers being wrongly removed, default: 1/<trueKmer>")
    ("deNoise", po::value<int>()->default_value(-1), "number of rounds of deNoise, when specified, the <fr> is ignored")
    ("endDeNoise", po::bool_switch()->default_value(false), "call deNoise after processing all the k-mers (not counted into the <deNoise>)")
    (",t", po::value<int>()->default_value(16), "number of threads")
    ("format,f", po::value<char>()->required(), "format of the input: g(gzip); b(bzip2); f(plain fastq)")
    ("input,i", po::value<string>()->required(), "a file containing list of input file name(s), should be in the same directory as the fastq file(s)")
    ("output,o", po::value<string>(), "output file name");
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  
  //if (argc==1 || vm.count("help") || (!vm.count("alpha") && !vm.count("errorProfile"))) {
  if (argc==1 || vm.count("help")) {
    cerr << endl << desc << "\n";
    exit(0);
  }
  if (vm["alpha"].as<double>()==-1 && vm["errorProfile"].as<string>()=="") {
    cerr<<endl<<"Please specify either <alpha> or <errorProfile>"<<endl<<endl;
    cerr << desc << "\n";
    exit(0);
  }
  po::notify(vm);

  return vm;
}

int main(int argc, char* argv[]){
  auto opts = get_opts(argc, argv);

  time_t start_time;
  vector<string> files;
  
  string flist = opts["input"].as<string>();
  string file_prefix="";
  auto pos = flist.find_last_of("/\\"); 
  if(pos != string::npos){
    file_prefix = flist.substr(0, pos+1);
  }

  fstream fin;
  fin.open(flist, ios::in);
  string line;
  if(fin.is_open()){
    while(!fin.eof()){
      getline(fin, line);
      if(line.empty()){
        continue;
      }
      files.push_back(file_prefix+line);
    }
    fin.close();
  }else{
    cerr<<"Failed to open file: "<<flist<<endl;
    return 0;
  }

  uint32_t seed = 2038074761;
  uint16_t thread_num = opts["-t"].as<int>();
  int K = opts["-k"].as<int>();
  double alpha = opts["alpha"].as<double>();
  string errorProfile = opts["errorProfile"].as<string>();
  uint64_t n_true_kmers = opts["trueKmer"].as<uint64_t>();
  uint64_t total_kmers = opts["-N"].as<uint64_t>();
  int num_deNoise = opts["deNoise"].as<int>();
  double fr = opts["fr"].as<double>();
  bool end_deNoise = opts["endDeNoise"].as<bool>();
  uint64_t num_slots;
  uint64_t n_distinct_elts_for_DeNoise;
  
  uint64_t num_true_kmers;
  uint64_t num_false_kmers;
  if(alpha==-1){
    double tmp = true2falseKmer_DP(errorProfile, K);
    num_true_kmers = total_kmers * tmp/(1+tmp);
    num_false_kmers = total_kmers - num_true_kmers;
  }else{
    num_true_kmers = total_kmers * pow(1-alpha, K);
    num_false_kmers = total_kmers - num_true_kmers;
  }

  if(num_deNoise < 0 ){
    if(!fr){
      fr = 1.0/n_true_kmers; 
    }
    num_deNoise = mean_CDF2deNoise(num_true_kmers/n_true_kmers, fr);    
  }
  
  int ave_num_slot_for_encode=0;
  uint64_t tmp=num_true_kmers/n_true_kmers+1;
  while(tmp){
    tmp >>= 7;
    ave_num_slot_for_encode++;
  }
  num_slots = n_true_kmers * (ave_num_slot_for_encode + (double)3/2) + num_false_kmers * 10 /((num_deNoise+1) * 9);

  uint64_t qb;
  uint64_t hb;
  qb=1;
  uint64_t base=2;
  while(base < num_slots){
    qb++;
    base <<= 1;
  }
  uint64_t num_deNoise_upper_bound, num_deNoise_lower_bound, num_slots_tmp;
  num_deNoise_upper_bound = num_deNoise_lower_bound = num_deNoise;
 
  num_slots_tmp = num_slots;
  while(num_deNoise && num_slots_tmp < (1ULL<<qb)){//need a little bit more space?
    num_deNoise--;
    num_slots_tmp = n_true_kmers * (ave_num_slot_for_encode + (double)3/2) + num_false_kmers * 10 / ((num_deNoise+1) * 9);
  }
  if(num_slots_tmp >= (1ULL<<qb)){
    num_deNoise++;
  }
  n_distinct_elts_for_DeNoise = n_true_kmers + num_false_kmers/(num_deNoise+1);  
  
  num_deNoise_lower_bound = num_deNoise;
  num_slots_tmp = n_true_kmers * (ave_num_slot_for_encode + (double)3/2);
  if(num_slots_tmp > (1ULL<<(qb-1))){
    num_deNoise_upper_bound = 0;
  }else{
    num_slots_tmp = num_slots;
    //cerr<<"finding upper boud..."<<endl;
    while(num_slots_tmp >= (1ULL<<(qb-1))){
      //cerr<<num_deNoise_upper_bound<<" ";
      num_deNoise_upper_bound++;
      num_slots_tmp = n_true_kmers * (ave_num_slot_for_encode + (double)3/2) + num_false_kmers * 10 / ((num_deNoise_upper_bound+1) * 9);
    }
    //cerr<<endl;
    if(num_slots_tmp < (1ULL<<(qb-1))){
      num_deNoise_upper_bound--;
    }
  }
  
  hb = qb+8;
  
  string output_file;
  if(!opts.count("output")){
    output_file = "k"+to_string(K)+".t"+to_string(thread_num)+".s"+to_string(qb)+".ser";
  }else{
    output_file = opts["output"].as<string>();
  }
  
  FILE_MODE ftype;
  char tmp_str = opts["format"].as<char>();
  if(tmp_str == 'g'){
    ftype = FILE_MODE::GZIP;
  }else if(tmp_str == 'b'){
    ftype = FILE_MODE::BZIP2;
  }else if(tmp_str == 'f'){
    ftype = FILE_MODE::TEXT;
  }else{
    cerr<<"Unrecognized file type "<<tmp_str<<endl;
    cerr<<"run following to get help"<<endl;
    cerr<<"\t"<<argv[0]<<" --help"<<endl;
    return 0;
  }

  cerr<<"CQF-deNoise settings:"<<endl
    <<"qb: "<<qb<<endl
    <<"hb: "<<hb<<endl
    <<"thread_num: "<<thread_num<<endl
    <<"K: "<<K<<endl
    <<"number of true k-mers: "<<n_true_kmers<<endl
    <<"tolerable wrong removal rate: "<<fr<<endl
    <<"number of deNoise rounds: "<<num_deNoise<<endl
    <<"deNoise after processing all k-mers: "<<(end_deNoise?"true":"false")<<endl
    <<"number of unique k-mers triggering deNoise: "<<n_distinct_elts_for_DeNoise<<endl;
    /*
    cerr<<"#true_kmers(unique)/#true_kmers(in total): "<<n_true_kmers<<"/"<<num_true_kmers<<"("<<num_true_kmers/n_true_kmers<<")"<<endl;
    cerr<<"#false_kmers: "<<num_false_kmers<<endl;
    cerr<<"#distinct_kmers_for_deNoise: "<<n_distinct_elts_for_DeNoise<<endl;
    cerr<<"Desired wrong removal rate: "<<fr<<endl;
    cerr<<"#deNoise rounds: "<<num_deNoise<<endl;
    cerr<<"Final wrong removal rate: "<<cdfpoi_positive(num_deNoise, num_true_kmers/n_true_kmers)<<endl;
    */
    cerr<<"#deNoise rounds leading to the same size of CQF: ["<<num_deNoise_lower_bound<<", "<<(num_deNoise_upper_bound==0?"+oo":to_string(num_deNoise_upper_bound))<<"]"<<endl;
    cerr<<"wrong removal rate leading to the same #deNoise rounds: ["<<cdfpoi_positive(num_deNoise_lower_bound, num_true_kmers/n_true_kmers)<<", ";
    if(num_deNoise_upper_bound == 0){
      cerr<<"1";
    }else{
      cerr<<cdfpoi_positive(num_deNoise_upper_bound, num_true_kmers/n_true_kmers);
    }
    cerr<<"]"<<endl<<endl;

  CQF_mt cqf_mt(qb, hb, thread_num, seed);
  
  start_time = time(NULL);
  cerr<<currentDateTime()<<endl;
  cerr<<"Start to build K-mer spectrum..."<<endl;
  cqf_mt.build_KmerSpectrum(files, FILE_TYPE::FASTQ, ftype, K, n_true_kmers, n_distinct_elts_for_DeNoise, num_deNoise, end_deNoise, fr);
  cqf_mt.save(output_file);
  cerr<<"Finished building K-mer spectrum!"<<endl;
  double seconds = difftime(time(NULL), start_time);
  cerr<<"Time for building K-mer spectrum: "<<seconds<<" seconds."<<endl;

  return 0;
}


