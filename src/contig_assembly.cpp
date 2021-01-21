/*
 * ============================================================================
 *  Filename:  contig_assembly.cpp
 *  
 *  Codes by
 *    "SH-assembly"
 *        Authors: Christina SHI <hshi@cse.cuhk.edu.hk>
 *                 Kevin Yip <kevinyip@cse.cuhk.edu.hk>
 *        
 * ============================================================================
 */

#define GRAPH_TRAVERSE

#include "cqf/CQF_mt.h"
#include "core/unitig_graph.h"

//typedef tbb::concurrent_unordered_set<DNAString, DNAString_hasher, DNAString_equalto> unordered_set_mt;
//using namespace boost::program_options;
//namespace po = boost::program_options;
//using boost::program_options::variables_map;

Params get_opts(int argc, char* argv[]){
  namespace po = boost::program_options;
  po::options_description desc(string(argv[0])+"  <options>\nOptions:");
  desc.add_options() 
    ("help,h", "print help messages") 
    (",k", po::value<int>()->required(), "k-mer size") 
    ("input,i", po::value<string>()->required(), "a file containing a list of read file name(s), should be absolute address if the read files are not in the running directory")
    ("format,f", po::value<char>()->default_value('f'), "format of the input: g(gzip); b(bzip2); f(plain fastq)")
    ("cqf,c", po::value<string>()->required(), "the counting quotient filter built with the same 'k'")
    ("abundance_min,s", po::value<int>()->default_value(2), "minimum coverage of k-mers used to extend the assembly") 
    ("solid_abundance_min,x", po::value<int>()->default_value(2), "minimum coverage of a solid k-mer to start the assembly")
    ("solid_abundance_max,X", po::value<int>()->default_value(1000000), "maximum coverage of a solid k-mer to start the assembly")
    //("dis_deviation_max,d", po::value<int>()->default_value(10), "contigs are connected if the inner distance between them is less than (distance suggested by aligned read) + dis_deviation_max.")
    (",t", po::value<int>()->default_value(16), "number of threads")
    ("output,o", po::value<string>()->default_value("unitigs.fa"), "output contig file name (fasta)");
    
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  
  if (argc==1 || vm.count("help")) {
    cerr << desc << "\n";
    exit(0);
  }
  
  po::notify(vm);

  Params options;
  options.K = vm["-k"].as<int>();
  options.readFileList = vm["input"].as<string>();
  options.cqfFile = vm["cqf"].as<string>();
  options.kmer_abundance_min = vm["abundance_min"].as<int>();
  options.solid_kmer_abundance_min = vm["solid_abundance_min"].as<int>();
  options.solid_kmer_abundance_max = vm["solid_abundance_max"].as<int>();
  //options.dis_deviation_max = vm["dis_deviation_max"].as<int>();
  options.thread_num = vm["-t"].as<int>();
  options.output = vm["output"].as<string>();
  switch(vm["format"].as<char>()){
    case 'g':
      options.fmode = FILE_MODE::GZIP;
      break;
    case 'b':
      options.fmode = FILE_MODE::BZIP2;
      break;
    case 'f':
      options.fmode = FILE_MODE::TEXT;
      break;
    default:
      std::cerr<<"[Error] unrecognized file format "<<vm["format"].as<char>()<<endl;
      break;
  }
  cout<<options;
  return options;
}

struct WorkQueue;
struct AssemblyInfo{
  volatile atomic<size_t> seq_num_;
  volatile atomic<size_t> palindrome_seq_num_;
  volatile atomic<size_t> seq_len_in_total_;
  AssemblyInfo(){
    seq_num_ = 0;
    palindrome_seq_num_ = 0;
    seq_len_in_total_ = 0;
  }
};
class WorkQueue2{
public:
  volatile atomic<size_t> start_;//start index of contigs to handle.
  volatile atomic<size_t> step_;//step size.
  volatile atomic<size_t> count_;//next index of non-empty contigs.
  AssemblyInfo assemblyInfo_;
  //volatile atomic<size_t> end;
  boost::mutex mut_;

  //WorkQueue2(): start_(0), step_(10), count_(1){}
  WorkQueue2(size_t start=0, size_t step = 10, size_t count = 1): start_(start), step_(step), count_(count){}
  bool get_work(const concurrent_vector<Contig>& contigs, int& start, int& end, int& counter){
    if(contigs.size() <= start_){
      return false;
    }
    mut_.lock();
    counter = count_;
    start = start_;
    end = min(start_+step_, contigs.size());
    int tmp =0;
    for(int x = start; x < end; x++){
      if(contigs[x].seq.length()!=0){
        tmp++;
      }
    }
    start_ = end;
    count_ += tmp;
    mut_.unlock();
    return true;
  }
};
class WorkQueue3 : public WorkQueue2{
public:
  queue<size_t> count2output;//maintain the output order.
  boost::mutex write_mut_;
  using WorkQueue2::WorkQueue2;
  //WorkQueue3(size_t start=0, size_t step = 10, size_t count = 1): start_(start), step_(step), count_(count){}
  bool get_work(const concurrent_vector<Contig>& contigs, int& start, int& end, int& counter){
    if(contigs.size() <= start_){
      return false;
    }
    mut_.lock();
    while(start_ < contigs.size() && contigs[start_].seq.dna_base_num() == 0){
      start_++;
    }
    if(start_ == contigs.size()){
      return false;
    }
    counter = count_;
    start = start_;
    //end = start_;
    //end = min(start_+step_, contigs.size());
    int tmp =0;
    for(end = start; end < contigs.size() && tmp < step_; end++){
      if(contigs[end].seq.length()!=0){
        tmp++;
      }
    }
    start_ = end;
    count_ += tmp;
    mut_.unlock();

    write_mut_.lock();
    count2output.push(counter);
    write_mut_.unlock();
    return true;
  }
};

//void find_unitigs_mt_master(CQF_mt& cqf, const vector<string>& seqFiles, const Params& params, vector_mt<string>& contigs);
//void find_unitigs_mt_master(CQF_mt& cqf, seqFile_batch& seqFiles, const Params& options, vector_mt<Contig>& unitigs);
//void find_unitigs_mt_master(CQF_mt& cqf, seqFile_batch& seqFiles, const Params& options, concurrent_vector<Contig>& unitigs);
void find_unitigs_mt_master(CQF_mt& cqf, seqFile_batch& seqFiles, const Params& options, concurrent_vector<Contig>& unitigs, hash_map_mt& startKmer2unitig);

//void find_unitigs_mt_worker(CQF_mt& cqf, const Params& params, vector_mt<string>& contigs, unordered_map_mt<string, int>& startKmer2unitig, WorkQueue* work_queue);
//void find_unitigs_mt_worker(CQF_mt& cqf, const Params& options, vector_mt<Contig>& contigs, unordered_map_mt<string, int>& startKmer2unitig, WorkQueue* work_queue);

//void find_unitigs_mt_worker(CQF_mt& cqf, seqFile_batch& seqFiles, const Params& options, concurrent_vector<Contig>& contigs, unordered_set_mt& startKmer2unitig, WorkQueue* work_queue);
void find_unitigs_mt_worker(CQF_mt& cqf, seqFile_batch& seqFiles, const Params& options, concurrent_vector<Contig>& contigs, hash_map_mt& startKmer2unitig, WorkQueue* work_queue);

//void get_unitig_forward(CQF_mt& cqf, const Params& params, vector_mt<string>& contigs, unordered_map_mt<string, int>& startKmer2unitig, WorkQueue* work_queue, int contig_id);
//void get_unitig_forward(CQF_mt& cqf, const Params& options, vector_mt<Contig>& contigs, unordered_map_mt<string, int>& startKmer2unitig, WorkQueue* work_queue, const int& contig_id);

//void get_unitig_forward(CQF_mt& cqf, const Params& options, vector_mt<Contig>& contigs, unordered_set_mt& startKmer2unitig, WorkQueue* work_queue, const int& contig_id);

//void get_unitig_forward(CQF_mt& cqf, const Params& options, concurrent_vector<Contig>& contigs, unordered_set_mt& startKmer2unitig, WorkQueue* work_queue, concurrent_vector<Contig>::iterator& contigIter);
void get_unitig_forward(CQF_mt& cqf, const Params& options, concurrent_vector<Contig>& contigs, hash_map_mt& startKmer2unitig, WorkQueue* work_queue, concurrent_vector<Contig>::iterator& contigIter);

//void get_unitig_forward(CQF_mt& cqf, const Params& options, vector_mt<Contig>& contigs, unordered_map_mt<string, int>& startKmer2unitig, WorkQueue* work_queue, Contig& contig);

//void get_unitig_forward(CQF_mt& cqf, const Params& options, vector_mt<Contig>& contigs, unordered_set_mt& startKmer2unitig, WorkQueue* work_queue, Contig& contig);

//void get_unitig_backward(CQF_mt& cqf, const Params& params, vector_mt<string>& contigs, unordered_map_mt<string, int>& startKmer2unitig, WorkQueue* work_queue, int contig_id);
//void get_unitig_backward(CQF_mt& cqf, const Params& options, vector_mt<Contig>& contigs, unordered_map_mt<string, int>& startKmer2unitig, WorkQueue* work_queue, int contig_id);

void check_unitig(const Params& options, concurrent_vector<Contig>& contigs, const hash_map_mt& startKmer2unitig, WorkQueue2* work_queue);

void track_kmer_worker(const Params& options, const concurrent_vector<Contig>& contigs, hash_map_mt& startKmer2unitig, WorkQueue2* work_queue);
void build_graph_worker(const Params& options, const concurrent_vector<Contig>& contigs, const hash_map_mt& startKmer2unitig, concurrent_vector<UnitigNode>& unitigNodes, WorkQueue2* work_queue);
void build_graph_worker_continue(const Params& options, const concurrent_vector<Contig>& contigs, const hash_map_mt& startKmer2unitig, concurrent_vector<UnitigNode>& unitigNodes, WorkQueue2* work_queue, const int& last_batch_end_idx);

void extend_worker(seqFile_batch& seqFiles, const Params& options, concurrent_vector<Contig>& contigs, hash_map_mt& startKmer2unitig, concurrent_vector<UnitigNode>& unitigNodes, concurrent_vector<Contig>& more_contigs, hash_map_mt& more_startKmer2contigs);

int main(int argc, char* argv[]){
  // /*test of DNA string 
  if(false){
  // string tmp_s("ATGC");
  // DNAString dna_seq = tmp_s;
  // cout<<tmp_s<<endl;
  // cout<<dna_seq.pop()<<endl;
  // dna_seq.append('T');
  // cout<<dna_seq<<endl;
  // return 0;
  // cout<<dna_seq.substr(0, 5)<<endl;
  // cout<<dna_seq.substr(7,5)<<endl;
  // if(dna_seq.substr(0,5) == dna_seq.substr(7,5)){
  //   cout<<"Equal operator works"<<endl;
  // }
  // dna_seq.RC();
  // cout<<dna_seq<<endl;
  // //cout<<dna_seq<<endl;
  // //cout<<dna_seq.dna_base_num()<<endl;
  
  // DNAString dna_seq = string("AAATT");
  // cout<<dna_seq.is_palindrome()<<endl;
  // cout<<dna_seq.is_hairpin()<<endl;
  
  // unordered_map<DNAString, int> kmerSet;
  
  // DNAString dna_seq = string("AAATTT");
  
  // char DNA_bases[4]={'A','C','G','T'};
  // DNAString fix = dna_seq.substr(0, 5);
  // DNAString dna_seq1;
  // for(int x = 0; x < 4; x++){
  //   dna_seq1 = dna_seq.substr(0,5).append(DNA_bases[x]);
  //   kmerSet[dna_seq1]=0;
  //   cout<<"Insert "<<dna_seq1<<endl;
  // }

  // for(int x = 0; x < 4; x++){
  //   dna_seq1 = fix;
  //   dna_seq1.append(DNA_bases[x]);
  //   auto it = kmerSet.find(dna_seq1);
  //   if(it != kmerSet.end()){
  //     cout<<dna_seq1<<" found."<<endl;
  //     cout<<it->first<<": "<<it->second<<endl;
  //   }
  // }
  
  // string seqs[4] = {"AAAAAA", "T", "GGGGGGGGGG", "CC"};
  // for (auto seq : seqs){
  //   DNAString dna = string(seq);
  //   cout<<dna<<" is simple "<<dna.is_simple()<<endl;
  // }
    return 0;
  }

  Params options = get_opts(argc, argv);
  
  vector<string> seqFileNames;
  ifstream fin;
  fin.open(options.readFileList, ios::in);
  string line;
  while(getline(fin, line)){
    if(line.empty())
      continue;
    seqFileNames.push_back(line);
  }
  FILE_TYPE ftype=FILE_TYPE::FASTQ;
  seqFile_batch seqFiles(seqFileNames, ftype, options.fmode);   

//*
  DisplayCurrentDateTime(); 
  cout<<"[CQF] load cqf from disk"<<endl;
  CQF_mt cqf_mt;
  cqf_mt.load(options.cqfFile);
  cout<<"[CQF] cqf loaded!"<<endl;

  DisplayCurrentDateTime();
  cout<<"[Unitig] find unitigs"<<endl<<std::flush;
  concurrent_vector<Contig> contigs;
  contigs.resize(1);
  hash_map_mt startKmer2unitig;
  find_unitigs_mt_master(cqf_mt, seqFiles, options, contigs, startKmer2unitig);

  //Make sure that there is no duplicate seq.
  DisplayCurrentDateTime();
  cout<<"[Check] make sure that there is no duplicate seq."<<endl<<std::flush;
  WorkQueue2* work_queue4check_unitig = new WorkQueue2(1, 10, 1);
  boost::thread_group prod_threads;
  for(int t = 0; t<options.thread_num; t++){
    prod_threads.add_thread(new boost::thread(check_unitig, boost::ref(options), boost::ref(contigs), boost::ref(startKmer2unitig), work_queue4check_unitig));
  }
  prod_threads.join_all();
  free(work_queue4check_unitig);
  // vector<bool> contig_notDup(contigs.size(), false);
  // for(auto kmer2unitig : startKmer2unitig){
  //   contig_notDup[kmer2unitig.second]=true;
  // }
  // for(int x = 1; x < contigs.size(); x++){
  //   if(contigs[x].seq.dna_base_num()>0 && !contig_notDup[x]){
  //     contigs[x].clear();
  //   }
  // }

  WorkQueue2* work_queue4trace_kmer = new WorkQueue2(1, 10, 1);
  // boost::thread_group prod_threads;
  for(int t = 0; t<options.thread_num; t++){
    prod_threads.add_thread(new boost::thread(track_kmer_worker, boost::ref(options), boost::ref(contigs), boost::ref(startKmer2unitig), work_queue4trace_kmer));
  }
  prod_threads.join_all();
  size_t non_duplicate_contig_num = work_queue4trace_kmer->assemblyInfo_.seq_num_;
  cout<<"[Unitig] "<<work_queue4trace_kmer->assemblyInfo_.seq_num_<<" unitigs reported of length "<<work_queue4trace_kmer->assemblyInfo_.seq_len_in_total_<<" bp in total"<<endl;
  cout<<"[Unitig] among them, there are "<<work_queue4trace_kmer->assemblyInfo_.palindrome_seq_num_<<" palindromes."<<endl;
  free(work_queue4trace_kmer);

//*/

  //check whether can reconstruct the full sequences from the graph.
/*
  if(false){
    string refFile = "/public/hshi/tools/SH-assembly/test/case1/genome10K.fasta";
    ifstream refFin;
    refFin.open(refFile, ios::in);
    string refSeq="";
    getline(refFin, line);
    while(getline(refFin, line)){
      refSeq += line;
    }
    to_upper_DNA(refSeq);
    cout<<"[Test] length of reference sequence: "<<refSeq.length()<<" bp."<<endl;
    int constructed_len=0;
    DNAString kmer=refSeq.substr(0,options.K);
    // while(true){
    //   auto it = startKmer2unitig.find(kmer);
    //   if(it == startKmer2unitig.end()){
    //     cout<<"[Test] kmer not found: "<<kmer<<endl;  
    //     if(constructed_len != 0){
    //       cout<<"[Test] last kmer in the traversed sequence: "<<refSeq.substr(constructed_len-1, options.K)<<endl;
    //       cout<<"[Test] traversal of reference seq broken at "<<constructed_len+options.K-1<<endl;
    //     }
    //     break;
    //   }
    //   //get the real ID
    //   auto id = abs(it->second);
    //   while(id2id_afterRemoveDuplicate[id] < abs(it->second)){
    //     id++;
    //   }
    //   constructed_len += (contigs[id].seq.length() - options.K +1);
    //   cout<<"[Test] find kmer: "<<kmer<<" in contig "<<it->second<<" of length "<<contigs[id].seq.length()<<" bp."<<endl;
    //   cout<<"[Test] contig "<< it->second<<": ";
    //   if(it->second < 0)
    //     cout<<contigs[id].seq.get_RC();
    //   else
    //     cout<<contigs[id].seq;
    //   cout<<" median_abundance: "<<contigs[id].median_abundance<<endl;
    //   kmer = refSeq.substr(constructed_len, options.K);
    //   if(constructed_len+options.K-1 >= refSeq.length()){
    //     cout<<"[Test] reference seq is encoded in the DBG of unitigs"<<endl;
    //     break;
    //   }
    // }

    // //check whether k-mer are in the CQF with proper coverage
    // kmer = refSeq.substr(constructed_len, options.K);
    // uint64_t kmer_hash, kmer_RC_hash, kmer_count;
    // kmer_hash = NTPC64(kmer, options.K, kmer_hash, kmer_RC_hash);
    // kmer_count = cqf_mt.count(kmer_hash%cqf_mt.qf->metadata->range);
    // cout<<"[Test] count: "<<kmer_count<<" of kmer: "<<kmer<<endl;
  }
*/

/*
  //load unitigs from file
  cout<<"[Msg] load unitig graph from file"<<endl;
  concurrent_vector<Contig> contigs;
  hash_map_mt startKmer2unitig;
  concurrent_vector<UnitigNode> unitigNodes;
  //if(!load_unitig_graph(options, options.output, contigs, startKmer2unitig, unitigNodes)){
  if(!load_unitig_graph(options.output, contigs)){
    cerr<<"[Error] failed loading unitig graph."<<endl;
    return 0;
  }
  cout<<"[Msg] "<<contigs.size()-1<<" unitig are loaded."<<endl;
  
  WorkQueue2* work_queue4trace_kmer = new WorkQueue2(1, 10, 1);
  boost::thread_group prod_threads;
  for(int t = 0; t<options.thread_num; t++){
    prod_threads.add_thread(new boost::thread(track_kmer_worker, boost::ref(options), boost::ref(contigs), boost::ref(startKmer2unitig), work_queue4trace_kmer));
  }
  prod_threads.join_all();
  free(work_queue4trace_kmer);

  unitigNodes.resize(contigs.size());
  WorkQueue2* work_queue4build_graph = new WorkQueue2(1, 10, 1);
  //boost::thread_group prod_threads;
  for(int t = 0; t<options.thread_num; t++){
    prod_threads.add_thread(new boost::thread(build_graph_worker, boost::ref(options), boost::ref(contigs), boost::ref(startKmer2unitig), boost::ref(unitigNodes), work_queue4build_graph));
  }
  prod_threads.join_all();
  free(work_queue4build_graph);
*/

/*
//trace the contigs produced by Minia in the unitig graph
if(true){
  //load unitig graph produced by Minia
  cout<<"[Msg] load unitig graph from file for Minia"<<endl;
  concurrent_vector<Contig> minia_unitigs;
  hash_map_mt minia_startKmer2unitig;
  concurrent_vector<UnitigNode> minia_unitigNodes;
  string minia_unitig_file = "../Minia-contig-500G/files.list.unitigs.fa";
  //if(!load_unitig_graph(options, minia_unitig_file, minia_unitigs, minia_startKmer2unitig, minia_unitigNodes)){
  if(!load_unitig_graph(minia_unitig_file, minia_unitigs)){
    cerr<<"[Error] failed loading unitig graph."<<endl;
    return 0;
  }
  cout<<"[Msg] "<<minia_unitigs.size()-1<<" unitig are loaded."<<endl;
  
  work_queue4trace_kmer = new WorkQueue2(1, 10, 1);
  for(int t = 0; t<options.thread_num; t++){
    prod_threads.add_thread(new boost::thread(track_kmer_worker, boost::ref(options), boost::ref(minia_unitigs), boost::ref(minia_startKmer2unitig), work_queue4trace_kmer));
  }
  prod_threads.join_all();
  free(work_queue4trace_kmer);

  minia_unitigNodes.resize(minia_unitigs.size());
  work_queue4build_graph = new WorkQueue2(1, 10, 1);
  for(int t = 0; t<options.thread_num; t++){
    prod_threads.add_thread(new boost::thread(build_graph_worker, boost::ref(options), boost::ref(minia_unitigs), boost::ref(minia_startKmer2unitig), boost::ref(minia_unitigNodes), work_queue4build_graph));
  }
  prod_threads.join_all();
  free(work_queue4build_graph);

  string header, seq;
  DNAString kmer;
  vector<int> unitig2SH_contig(contigs.size(), -1);
  vector<int> SH_contig_lens(1, 0);
  size_t matched_seg_id = 0;
  fstream SH_fin;
  string SH_contig_file = "unitigs.contigs.fa";
  SH_fin.open(SH_contig_file, ios::in);
  while(getline(SH_fin, header)){
    getline(SH_fin, seq);
    if(seq.length() < 500){
      continue;
    }
    matched_seg_id++;
    SH_contig_lens.push_back(seq.length());
    int seq_idx = 0;
    hash_map_mt::const_accessor const_access;
    kmer = seq.substr(0, options.K);
    size_t last_matched_loc;
    bool isGap = true;
    while(seq_idx <= seq.length() - options.K){
      if(startKmer2unitig.find(const_access, kmer)){
        // if(isGap){
        //   isGap = false;
        //   last_matched_loc = seq_idx;
        //   matched_seg_id++;
        // }
        unitig2SH_contig[abs(const_access->second)] = matched_seg_id;
        seq_idx += contigs[abs(const_access->second)].seq.length() - options.K + 1;
        if(seq_idx <= seq.length() -options.K){
          kmer = seq.substr(seq_idx, options.K);
        }
        const_access.release();
      }else{
        // if(!isGap){
        //   isGap = true;
        // }
        seq_idx++;
        if(seq_idx <= seq.length() - options.K){
          kmer.pop();
          kmer.append(seq[seq_idx+options.K-1]);
        }
      }
    }
  }
  SH_fin.close();

  int gap_num = 0;
  vector<size_t> withoutGap_contig_lens;
  vector<size_t> withGap_contig_lens;

  string Minia_contig = "../Minia-contig-500G/minia.fa";
  fstream minia_fin;
  minia_fin.open(Minia_contig, ios::in);
  while(getline(minia_fin, header)){
    getline(minia_fin, seq);
    if(seq.length() < 500){
      continue;
    }
    withoutGap_contig_lens.push_back(seq.length());
    int seq_idx = 0;
    vector<array<int, 2>> path;
    hash_map_mt::const_accessor const_access;
    kmer = seq.substr(0, options.K);
    size_t last_matched_loc;
    bool isGap = true;
    while(seq_idx <= seq.length() - options.K){
      if(startKmer2unitig.find(const_access, kmer)){
        if(isGap && unitig2SH_contig[abs(const_access->second)] != -1){
          isGap = false;
          last_matched_loc = seq_idx;
          matched_seg_id = unitig2SH_contig[abs(const_access->second)];
        }
        if(matched_seg_id != -1 && matched_seg_id != unitig2SH_contig[abs(const_access->second)]){
          bool isAlternative = false;
          if(unitig2SH_contig[abs(const_access->second)] == -1){//could be due to alternative paths
            auto last_matched_node = (*path.begin())[1];
            if(last_matched_node < 0){
              for(auto ele : unitigNodes[-last_matched_node].beforeNodes){
                if(unitig2SH_contig[abs(ele.toNode)] == matched_seg_id){
                  isAlternative = true;
                  break;
                }
              }
            }else{
              for(auto ele : unitigNodes[last_matched_node].afterNodes){
                if(unitig2SH_contig[abs(ele.toNode)] == matched_seg_id){
                  isAlternative = true;
                  break;
                }
              }
            }
          }
          if(!isAlternative){
            cout<<"There is a breakpoint, where "<<(*path.rbegin())[1]<<"("<<contigs[abs((*path.rbegin())[1])].median_abundance<<") should connect to "<<const_access->second<<"("<<contigs[abs(const_access->second)].median_abundance<<")"<<endl;
            //search in the unitig graph of Minia
            DNAString kmer_tmp;
            kmer_tmp = RC_DNA(seq.substr(seq_idx-1, options.K));
            hash_map_mt::const_accessor const_access_tmp;
            if(minia_startKmer2unitig.find(const_access_tmp, kmer_tmp)){
              cout<<"The corresponding node in Minia's unitig graph is "<<-const_access_tmp->second<<endl;
            }else{
              cout<<"Why?"<<endl;
            }
            //to proceed
            if(unitig2SH_contig[abs(const_access->second)] == -1){
              isGap = true;
              matched_seg_id = -1;
            }else{
              last_matched_loc = seq_idx;
              matched_seg_id = unitig2SH_contig[abs(const_access->second)];
            }
          }
        }
        path.push_back({seq_idx, const_access->second});
        seq_idx += contigs[abs(const_access->second)].seq.length() - options.K + 1;
        if(seq_idx <= seq.length() -options.K){
          kmer = seq.substr(seq_idx, options.K);
        }
        const_access.release();
      }else{
        if(!isGap){
          isGap = true;
          withGap_contig_lens.push_back(seq_idx - last_matched_loc);
        }
        seq_idx++;
        if(seq_idx <= seq.length() - options.K){
          kmer.pop();
          kmer.append(seq[seq_idx+options.K-1]);
        }
      }
    }
    if(!isGap){
      withGap_contig_lens.push_back(seq.length() - last_matched_loc);
    }else{
      gap_num ++;
    }
    //output the path
    // cout<<"Alignment of seq: "<<header<<endl;
    // int start, end; start = 0;
    // for(int x = 0; x < path.size(); x++){
    //   if(path[x][0] > start){
    //     cout<<start<<"-"<<path[x][0]<<"(gap) ";
    //     start = path[x][0];
    //   }else if(path[x][0] < start){
    //     cout<<"[error] unexpected."<<endl;
    //   }
    //   if(path[x][0] == start){
    //     end = contigs[abs(path[x][1])].seq.length() - options.K + 1 + start;
    //     cout<<start<<"-"<<std::min(int(seq.length()), end)<<"(match) ";
    //     start = end;
    //   }
    // }
    // cout<<endl;
  }
  minia_fin.close();
  gap_num += (withGap_contig_lens.size() - withoutGap_contig_lens.size());
  cout<<"[Info] "<<gap_num<<" gaps are found inside "<<withoutGap_contig_lens.size()<<" contigs."<<endl
  <<"[Info] Minia assembly NG50 :"<<NGx(withoutGap_contig_lens, 50, 3209286105)<<endl
  <<"[Info] After removing the section not in our unitig graph, NG50: "<<NGx(withGap_contig_lens, 50, 3209286105)<<endl;
}  
*/

//*
  //build unitig graph
  cout<<"[Unitig] build unitig graph."<<endl;
  DisplayCurrentDateTime();
  concurrent_vector<UnitigNode> unitigNodes(non_duplicate_contig_num+1);//store the connection info of each node
  WorkQueue2* work_queue4build_graph = new WorkQueue2(1, 10, 1);
  for(int t = 0; t<options.thread_num; t++){
    prod_threads.add_thread(new boost::thread(build_graph_worker, boost::ref(options), boost::ref(contigs), boost::ref(startKmer2unitig), boost::ref(unitigNodes), work_queue4build_graph));
  }
  prod_threads.join_all();
  free(work_queue4build_graph);

  startKmer2unitig.clear();

  //output unitig graph to file
  cout<<"[Dump] save the unitig graph to file."<<endl;
  DisplayCurrentDateTime();
  ofstream fout;
  fout.open(options.output, ios::out);
  size_t noduplicate_contig_num=0;
  for(int contig_id = 1; contig_id<contigs.size(); contig_id++){
    if(contigs[contig_id].seq.dna_base_num()==0){
      continue;
    }
    fout<<">"<<noduplicate_contig_num<<" LN:i:"<<contigs[contig_id].seq.length()<<" KC:i:"<<contigs[contig_id].median_abundance*(contigs[contig_id].seq.length()-options.K+1)<<" km:f:"<<contigs[contig_id].median_abundance;
    for(auto tmp : unitigNodes[noduplicate_contig_num+1].afterNodes){
      if(tmp>0){
        fout<<" L:+:"<<tmp - 1<<":+";
      }else{
        fout<<" L:+:"<<-tmp - 1<<":-";
      }
    }
  
    for(auto tmp : unitigNodes[noduplicate_contig_num+1].beforeNodes){
      if(tmp>0){
        fout<<" L:-:"<<tmp-1<<":+";
      }else{
        fout<<" L:-:"<<-tmp-1<<":-";
      }
    }
    fout<<endl<<contigs[contig_id].seq<<endl;
    noduplicate_contig_num++;
  }
  fout.close();
//*/

//pseudo-align reads and extend
/*
if(false){
  // //remove the duplicates
  // cout<<"[Post-process] resize contigs by removing duplicates."<<endl;
  // DisplayCurrentDateTime();
  // size_t noduplicate_contig_num = 1;
  // for(int idx = 1; idx < contigs.size(); idx++){
  //   if(contigs[idx].seq.dna_base_num()==0){
  //     continue;
  //   }
  //   if(idx > noduplicate_contig_num){
  //     contigs[noduplicate_contig_num] = contigs[idx];
  //   }
  //   noduplicate_contig_num++;
  // }
  // contigs.resize(noduplicate_contig_num);

  //reset the median_abudances of each unitig to 0
  for(auto& contig : contigs){
    contig.median_abundance = 0;
  }

  //pseudo-align reads to the compacted de Bruijn graph
  cout<<"[Align] pseudo align reads to compacted DBG..."<<endl;
  DisplayCurrentDateTime();
  seqFile_batch seqFiles_1(seqFileNames, ftype, options.fmode);
  concurrent_vector<Contig> more_contigs;//(1);
  hash_map_mt more_startKmer2contigs;
  for(int t = 0; t<options.thread_num; t++){
    prod_threads.add_thread(new boost::thread(extend_worker, boost::ref(seqFiles_1), boost::ref(options), boost::ref(contigs), boost::ref(startKmer2unitig), boost::ref(unitigNodes), boost::ref(more_contigs), boost::ref(more_startKmer2contigs)));
  }
  prod_threads.join_all();

  //add the new contigs
  int more_contigs_len_total=0;
  size_t noduplicate_contig_num = contigs.size();
  int single_duplicate, both_duplicate; single_duplicate = both_duplicate = 0;
  int duplicate_nonequal = 0;
  for(auto& contig : more_contigs){
    DNAString startKmer, RC_startKmer;
    startKmer = contig.seq.substr(0, options.K);
    RC_startKmer = contig.seq.substr(contig.seq.length() - options.K); RC_startKmer.RC();
    hash_map_mt::accessor access1, access2;
    hash_map_mt::const_accessor const_access1, const_access2;
    bool left_dup, right_dup; left_dup = right_dup = false;
    if(startKmer2unitig.find(const_access1, startKmer)){
      left_dup = true;
      if(const_access1->second < 0){
        if(contig.seq != contigs[-const_access1->second].seq.get_RC()){
          duplicate_nonequal++;
        }
      }else{
        if(contig.seq != contigs[const_access1->second].seq){
          duplicate_nonequal++;
        }
      }
      const_access1.release();
    }
    if(startKmer2unitig.find(const_access2, RC_startKmer)){
      right_dup = true;
      if(!left_dup){
        if(const_access2->second < 0){
          if(contig.seq != contigs[-const_access2->second].seq){
            duplicate_nonequal++;
          }
        }else{
          if(contig.seq != contigs[const_access2->second].seq.get_RC()){
            duplicate_nonequal++;
          }
        }
      }
      const_access2.release();
    }
    if(left_dup && right_dup){
      both_duplicate ++;
    }else if(left_dup || right_dup){
      single_duplicate ++;
    }else{
      startKmer2unitig.insert(access1, startKmer);
      access1->second = noduplicate_contig_num;
      access1.release();
      startKmer2unitig.insert(access2, RC_startKmer);
      access2->second = -noduplicate_contig_num;
      access2.release();
      contigs.push_back(contig);
      noduplicate_contig_num ++;
      more_contigs_len_total += contig.seq.length();
    }
  }
  cout<<"[Info] among "<<more_contigs.size()<<" novel contigs, "<<both_duplicate+single_duplicate<<" are duplicates with "
  <<single_duplicate<<" having duplicate starting k-mer at one end, and "
  <<both_duplicate<<" having duplicate starting k-mers at both ends. "
  <<"There are "<<duplicate_nonequal<<" contigs identified as duplicates, but not completely equal to its duplicates."<<endl;
  cout<<"[Info] "<<more_contigs.size() - both_duplicate - single_duplicate<<" non-duplicate novel contigs of length "<<more_contigs_len_total<<" bp in total are introduced based on pseudo-align."<<endl;

  //build the graph for the novel contigs
  int last_batch_end_idx = unitigNodes.size();
  WorkQueue2* work_queue4build_graph = new WorkQueue2(last_batch_end_idx, 10, last_batch_end_idx);
  unitigNodes.resize(contigs.size());
  boost::thread_group prod_threads;
  for(int t = 0; t<options.thread_num; t++){
    prod_threads.add_thread(new boost::thread(build_graph_worker_continue, boost::ref(options), boost::ref(contigs), boost::ref(startKmer2unitig), boost::ref(unitigNodes), work_queue4build_graph, last_batch_end_idx));
  }
  prod_threads.join_all();
  free(work_queue4build_graph);

  //output new graph
  vector<int> idx2newidx(contigs.size());
  int unitigRemoved = 0;
  int unitigRemovedLen=0;
  for(int x = 1, y=0; x < contigs.size(); x++){
    if(contigs[x].median_abundance > 0){
      idx2newidx[x] = y++;
    }else{
      unitigRemoved++;
      unitigRemovedLen += contigs[x].seq.length();
    }
  }
  cout<<"[Prune] "<<unitigRemoved<<" unitigs of length "<<unitigRemovedLen<<" bp with zero abundance are removed."<<endl;
  //output unitig graph to file
  DisplayCurrentDateTime();
  cout<<"[Dump] save the pruned unitig graph to file."<<endl;
  ofstream fout;
  string new_output_file = "unitigs.new.fa";
  fout.open(new_output_file, std::ios::out);
  // size_t 
  noduplicate_contig_num=0;
  for(int contig_id = 1; contig_id<contigs.size(); contig_id++){
    if(contigs[contig_id].median_abundance == 0){
      continue;
    }
    fout<<">"<<noduplicate_contig_num<<" LN:i:"<<contigs[contig_id].seq.length()<<" KC:i:"<<contigs[contig_id].median_abundance*(contigs[contig_id].seq.length()-options.K+1)<<" km:f:"<<contigs[contig_id].median_abundance;
    for(auto tmp : unitigNodes[contig_id].afterNodes){
      //if(tmp.coverage ==0 || contigs[abs(tmp.toNode)].median_abundance == 0) continue;
      if(contigs[abs(tmp.toNode)].median_abundance == 0) continue;
      if(tmp.toNode>0){
        fout<<" L:+:"<< idx2newidx[tmp.toNode]<<":+";
      }else{
        fout<<" L:+:"<< idx2newidx[-tmp.toNode]<<":-";
      }
    }
    for(auto tmp : unitigNodes[contig_id].beforeNodes){
      //if(tmp.coverage == 0 || contigs[abs(tmp.toNode)].median_abundance == 0) continue;
      if(contigs[abs(tmp.toNode)].median_abundance == 0) continue;
      if(tmp.toNode>0){
        fout<<" L:-:"<<idx2newidx[tmp.toNode]<<":+";
      }else{
        fout<<" L:-:"<<idx2newidx[-tmp.toNode]<<":-";
      }
    }
    fout<<endl<<contigs[contig_id].seq<<endl;
    noduplicate_contig_num++;
  }
  fout.close();
}
*/
  cout<<"[Dump] saved!"<<endl;
  DisplayCurrentDateTime();  
  //contig_summary(contigs);
  return 0;
}

//use 1-based counter, e.g. the first work has work_id 1.
/*
struct WorkQueue{
  volatile atomic<uint32_t> next_work;
  volatile atomic<uint32_t> total_work;
  volatile atomic<uint32_t> work_done;
  volatile atomic<bool> master_done;
  boost::mutex mut;

  WorkQueue(){
    next_work = 1;
    total_work = 0;
    work_done = 0;
    master_done = false;
  }
  WorkQueue(uint32_t n, uint32_t t){
    assert(n>0 && t>0);
    next_work = n;
    total_work = t;
    work_done = n-1;
    master_done = false;
  }
  
  bool get_next_work(uint32_t& work_id){
    boost::unique_lock<boost::mutex> lock(mut);
    if(next_work <= total_work){
      work_id = next_work;
      next_work++; 
      lock.unlock();
      return true;
    }else{
      lock.unlock();
      return false;
    }
  }
  void report_work_done(int num=1){
    work_done+=num;
  }

  void add_work(uint32_t work_num){
    total_work += work_num;
  }
  void add_skip_work(uint32_t work_num){
    boost::unique_lock<boost::mutex> lock(mut);
    next_work += work_num;
    total_work += work_num;
    work_done += work_num;
    lock.unlock();
  }
};
*/

struct WorkQueue{
  concurrent_queue<concurrent_vector<Contig>::iterator> jobQueue;
  volatile atomic<uint32_t> total_work;
  volatile atomic<uint32_t> work_done;
  volatile atomic<bool> master_done;
  //boost::mutex mut;

  WorkQueue(){
    total_work = 0;
    work_done = 0;
    master_done = false;
  }
  //WorkQueue(uint32_t t){
  //  assert(t>0);
    //next_work = n;
  //   total_work = t;
  //   work_done = ;
  //   master_done = false;
  // }
  
  bool get_next_work(concurrent_vector<Contig>::iterator& jobIter){
    if(jobQueue.try_pop(jobIter)){
      return true;
    }else{
      return false;
    }
  }
  void report_work_done(int num=1){
    work_done+=num;
  }

  void add_work(concurrent_vector<Contig>::iterator jobIter){
    jobQueue.push(jobIter);
    total_work ++;
  }
};

/*
void track_kmer_worker(const Params& options, const concurrent_vector<Contig>& contigs, hash_map_mt& startKmer2unitig, WorkQueue2* work_queue){
  int start, end, counter, totalUnitigs_len, seq_num, palindrome_contig_num;
  totalUnitigs_len = seq_num = palindrome_contig_num = 0;
  DNAString first_kmer,last_kmer_RC;
  while(work_queue->get_work(contigs, start, end, counter)){
    for(size_t contig_id = start; contig_id < end; contig_id++){
      if(contigs[contig_id].seq.dna_base_num()!=0){
        seq_num++;
        totalUnitigs_len += contigs[contig_id].seq.length();
        first_kmer = contigs[contig_id].seq.substr(0,options.K);
        last_kmer_RC = contigs[contig_id].seq.substr(contigs[contig_id].seq.length()-options.K).RC();
        
        if(first_kmer == last_kmer_RC){
          palindrome_contig_num++;
          //startKmer2unitig[first_kmer] = noduplicate_contig_num;
          hash_map_mt::accessor access;
          if(startKmer2unitig.find(access, first_kmer)){
            access->second = counter;
            access.release();
          }else{
            cerr<<"[Error] kmer not found!"<<endl;
          }
          //dont't know how to do  
          //cout<<"[Warning] contig "<<contig_id<<" seems like a palindrome seq (+ and - strand): "<<contigs[contig_id].seq<<endl;
        }else{
          hash_map_mt::accessor access;
          if(startKmer2unitig.find(access, first_kmer)){
            access->second = counter;
            access.release();
          }else{
            cerr<<"[Error] kmer not found!"<<endl;
          }
          if(startKmer2unitig.find(access, last_kmer_RC)){
            access->second = -counter;
            access.release();
          }else{
            cerr<<"[Error] kmer not found!"<<endl;
          }
          //startKmer2unitig[first_kmer] = noduplicate_contig_num;
          //startKmer2unitig[last_kmer_RC] = -noduplicate_contig_num;
        }
        counter++;
      }
    }
  }
  work_queue->assemblyInfo_.seq_num_ += seq_num;
  work_queue->assemblyInfo_.palindrome_seq_num_ += palindrome_contig_num;
  work_queue->assemblyInfo_.seq_len_in_total_ += totalUnitigs_len;
}
*/
void check_unitig(const Params& options, concurrent_vector<Contig>& contigs, const hash_map_mt& startKmer2unitig, WorkQueue2* work_queue){
  int start, end, counter;
  DNAString first_kmer;
  while(work_queue->get_work(contigs, start, end, counter)){
    for(size_t contig_id = start; contig_id < end; contig_id++){
      if(contigs[contig_id].seq.dna_base_num()!=0){
        first_kmer = contigs[contig_id].seq.substr(0,options.K);
        hash_map_mt::const_accessor access;
        if(startKmer2unitig.find(access, first_kmer)){
          if(access->second != contig_id){
            contigs[contig_id].clear();
          }
          access.release();
        }else{
          cerr<<"[Error] kmer not found!"<<endl;
        }
      }
    }
  }
}

void track_kmer_worker(const Params& options, const concurrent_vector<Contig>& contigs, hash_map_mt& startKmer2unitig, WorkQueue2* work_queue){
  int start, end, counter, totalUnitigs_len, seq_num, palindrome_contig_num;
  totalUnitigs_len = seq_num = palindrome_contig_num = 0;
  DNAString first_kmer,last_kmer_RC;
  while(work_queue->get_work(contigs, start, end, counter)){
    for(size_t contig_id = start; contig_id < end; contig_id++){
      if(contigs[contig_id].seq.dna_base_num()!=0){
        seq_num++;
        totalUnitigs_len += contigs[contig_id].seq.length();
        first_kmer = contigs[contig_id].seq.substr(0,options.K);
        last_kmer_RC = contigs[contig_id].seq.substr(contigs[contig_id].seq.length()-options.K).RC();
        
        if(first_kmer == last_kmer_RC){
          palindrome_contig_num++;
          hash_map_mt::accessor access;
          // startKmer2unitig.insert(access, first_kmer);
          // access->second = counter;
          // access.release();
          if(startKmer2unitig.find(access, first_kmer)){
            access->second = counter;
            access.release();
          }else{
            cerr<<"[Error] kmer not found!"<<endl;
          }
          //dont't know how to do  
          //cout<<"[Warning] contig "<<contig_id<<" seems like a palindrome seq (+ and - strand): "<<contigs[contig_id].seq<<endl;
        }else{
          hash_map_mt::accessor access;
          // startKmer2unitig.insert(access, first_kmer);
          // access->second = counter;
          // access.release();
          // startKmer2unitig.insert(access, last_kmer_RC);
          // access->second = -counter;
          // access.release();
          if(startKmer2unitig.find(access, last_kmer_RC)){
            access->second = -counter;
            access.release();
          }else{
            cerr<<"[Error] kmer not found!"<<endl;
          }
          if(startKmer2unitig.find(access, first_kmer)){
            access->second = counter;
            access.release();
          }else{
            cerr<<"[Error] kmer not found!"<<endl;
          }
        }
        counter++;
      }
    }
  }
  work_queue->assemblyInfo_.seq_num_ += seq_num;
  work_queue->assemblyInfo_.palindrome_seq_num_ += palindrome_contig_num;
  work_queue->assemblyInfo_.seq_len_in_total_ += totalUnitigs_len;
}

void build_graph_worker(const Params& options, const concurrent_vector<Contig>& contigs, const hash_map_mt& startKmer2unitig, concurrent_vector<UnitigNode>& unitigNodes, WorkQueue2* work_queue){
  DNAString kmer_fix, kmer;
  char bases[4]={'A','C','G','T'};
  int start, end, counter;//, first_count;
  while(work_queue->get_work(contigs, start, end, counter)){
    //first_count = counter; //mark the start of the counter.
    //stringstream buf;
    for(int contig_id = start; contig_id<end; contig_id++){
      if(contigs[contig_id].seq.dna_base_num()==0){
        continue;
      }
      //buf<<">"<<counter-1<<" LN:i:"<<contigs[contig_id].seq.length()<<" KC:i:"<<contigs[contig_id].median_abundance*(contigs[contig_id].seq.length()-options.K+1)<<" km:f:"<<contigs[contig_id].median_abundance;
      kmer_fix=contigs[contig_id].seq.substr(contigs[contig_id].seq.length()-options.K+1, options.K-1);

      for(int x=0; x<4; x++){
        kmer = kmer_fix; kmer.append(bases[x]);
        hash_map_mt::const_accessor const_access;
        if(startKmer2unitig.find(const_access, kmer)){
          auto tmp = const_access->second;
          unitigNodes[counter].afterNodes.push_back(tmp);
          // if(tmp>0){
          //   buf<<" L:+:"<<tmp-1<<":+";
          //   //if(isHairpin[tmp]){
          //   //  fout<<" L:+:"<<tmp-1<<":-";
          //   //}
          // }else{
          //   buf<<" L:+:"<<-tmp-1<<":-";
          // }
          const_access.release();
        }
      }
      // if(unitigNodes[counter].afterNodes.size()){
      //   sort(unitigNodes[counter].afterNodes.begin(), unitigNodes[counter].afterNodes.end(), [](const Edge& lhs, const Edge& rhs){
      // return abs(lhs.toNode) < abs(rhs.toNode);
      //   });
      // }

      kmer_fix=contigs[contig_id].seq.substr(0, options.K-1).RC();
      for(int x=3; x >= 0; x--){
        kmer = kmer_fix; kmer.append(bases[x]);
        hash_map_mt::const_accessor const_access;
        if(startKmer2unitig.find(const_access, kmer)){
          auto tmp = const_access->second;
          unitigNodes[counter].beforeNodes.push_back(tmp);
          // if(tmp>0){
          //   buf<<" L:-:"<<tmp-1<<":+";
          //   //if(isHairpin[tmp]){
          //   //  fout<<" L:-:"<<tmp-1<<":-";
          //   //}
          // }else{
          //   buf<<" L:-:"<<-tmp-1<<":-";
          // }
          const_access.release();
        }
      }
      // if(unitigNodes[counter].beforeNodes.size()){
      //   //fix the bug in Minia
      //   sort(unitigNodes[counter].beforeNodes.begin(), unitigNodes[counter].beforeNodes.end(), [](const Edge& lhs, const Edge& rhs){
      //     return abs(lhs.toNode) < abs(rhs.toNode);
      //   });
      // }
      //buf<<endl<<contigs[contig_id].seq<<endl;
      counter++;
    }
    // while(work_queue->count2output.front() != first_count){
    //   boost::this_thread::sleep_for(boost::chrono::milliseconds(20));
    // }
    // work_queue->write_mut_.lock();
    // fout<<buf.str();
    // work_queue->count2output.pop();
    // work_queue->write_mut_.unlock();
  }
}

/*
void build_graph_worker_continue(const Params& options, const concurrent_vector<Contig>& contigs, const hash_map_mt& startKmer2unitig, concurrent_vector<UnitigNode>& unitigNodes, WorkQueue2* work_queue, const int& last_batch_end_idx){
  DNAString kmer_fix, kmer;
  char bases[4]={'A','C','G','T'};
  int start, end, counter;//, first_count;
  while(work_queue->get_work(contigs, start, end, counter)){
    //first_count = counter; //mark the start of the counter.
    //stringstream buf;
    for(int contig_id = start; contig_id<end; contig_id++){
      if(contigs[contig_id].seq.dna_base_num()==0){
        continue;
      }
      //buf<<">"<<counter-1<<" LN:i:"<<contigs[contig_id].seq.length()<<" KC:i:"<<contigs[contig_id].median_abundance*(contigs[contig_id].seq.length()-options.K+1)<<" km:f:"<<contigs[contig_id].median_abundance;
      kmer_fix=contigs[contig_id].seq.substr(contigs[contig_id].seq.length()-options.K+1, options.K-1);

      for(int x=0; x<4; x++){
        kmer = kmer_fix; kmer.append(bases[x]);
        hash_map_mt::const_accessor const_access;
        if(startKmer2unitig.find(const_access, kmer)){
          auto tmp = const_access->second;
          unitigNodes[counter].afterNodes.push_back(Edge(tmp));
          if(abs(tmp) < last_batch_end_idx){
            if(tmp > 0){
              unitigNodes[tmp].beforeNodes.push_back(Edge(-counter));
            }else{
              unitigNodes[-tmp].afterNodes.push_back(Edge(-counter));
            }
          }
          // if(tmp>0){
          //   buf<<" L:+:"<<tmp-1<<":+";
          //   //if(isHairpin[tmp]){
          //   //  fout<<" L:+:"<<tmp-1<<":-";
          //   //}
          // }else{
          //   buf<<" L:+:"<<-tmp-1<<":-";
          // }
          const_access.release();
        }
      }
      kmer_fix=contigs[contig_id].seq.substr(0, options.K-1).RC();
      for(int x=0; x<4; x++){
        kmer = kmer_fix; kmer.append(bases[x]);
        hash_map_mt::const_accessor const_access;
        if(startKmer2unitig.find(const_access, kmer)){
          auto tmp = const_access->second;
          unitigNodes[counter].beforeNodes.push_back(Edge(tmp));
          if(abs(tmp) < last_batch_end_idx){
            if(tmp > 0){
              unitigNodes[tmp].beforeNodes.push_back(Edge(counter));
            }else{
              unitigNodes[-tmp].afterNodes.push_back(Edge(counter));
            }
          }
          // if(tmp>0){
          //   buf<<" L:-:"<<tmp-1<<":+";
          //   //if(isHairpin[tmp]){
          //   //  fout<<" L:-:"<<tmp-1<<":-";
          //   //}
          // }else{
          //   buf<<" L:-:"<<-tmp-1<<":-";
          // }
          const_access.release();
        }
      }
      //buf<<endl<<contigs[contig_id].seq<<endl;
      counter++;
    }
    // while(work_queue->count2output.front() != first_count){
    //   boost::this_thread::sleep_for(boost::chrono::milliseconds(20));
    // }
    // work_queue->write_mut_.lock();
    // fout<<buf.str();
    // work_queue->count2output.pop();
    // work_queue->write_mut_.unlock();
  }
}
*/

/*
void processDataChunk_extend(const Params& options, concurrent_vector<Contig>& contigs, hash_map_mt& startKmer2unitig, concurrent_vector<UnitigNode>& unitigNodes, concurrent_vector<Contig>& more_contigs, hash_map_mt& more_startKmer2contigs, chunk& dataChunk);
void extend_worker(seqFile_batch& seqFiles, const Params& options, concurrent_vector<Contig>& contigs, hash_map_mt& startKmer2unitig, concurrent_vector<UnitigNode>& unitigNodes, concurrent_vector<Contig>& more_contigs, hash_map_mt& more_startKmer2contigs){
  chunk dataChunk;
  while(seqFiles.getDataChunk(dataChunk)){
    processDataChunk_extend(options, contigs, startKmer2unitig, unitigNodes, more_contigs, more_startKmer2contigs, dataChunk);
  }
}
bool insert_or_replace(hash_map_mt& startKmer2unitig, const DNAString& dnastr, const size_t idx);
bool is_connected(const Params& options, const concurrent_vector<Contig>& contigs, const concurrent_vector<UnitigNode>& unitigNodes, const int startNode, const int endNode, const int dis_in_read);
*/
// void processDataChunk_extend(const Params& options, concurrent_vector<Contig>& contigs, hash_map_mt& startKmer2unitig, concurrent_vector<UnitigNode>& unitigNodes, concurrent_vector<Contig>& more_contigs, hash_map_mt& more_startKmer2contigs, chunk& dataChunk){
//   std::string line, seq;
//   while(dataChunk.readLine(line)){
//     if(line.empty()){
//       continue;
//     }else if(line[0] != '@'){
//       continue;
//     }
//     if(!dataChunk.readLine(seq)){
//       break;
//     }
//     //seq = "AAAATTTATAGCTAACGCCAAATTTCTTTGGGTCAGTTTCAATGTTTACCTCAAACTTGGGTAATTAAACCGACTTGAACGGCGCGTCTCTGCGCCTAAC";
//     if(seq.length()<options.K){
//       dataChunk.skipLines(2);//skip two lines after the sequence line
//       continue;
//     }
//     DNAString kmer;
//     vector<array<int, 2>> path;
    
//     int seq_idx = 0;
//     for(int x = seq_idx; x < (seq_idx + options.K) && x < seq.length(); x++){
//       if(seq[x] == 'n' || seq[x]=='N'){
//         seq_idx = x+1;
//         //break;
//       }
//     }
//     if(seq_idx <= seq.length()-options.K){
//       kmer = seq.substr(seq_idx, options.K);
//     }//multithread_io::mtx.lock();
//     while(seq_idx <= seq.length()-options.K){ //cout<<seq_idx<<":";
//       hash_map_mt::const_accessor const_access;
//       if(startKmer2unitig.find(const_access, kmer)){
//         auto tmp = const_access->second;
//         const_access.release();
//         //see whether the RC(k-mer before) is also recorded
//         if(seq_idx > 0 && seq[seq_idx-1] != 'n' && seq[seq_idx-1] != 'N'){
//           DNAString kmer_tmp = seq.substr(seq_idx-1, options.K);
//           kmer_tmp.RC();
//           if(startKmer2unitig.find(const_access, kmer_tmp)){
//             auto contig_idx  = const_access->second;
//             path.push_back(array<int, 2>{-seq_idx, contig_idx});//-seq_idx indicate RC of k-mer ends before seq_idx 
//             const_access.release();
//           }
//         }
//         path.push_back(array<int, 2>{seq_idx, tmp});
//         //jump to next k-mer
//         seq_idx += (contigs[abs(tmp)].seq.dna_base_num() - options.K + 1);
//       }else if(seq_idx+options.K < seq.length()){
//         if(seq[seq_idx + options.K] != 'n' && seq[seq_idx + options.K] != 'N'){
//           kmer.pop().append(seq[seq_idx + options.K]);
//           seq_idx++;
//           continue;
//         }else{
//           seq_idx = seq_idx + options.K + 1;
//         }
//       }else{
//         break;
//       }
//       if(seq_idx + options.K - 1 >= seq.length()){
//         break;      
//       }
//       for(int x = seq_idx; x < (seq_idx + options.K) && x < seq.length(); x++){
//         if(seq[x] == 'n' || seq[x]=='N'){
//           x++;
//           while(x < seq.length() && (seq[x] == 'n' || seq[x]=='N')){
//             x++;
//           }
//           seq_idx = x;
//         }
//       }
//       if(seq_idx <= seq.length()-options.K){
//         kmer = seq.substr(seq_idx, options.K);
//       }
//     }//cout<<"seq_idx"<<endl;
    
//     //Process the alignment
//     //add node coverage and connections
//     /*
//     for(int x = 0; x < path.size(); x++){
//       if(path[x][0]<0){
//         int A, B;
//         A = abs(path[x][1]);
//         B = abs(path[x+1][1]);
//         if(path[x][1] < 0){
//           for(auto& node : unitigNodes[A].afterNodes){
//             if(abs(node.toNode) == B){
//               node.coverage++;
//             }
//           }
//         }else{
//           for(auto& node : unitigNodes[A].beforeNodes){
//             if(abs(node.toNode) == B){
//               node.coverage++;
//             }
//           }
//         }
//         if(path[x+1][1] < 0){
//           for(auto& node : unitigNodes[B].afterNodes){
//             if(abs(node.toNode) == A){
//               node.coverage++;
//             }
//           }
//         }else{
//           for(auto& node : unitigNodes[B].beforeNodes){
//             if(abs(node.toNode) == A){
//               node.coverage++;
//             }
//           }
//         }
//         if(x==0 || path[x-1][1] != -path[x][1]){
//           ++(contigs[abs(path[x][1])].median_abundance);
//         }
//       }else{
//         ++(contigs[abs(path[x][1])].median_abundance);
//       }
//     }
//     */
//     for(int x = 0; x < path.size(); x++){
//       if(path[x][0] >= 0){
//         ++(contigs[abs(path[x][1])].median_abundance);
//       }else if(path[x-1][1] != -path[x][1]){
//         ++(contigs[abs(path[x][1])].median_abundance);
//       }
//     }

//     //trying to connect unitigs
//     // int new_contig_num = 0;
//     //*
//     if(path.size()){
//       int last_contig_idx = 0;
//       int seq_toMatch = 0;//next sequence to match
//       int path_idx = 0;
//       if(path[0][0] >= 0){
//         path_idx = 1;
//       }
//       while(path_idx < path.size()){
//         if(path[path_idx][0] < 0){
//           path_idx += 2;
//         }else{
//           last_contig_idx = path[path_idx-1][1];
//           seq_toMatch = path[path_idx-1][0] + contigs[abs(last_contig_idx)].seq.length() - options.K + 1;

//           if(seq_toMatch >= path[path_idx][0]){
//             path_idx ++;
//             continue;
//           }

//           // DNAString startKmer;
//           // startKmer = seq.substr(seq_toMatch, options.K);
//           // hash_map_mt::const_accessor const_access;
//           // if(more_startKmer2contigs.find(const_access, startKmer)){
//           //   path_idx ++;
//           //   const_access.release();
//           //   continue;
//           // }

//           if(!is_connected(options, contigs, unitigNodes, last_contig_idx, path[path_idx][1], path[path_idx][0] - seq_toMatch)){
//             DNAString startKmer;
//             startKmer = seq.substr(seq_toMatch, options.K);
//             hash_map_mt::const_accessor const_access;
//             if(more_startKmer2contigs.find(const_access, startKmer)){
//               path_idx ++;
//               const_access.release();
//               continue;
//             }
//             //see whether already exist
//             // bool isDuplicate=false;
//             // DNAString startKmer, RC_startKmer;
//             DNAString RC_startKmer;
//             // startKmer = seq.substr(seq_toMatch, options.K);    
//             // hash_map_mt::const_accessor const_access;
//             // if(more_startKmer2contigs.find(const_access, startKmer)){
//             //   seq_idx ++;
//             //   const_access.release();
//             //   continue;
//             // }
//             RC_startKmer = seq.substr(path[path_idx][0]-1, options.K); 
//             RC_startKmer.RC();
//             // if(more_startKmer2contigs.find(const_access, RC_startKmer)){
//             //   seq_idx ++;
//             //   const_access.release();
//             //   continue;
//             // }

//             //int contig_id_new = contigs.size();
//             auto it = more_contigs.push_back(Contig(seq.substr(seq_toMatch, path[path_idx][0] - seq_toMatch + options.K - 1), 1));
//             hash_map_mt::accessor access;
//             if(more_startKmer2contigs.insert(access, startKmer)){
//               access.release();
//             }
//             if(more_startKmer2contigs.insert(access, RC_startKmer)){
//               access.release();
//             }
//             //int contig_id_new = it - more_contigs.begin();
//             //unitigNodes.push_back(UnitigNode());
//             //trace kmer
//             // DNAString kmer;
//             // kmer = it->seq.substr(0, options.K);
//             // hash_map_mt::accessor access;
//             // if(startKmer2unitig.insert(access, kmer)){
//             //   access->second = contig_id_new;
//             //   access.release();
//             // }else{
//             //   cerr<<"[Error] failed to insert k-mer "<<kmer<<" for contig "<<contig_id_new<<endl;
//             // }
//             // kmer = it->seq.substr(it->seq.length() - options.K).get_RC();
//             // if(startKmer2unitig.insert(access, kmer)){
//             //   access->second = -contig_id_new;
//             //   access.release();
//             // }else{
//             //   cerr<<"[Error] failed to insert k-mer "<<kmer<<" for contig "<<-contig_id_new<<endl;
//             // }
//             // //add connections
//             // if(last_contig_idx < 0){
//             //   unitigNodes[-last_contig_idx].beforeNodes.push_back(Edge(contig_id_new, 1));
//             //   unitigNodes[contig_id_new].beforeNodes.push_back(Edge(-last_contig_idx, 1));
//             // }else{
//             //   unitigNodes[last_contig_idx].afterNodes.push_back(Edge(contig_id_new, 1));
//             //   unitigNodes[contig_id_new].beforeNodes.push_back(Edge(-last_contig_idx, 1));
//             // }
//           }
//           path_idx ++;
//         }
//       }

//       // if(path[0][0] >= 0){
//       //   last_contig_idx = path[0][1];
//       //   seq_toMatch = path[0][0] + contigs[abs(last_contig_idx)].seq.length() - options.K + 1;
//       //   path_idx = 1;
//       // }else{
//       //   last_contig_idx = path[1][1];
//       //   seq_toMatch = path[1][0] + contigs[abs(last_contig_idx)].seq.length() - options.K + 1;
//       //   path_idx = 2;
//       // }
//       // //int last_matched_loc = -1;
//       // while(path_idx < path.size()){
//       //   if(path[path_idx][0] < 0){
//       //     if(path[path_idx-1][1] == -path[path_idx][1]){
//       //       path_idx++;
//       //       continue;
//       //     }else{//not found
//       //       int tmp = -path[path_idx][0] - (contigs[abs(path[path_idx][1])].seq.length() - options.K + 1);
//       //       if(tmp == seq_toMatch){
//       //         last_contig_idx  = path[path_idx][1];
//       //         seq_toMatch  = -path[path_idx][0];
//       //         path_idx++;
//       //       }else if(tmp > seq_toMatch){//GAP
//       //         //tmp is not the start of path[path_idx][1], which could be caused by sequencing error at the start k-mer.
//       //         last_contig_idx = path[path_idx][1];
//       //         seq_toMatch = -path[path_idx][0];// + contigs[abs(last_contig_idx)].seq.length() - options.K + 1;
//       //         path_idx++;
//       //         // cout<<"[Warning] unexpected case "<<endl;
//       //       }else{//TODO: what if tmp < seq_toMatch, which is unexpected
//       //         last_contig_idx = path[path_idx][1];
//       //         seq_toMatch = -path[path_idx][0];// + contigs[abs(last_contig_idx)].seq.length() - options.K + 1;
//       //         path_idx++;
//       //         // cout<<"[Warning] unexpected case "<<endl;
//       //       }
//       //     }
//       //   }else{
//       //     if(seq_toMatch == path[path_idx][0]){
//       //       last_contig_idx = path[path_idx][1];
//       //       seq_toMatch += contigs[abs(last_contig_idx)].seq.length() - options.K + 1;
//       //       path_idx++;
//       //     }else if(seq_toMatch < path[path_idx][0]){//GAP
//       //       if(!is_connected(options, contigs, unitigNodes, last_contig_idx, path[path_idx][1], path[path_idx][0] - seq_toMatch)){//introduce a new unitig
//       //         //new_contig_num++;
              
//       //         mt_log(std::to_string(new_contig_num)+" new contig added.\n");
//       //         mt_log("between node "+std::to_string(last_contig_idx)+" in "+ std::to_string(path_idx-1)+" and next node.\n");
//       //         multithread_io::mtx.lock();
//       //         cout<<"#sequence "<<seq.length()<<" bp = ";
//       //         for(auto ele : path){
//       //           cout<<"\t"<<ele[0]<<":"<<ele[1]<<"("<<contigs[abs(ele[1])].seq.length()<<" bp)\n";
//       //           // <<"\nKmer:  ";
//       //           // if(ele[0]<0){
//       //           //   cout<<RC_DNA(seq.substr(-ele[0], options.K))<<endl;
//       //           // }else{
//       //           //   cout<<seq.substr(ele[0], options.K)<<endl;
//       //           // }
//       //           // if(ele[1] < 0){
//       //           //   cout<<"Contig:"<<contigs[-ele[1]].seq.get_RC();
//       //           // }else{
//       //           //   cout<<"Contig:"<<contigs[ele[1]].seq;
//       //           // }
//       //           // cout<<endl;
//       //         }
//       //         cout<<endl<<endl<<std::flush;
//       //         multithread_io::mtx.unlock();
              

//       //         int contig_id_new = contigs.size();
//       //         auto it = contigs.push_back(Contig(seq.substr(seq_toMatch, path[path_idx][0] - seq_toMatch + options.K - 1), 1));
//       //         unitigNodes.push_back(UnitigNode());
//       //         //trace kmer
//       //         DNAString kmer;
//       //         kmer = it->seq.substr(0, options.K);
//       //         hash_map_mt::accessor access;
//       //         if(startKmer2unitig.insert(access, kmer)){
//       //           access->second = contig_id_new;
//       //           access.release();
//       //         }else{
//       //           cerr<<"[Error] failed to insert k-mer "<<kmer<<" for contig "<<contig_id_new<<endl;
//       //         }
//       //         kmer = it->seq.substr(it->seq.length() - options.K).get_RC();
//       //         if(startKmer2unitig.insert(access, kmer)){
//       //           access->second = -contig_id_new;
//       //           access.release();
//       //         }else{
//       //           cerr<<"[Error] failed to insert k-mer "<<kmer<<" for contig "<<-contig_id_new<<endl;
//       //         }
//       //         //add connections
//       //         if(last_contig_idx < 0){
//       //           unitigNodes[-last_contig_idx].beforeNodes.push_back(Edge(contig_id_new, 1));
//       //           unitigNodes[contig_id_new].beforeNodes.push_back(Edge(-last_contig_idx, 1));
//       //         }else{
//       //           unitigNodes[last_contig_idx].afterNodes.push_back(Edge(contig_id_new, 1));
//       //           unitigNodes[contig_id_new].beforeNodes.push_back(Edge(-last_contig_idx, 1));
//       //         }
//       //       }
//       //       last_contig_idx = path[path_idx][1];
//       //       seq_toMatch = path[path_idx][0] + contigs[abs(last_contig_idx)].seq.length() - options.K + 1;
//       //       path_idx++;
//       //     }else{//TODO: what if seq_toMatch > path[path_idx][0], which is unexpected.
//       //       last_contig_idx = path[path_idx][1];
//       //       seq_toMatch = path[path_idx][0] + contigs[abs(last_contig_idx)].seq.length() - options.K + 1;
//       //       path_idx++;
//       //       // cout<<"[Warning] unexpected case "<<endl;
//       //     }
//       //   }
//       // }

//     }
//     //*/

//     //multithread_io::mtx.lock();
//     // cout<<"#sequence "<<seq.length()<<" bp = ";
//     // for(auto ele : path){
//     //   cout<<ele[0]<<":"<<ele[1]<<"("<<contigs[abs(ele[1])].seq.length()<<" bp)\t";
//     //   // <<"\nKmer:  ";
//     //   // if(ele[0]<0){
//     //   //   cout<<RC_DNA(seq.substr(-ele[0], options.K))<<endl;
//     //   // }else{
//     //   //   cout<<seq.substr(ele[0], options.K)<<endl;
//     //   // }
//     //   // if(ele[1] < 0){
//     //   //   cout<<"Contig:"<<contigs[-ele[1]].seq.get_RC();
//     //   // }else{
//     //   //   cout<<"Contig:"<<contigs[ele[1]].seq;
//     //   // }
//     //   // cout<<endl;
//     // }
//     // cout<<endl<<endl<<std::flush;
//     // multithread_io::mtx.unlock();
//     //skip two lines
//     dataChunk.skipLines(2);
//   }
//   //free data chunk
//   free(dataChunk.get_reads());
// }


/*
void find_unitigs_mt_master(CQF_mt& cqf, const vector<string>& seqFiles, const Params& params, vector_mt<string>& contigs){
  //initialized worker, assume contigs.size()>=1
  WorkQueue* work_queue = new WorkQueue(contigs.size(), contigs.size());
  unordered_map_mt<string, int> startKmer2unitig(1000);

  boost::thread_group prod_threads;
  for(int t = 0; t<params.thread_num; t++){
    prod_threads.add_thread(new boost::thread(find_unitigs_mt_worker, boost::ref(cqf), boost::ref(params), boost::ref(contigs), boost::ref(startKmer2unitig), work_queue));
  }
  
  string kmer, kmer_RC;
  uint64_t kmer_hash, kmer_RC_hash;
  uint64_t kmer_count, kmer_RC_count;
  size_t contig_id;
  
  for(auto file : seqFiles){
    cout<<"Processing "<<file<<endl;
    std::ifstream fin(file, std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
    inbuf.push(boost::iostreams::gzip_decompressor());
    inbuf.push(fin);
    //Convert streambuf to istream
    std::istream instream(&inbuf);
    //Iterate lines
    std::string line, seq;
    while(std::getline(instream, line)){
      //skip the header line
      if(line.empty()){
        continue;
      }else if(line[0] != '@'){
        continue;
      }
      if(!std::getline(instream, seq)){
        break;
      }
      //cout<<seq<<endl;
      if(seq.length()<params.K){
        continue;
      }
      //seq_RC = RC_DNA(seq);
      int seq_len = seq.length();
      int step = seq_len/3;
      for(int x = 0; x<=seq_len-params.K; x += step){
        kmer = seq.substr(x, params.K);
        kmer_RC = RC_DNA(kmer);//seq_RC.substr(seq_len-params.K-x, params.K);
        if(kmer.find_first_of("nN")!=string::npos){
          continue;
        }
        //kmer_hash = MurmurHash2(kmer);
        //kmer_RC_hash = MurmurHash2(kmer_RC);
        NTPC64(kmer.c_str(), params.K, kmer_hash, kmer_RC_hash);
        
        if(cqf.count_key_value_set_traveled(kmer_hash%cqf.qf->metadata->range, kmer_count) && cqf.count_key_value_set_traveled(kmer_RC_hash%cqf.qf->metadata->range, kmer_RC_count)){
          continue;
        }else if(kmer_count < params.CONTIG_min_cov || kmer_RC_count < params.CONTIG_min_cov){
          continue;
        }

        contig_id = contigs.size_mt();
        startKmer2unitig.insert_mt(kmer, contig_id);
        startKmer2unitig.insert_mt(kmer_RC, -contig_id);
        work_queue->add_skip_work(1);    
        contigs.push_back_mt(seq);

        get_unitig_forward(cqf, params, contigs, startKmer2unitig, work_queue, contig_id);
        get_unitig_backward(cqf, params, contigs, startKmer2unitig, work_queue, contig_id);

        //wait for the end of this run.
        while(work_queue->next_work != work_queue->total_work){
          std::this_thread::sleep_for(std::chrono::milliseconds(1000));//sleep for 1 seconds  
        }
      }
      //skip two lines
      std::getline(instream, line);
      std::getline(instream, line);
    }
  }
  work_queue->master_done = true;
  prod_threads.join_all();
}
*/
/*
void find_unitigs_mt_master(CQF_mt& cqf, seqFile_batch& seqFiles, const Params& options, vector_mt<Contig>& contigs){
  WorkQueue* work_queue = new WorkQueue();
  unordered_map_mt<string, int> startKmer2unitig(1000);

  boost::thread_group prod_threads;
  for(int t = 0; t<options.thread_num; t++){
    prod_threads.add_thread(new boost::thread(find_unitigs_mt_worker, boost::ref(cqf), boost::ref(options), boost::ref(contigs), boost::ref(startKmer2unitig), work_queue));
  }
  
  string kmer, kmer_RC;
  uint64_t kmer_hash, kmer_RC_hash;
  uint64_t kmer_count, kmer_RC_count;
  size_t contig_id;
  
  vector<Contig> master_contigs;
  chunk dataChunk;
  while(seqFiles.getDataChunk(dataChunk)){
    std::string line, seq;
    while(dataChunk.readLine(line)){
      //skip the header line
      if(line.empty()){
        continue;
      }else if(line[0] != '@'){
        continue;
      }
      if(!dataChunk.readLine(seq)){
        break;
      }
      //cout<<seq<<endl;
      if(seq.length()<options.K){
        continue;
      }
      //seq_RC = RC_DNA(seq);
      int seq_len = seq.length();
      int step = seq_len/3;
      for(int x = 0; x<=seq_len- options.K; x += step){
        kmer = seq.substr(x, options.K);
        kmer = to_upper_DNA(kmer);
        //kmer_RC = RC_DNA(kmer);//seq_RC.substr(seq_len-params.K-x, params.K);
        if(kmer.find_first_of("nN")!=string::npos){
          continue;
        }
        kmer_hash = NTPC64(kmer.c_str(), options.K, kmer_hash, kmer_RC_hash);

        if(cqf.count_key_value_set_traveled(kmer_hash%cqf.qf->metadata->range, kmer_count)){
          continue;
        }else if(kmer_count < options.kmer_abundance_min){
          continue;
        }

        //contig_id = contigs.push_back_mt(Contig(kmer, kmer_count));
        //work_queue->add_skip_work(1);
        //startKmer2unitig.insert_mt(kmer, contig_id); //may not be the start k-mer
        startKmer2unitig.insert_mt(kmer, 0);//fake id 
        //startKmer2unitig.insert_mt(kmer_RC, -contig_id); //may not be the end k-mer
        Contig master_contig(kmer, kmer_count);
        get_unitig_forward(cqf, options, contigs, startKmer2unitig, work_queue, master_contig);
        master_contig.seq = RC_DNA(master_contig.seq);
        get_unitig_forward(cqf, options, contigs, startKmer2unitig, work_queue, master_contig);
        master_contigs.push_back(master_contig);
        //std::this_thread::sleep_for(std::chrono::milliseconds(1000));
        break;
        //get_unitig_forward(cqf, options, contigs, startKmer2unitig, work_queue, contig_id);
        //contigs.set_mt(contig_id, Contig(RC_DNA(contigs[contig_id].seq), contigs[contig_id].median_abundance));
        //get_unitig_forward(cqf, options, contigs, startKmer2unitig, work_queue, contig_id);
        //get_unitig_backward(cqf, options, contigs, startKmer2unitig, work_queue, contig_id);

        //wait for the end of this run.
        //while(work_queue->work_done != work_queue->total_work){
        //  std::this_thread::sleep_for(std::chrono::milliseconds(1000));//sleep for 1 seconds  
        //}
      }
      //skip two lines
      dataChunk.skipLines(2);
    }
    //if load detected is less than half. start processing a new chunk
    while(work_queue->work_done != work_queue->total_work){
    //while(work_queue->total_work/work_queue->work_done) > ){
      if(work_queue->next_work > work_queue->total_work){
        cout<<"[test] stucked.."<<endl;
      }
      cout<<"[test] "<<work_queue->next_work<<":"<<work_queue->work_done<<":"<<work_queue->total_work<<endl;
      std::this_thread::sleep_for(std::chrono::milliseconds(1000));//sleep for 1 seconds  
    }
  }
  work_queue->master_done = true;
  prod_threads.join_all();

  contigs.insert(contigs.end(), master_contigs.begin(), master_contigs.end());
}
*/
/*
void processDataChunk(CQF_mt& cqf, const Params& options, concurrent_vector<Contig>& contigs, unordered_set_mt& startKmers, WorkQueue* work_queue, chunk& dataChunk){
  string kmer;
  uint64_t kmer_hash, kmer_RC_hash;
  uint64_t kmer_count, kmer_RC_count;
  std::string line, seq;
  while(dataChunk.readLine(line)){
    //skip the header line
    if(line.empty()){
      continue;
    }else if(line[0] != '@'){
      continue;
    }
    if(!dataChunk.readLine(seq)){
      break;
    }
    if(seq.length()<options.K){
      dataChunk.skipLines(2);//skip two lines after the sequence line
      continue;
    }

    int seq_len = seq.length();
    int middle = seq_len/2;

    if(middle <= seq_len- options.K){
      kmer = seq.substr(middle, options.K);
      if(kmer.find_first_of("nN")!=string::npos){
        continue;
      }
      kmer_hash = NTPC64(kmer.c_str(), options.K, kmer_hash, kmer_RC_hash);

      if(cqf.count_key_value_set_traveled(kmer_hash%cqf.qf->metadata->range, kmer_count)){
        continue;
      }else if(kmer_count < options.solid_kmer_abundance_min || kmer_count > options.solid_kmer_abundance_max){
        continue;
      }

      auto iter = contigs.push_back(Contig(kmer, kmer_count));
      to_upper_DNA(kmer);
      startKmers.insert(kmer);//fake id 

      get_unitig_forward(cqf, options, contigs, startKmers, work_queue, iter);
      iter->seq.RC();
      get_unitig_forward(cqf, options, contigs, startKmers, work_queue, iter);
    }
    //skip two lines
    dataChunk.skipLines(2);
  }
  //free data chunk
  free(dataChunk.get_reads());
}
*/

bool is_connected(const Params& options, const concurrent_vector<Contig>& contigs, const concurrent_vector<UnitigNode>& unitigNodes, const int startNode, const int endNode, const int dis_in_read){
  if(startNode < 0){
    if(unitigNodes[-startNode].beforeNodes.size()){
      return false;
    }
  }else{
    if(unitigNodes[startNode].afterNodes.size()){
      return false;
    }
  }
  if(endNode < 0){
    if(unitigNodes[-endNode].afterNodes.size()){
      return false;
    }
  }else{
    if(unitigNodes[endNode].beforeNodes.size()){
      return false;
    }
  }
  return true;
}

/*
bool is_connected(const Params& options, const concurrent_vector<Contig>& contigs, const concurrent_vector<UnitigNode>& unitigNodes, const int startNode, const int endNode, const int dis_in_read){
  int max_dis = dis_in_read + options.dis_deviation_max;
  //breadth-first search
  //map<int, int> toNodeDis_old, toNodeDis;
  vector<array<int, 2>> toNodeDis_old, toNodeDis;
  toNodeDis_old.push_back({startNode, 0});//distance to start node is 0
  while(toNodeDis_old.size()){
    for(auto ele : toNodeDis_old){
      if(ele[0] < 0){
        for(auto edge : unitigNodes[-ele[0]].beforeNodes){
          // if(abs(edge.toNode) == abs(endNode)){
          if(edge.toNode == endNode){
            return true;
          }
          int dis_tmp = ele[1] + (contigs[abs(edge.toNode)].seq.length() - options.K + 1);
          if(dis_tmp <= max_dis){
            toNodeDis.push_back({edge.toNode, dis_tmp});
            // auto it = toNodeDis.find(edge.toNode);
            // if(it == toNodeDis.end()){
            //   toNodeDis.insert(pair<int, int>(edge.toNode, dis_tmp));
            // }else{
            //   it->second = min(it->second, dis_tmp);
            // }
          }
        }
      }else{
        for(auto edge : unitigNodes[ele[0]].afterNodes){
          // if(abs(edge.toNode) == abs(endNode)){
          if(edge.toNode == endNode){
            return true;
          }
          int dis_tmp = ele[1] + (contigs[abs(edge.toNode)].seq.length() - options.K + 1);
          if(dis_tmp <= max_dis){
            toNodeDis.push_back({edge.toNode, dis_tmp});
            // auto it = toNodeDis.find(edge.toNode);
            // if(it == toNodeDis.end()){
            //   toNodeDis.insert(pair<int, int>(edge.toNode, dis_tmp));  
            // }else{
            //   it->second = min(it->second, dis_tmp);
            // }
          }
        }
      }
    }
    toNodeDis_old = std::move(toNodeDis);
    // toNodeDis_old = toNodeDis;
    // toNodeDis.clear();
  }
  return false;
}
*/

void processDataChunk(CQF_mt& cqf, const Params& options, concurrent_vector<Contig>& contigs, hash_map_mt& startKmer2unitig, WorkQueue* work_queue, chunk& dataChunk){
  string kmer;
  uint64_t kmer_hash, kmer_RC_hash;
  uint64_t kmer_count, kmer_RC_count;
  std::string line, seq;
  while(dataChunk.readLine(line)){
    //skip the header line
    if(line.empty()){
      continue;
    }else if(line[0] != '@'){
      continue;
    }
    if(!dataChunk.readLine(seq)){
      break;
    }
    if(seq.length()<options.K){
      dataChunk.skipLines(2);//skip two lines after the sequence line
      continue;
    }

    int seq_len = seq.length();
    int middle = seq_len/2 - options.K/2;

    if(seq_len>=options.K and middle <= seq_len- options.K){
      kmer = seq.substr(middle, options.K);
      to_upper_DNA(kmer);

      if(kmer.find_first_of("nN")!=string::npos){
        continue;
      }
      kmer_hash = NTPC64(kmer.c_str(), options.K, kmer_hash, kmer_RC_hash);

      if(cqf.count_key_value_set_traveled(kmer_hash%cqf.qf->metadata->range, kmer_count)){
        continue;
      }else if(kmer_count < options.solid_kmer_abundance_min || kmer_count > options.solid_kmer_abundance_max){
        continue;
      }

      auto iter = contigs.push_back(Contig(kmer, kmer_count));
      
      size_t contig_id = iter - contigs.begin(); //std::distance(contigs.begin(), iter);
      // bool is_dup = false;
      // for(int x = 0; x<4; x++){
      //   if(kmer == string(options.K, DNA::bases[x])){
      //     hash_map_mt::accessor access;
      //     if(startKmer2unitig.insert(access, DNAString(kmer))){
      //       access->second = contig_id;
      //     }else{
      //       is_dup = true;
      //     }
      //     break;
      //   }
      // }
      // //startKmers.insert(DNAString(kmer));//fake id 
      // if(is_dup){
      //   continue;
      // }

      get_unitig_forward(cqf, options, contigs, startKmer2unitig, work_queue, iter);

      hash_map_mt::const_accessor const_access;
      if(iter->seq.dna_base_num() > 0){
        if(startKmer2unitig.find(const_access, DNAString(kmer))){
          if(const_access->second > contig_id){
            iter->seq.RC(); const_access.release();
            get_unitig_forward(cqf, options, contigs, startKmer2unitig, work_queue, iter);    
          }else if(const_access->second < contig_id){
            iter->clear(); const_access.release();
          }
          // const_access.release();
        }else{
          iter->seq.RC();
          get_unitig_forward(cqf, options, contigs, startKmer2unitig, work_queue, iter);
        }
      }
      // if(iter->seq.dna_base_num() > 0){
      //   iter->seq.RC() ;
      //   get_unitig_forward(cqf, options, contigs, startKmer2unitig, work_queue, iter);
      // }
    }
    //skip two lines
    dataChunk.skipLines(2);
  }
  //free data chunk
  free(dataChunk.get_reads());
}
/*
void find_unitigs_mt_master(CQF_mt& cqf, seqFile_batch& seqFiles, const Params& options, concurrent_vector<Contig>& contigs){
  WorkQueue* work_queue = new WorkQueue();
  unordered_set_mt startKmers;

  boost::thread_group prod_threads;
  for(int t = 0; t<options.thread_num; t++){
    prod_threads.add_thread(new boost::thread(find_unitigs_mt_worker, boost::ref(cqf), boost::ref(seqFiles), boost::ref(options), boost::ref(contigs), boost::ref(startKmers), work_queue));
  }
  
  string kmer, kmer_RC;
  uint64_t kmer_hash, kmer_RC_hash;
  uint64_t kmer_count, kmer_RC_count;
  //size_t contig_id;
  
  //vector<Contig> master_contigs;
  chunk dataChunk;
  while(seqFiles.getDataChunk(dataChunk)){
    std::string line, seq;
    while(dataChunk.readLine(line)){
      //skip the header line
      if(line.empty()){
        continue;
      }else if(line[0] != '@'){
        continue;
      }
      if(!dataChunk.readLine(seq)){
        break;
      }
      //cout<<seq<<endl;
      if(seq.length()<options.K){
        dataChunk.skipLines(2);
        continue;
      }
      //seq_RC = RC_DNA(seq);
      int seq_len = seq.length();
      int middle = seq_len/2;

      if(middle <= seq_len- options.K){
        kmer = seq.substr(middle, options.K);
        //kmer = 
        to_upper_DNA(kmer);
        //kmer_RC = RC_DNA(kmer);//seq_RC.substr(seq_len-params.K-x, params.K);
        if(kmer.find_first_of("nN")!=string::npos){
          continue;
        }
        kmer_hash = NTPC64(kmer.c_str(), options.K, kmer_hash, kmer_RC_hash);

        if(cqf.count_key_value_set_traveled(kmer_hash%cqf.qf->metadata->range, kmer_count)){
          continue;
        }else if(kmer_count < options.solid_kmer_abundance_min || kmer_count > options.solid_kmer_abundance_max){
          continue;
        }

        //contig_id = contigs.push_back_mt(Contig(kmer, kmer_count));
        //work_queue->add_skip_work(1);
        auto iter = contigs.push_back(Contig(kmer, kmer_count));
        //startKmer2unitig.insert_mt(kmer, contig_id); //may not be the start k-mer
        startKmers.insert(DNAString(kmer));//fake id 
        //startKmer2unitig.insert_mt(kmer_RC, -contig_id); //may not be the end k-mer
        // Contig master_contig(kmer, kmer_count);
        // get_unitig_forward(cqf, options, contigs, startKmer2unitig, work_queue, master_contig);
        // master_contig.seq = RC_DNA(master_contig.seq);
        // get_unitig_forward(cqf, options, contigs, startKmer2unitig, work_queue, master_contig);
        // master_contigs.push_back(master_contig);
        //std::this_thread::sleep_for(std::chrono::milliseconds(1000));
        //break;
        get_unitig_forward(cqf, options, contigs, startKmers, work_queue, iter);
        //contigs.set(contig_id, Contig(RC_DNA(contigs[contig_id].seq), contigs[contig_id].median_abundance));
        iter->seq.RC();
        get_unitig_forward(cqf, options, contigs, startKmers, work_queue, iter);
        //get_unitig_backward(cqf, options, contigs, startKmer2unitig, work_queue, contig_id);

        //std::this_thread::sleep_for(std::chrono::milliseconds(300));
        //wait for the end of this run.
      }
      //skip two lines
      dataChunk.skipLines(2);
      //auto tmp_total = work_queue->total_work;
      //auto tmp_done = work_queue->work_done;
      while((work_queue->total_work - work_queue->work_done) > options.thread_num){
        //cout<<" "<<work_queue->work_done<<"/"<<work_queue->total_work<<"("<<work_queue->total_work-work_queue->work_done<<")"<<contigs.size()<<std::flush;
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));//sleep for 1 seconds  
      }
    }
    //free data chunk
    free(dataChunk.get_reads());

    //if load detected is less than half. start processing a new chunk
    //while(work_queue->work_done != work_queue->total_work){
    //while(work_queue->total_work/work_queue->work_done) > ){
      //if(work_queue->next_work > work_queue->total_work){
      //  cout<<"[test] stucked.."<<endl;
      //}
      //cout<<"[test] "<<work_queue->next_work<<":"<<work_queue->work_done<<":"<<work_queue->total_work<<endl;
      //std::this_thread::sleep_for(std::chrono::milliseconds(1000));//sleep for 1 seconds  
    //}
  }
  while(work_queue->total_work > work_queue->work_done){
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  }
  work_queue->master_done = true;
  prod_threads.join_all();

  //contigs.insert(contigs.end(), master_contigs.begin(), master_contigs.end());
  startKmers.clear();
}
*/

void find_unitigs_mt_master(CQF_mt& cqf, seqFile_batch& seqFiles, const Params& options, concurrent_vector<Contig>& contigs, hash_map_mt& startKmer2unitig){
  WorkQueue* work_queue = new WorkQueue();

  boost::thread_group prod_threads;
  for(int t = 0; t<options.thread_num; t++){
    prod_threads.add_thread(new boost::thread(find_unitigs_mt_worker, boost::ref(cqf), boost::ref(seqFiles), boost::ref(options), boost::ref(contigs), boost::ref(startKmer2unitig), work_queue));
  }
  
  string kmer, kmer_RC;
  uint64_t kmer_hash, kmer_RC_hash;
  uint64_t kmer_count, kmer_RC_count;
  size_t contig_id;
  
  //vector<Contig> master_contigs;
  chunk dataChunk;
  while(seqFiles.getDataChunk(dataChunk)){
    std::string line, seq;
    while(dataChunk.readLine(line)){
      //skip the header line
      if(line.empty()){
        continue;
      }else if(line[0] != '@'){
        continue;
      }
      if(!dataChunk.readLine(seq)){
        break;
      }
      //cout<<seq<<endl;
      if(seq.length()<options.K){
        dataChunk.skipLines(2);
        continue;
      }
      //seq_RC = RC_DNA(seq);
      int seq_len = seq.length();
      int middle = seq_len/2;

      if(middle <= seq_len- options.K){
        kmer = seq.substr(middle, options.K);
        //kmer = 
        to_upper_DNA(kmer);
        //kmer_RC = RC_DNA(kmer);//seq_RC.substr(seq_len-params.K-x, params.K);
        if(kmer.find_first_of("nN")!=string::npos){
          continue;
        }
        kmer_hash = NTPC64(kmer.c_str(), options.K, kmer_hash, kmer_RC_hash);

        if(cqf.count_key_value_set_traveled(kmer_hash%cqf.qf->metadata->range, kmer_count)){
          continue;
        }else if(kmer_count < options.solid_kmer_abundance_min || kmer_count > options.solid_kmer_abundance_max){
          continue;
        }

        //contig_id = contigs.push_back_mt(Contig(kmer, kmer_count));
        //work_queue->add_skip_work(1);
        auto iter = contigs.push_back(Contig(kmer, kmer_count));
        contig_id = iter - contigs.begin(); //std::distance(contigs.begin(), iter);
        //startKmer2unitig.insert_mt(kmer, contig_id); //may not be the start k-mer
        /*whether it is a simple sequence
        bool is_dup = false;
        for(int x = 0; x<4; x++){
          if(kmer == string(options.K, DNA::bases[x])){
            hash_map_mt::accessor access;
            if(startKmer2unitig.insert(access, DNAString(kmer))){
              access->second = contig_id;
            }else{
              is_dup = true;
            }
            access.release();
            break;
          }
        }
        //startKmers.insert(DNAString(kmer));//fake id 
        if(is_dup){
          continue;
        }
        */

        //startKmer2unitig.insert_mt(kmer_RC, -contig_id); //may not be the end k-mer
        // Contig master_contig(kmer, kmer_count);
        // get_unitig_forward(cqf, options, contigs, startKmer2unitig, work_queue, master_contig);
        // master_contig.seq = RC_DNA(master_contig.seq);
        // get_unitig_forward(cqf, options, contigs, startKmer2unitig, work_queue, master_contig);
        // master_contigs.push_back(master_contig);
        //std::this_thread::sleep_for(std::chrono::milliseconds(1000));
        //break;
        get_unitig_forward(cqf, options, contigs, startKmer2unitig, work_queue, iter);
        //contigs.set(contig_id, Contig(RC_DNA(contigs[contig_id].seq), contigs[contig_id].median_abundance));
        hash_map_mt::const_accessor const_access;
        if(iter->seq.dna_base_num() > 0){
          if(startKmer2unitig.find(const_access, DNAString(kmer))){
            if(const_access->second > contig_id){
              iter->seq.RC(); const_access.release();
              get_unitig_forward(cqf, options, contigs, startKmer2unitig, work_queue, iter);    
            }else if(const_access->second < contig_id){
              iter->clear(); const_access.release();
            }
            // const_access.release();
          }else{
            iter->seq.RC();
            get_unitig_forward(cqf, options, contigs, startKmer2unitig, work_queue, iter);
          }
        }
        //get_unitig_backward(cqf, options, contigs, startKmer2unitig, work_queue, contig_id);

        //std::this_thread::sleep_for(std::chrono::milliseconds(300));
        //wait for the end of this run.
      }
      //skip two lines
      dataChunk.skipLines(2);
      //auto tmp_total = work_queue->total_work;
      //auto tmp_done = work_queue->work_done;
      while((work_queue->total_work - work_queue->work_done) > options.thread_num){
        //cout<<" "<<work_queue->work_done<<"/"<<work_queue->total_work<<"("<<work_queue->total_work-work_queue->work_done<<")"<<contigs.size()<<std::flush;
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));//sleep for 1 seconds  
      }
    }
    //free data chunk
    free(dataChunk.get_reads());

    //if load detected is less than half. start processing a new chunk
    //while(work_queue->work_done != work_queue->total_work){
    //while(work_queue->total_work/work_queue->work_done) > ){
      //if(work_queue->next_work > work_queue->total_work){
      //  cout<<"[test] stucked.."<<endl;
      //}
      //cout<<"[test] "<<work_queue->next_work<<":"<<work_queue->work_done<<":"<<work_queue->total_work<<endl;
      //std::this_thread::sleep_for(std::chrono::milliseconds(1000));//sleep for 1 seconds  
    //}
  }
  while(work_queue->work_done < work_queue->total_work){//total_work on the right is important.
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  }
  work_queue->master_done = true;
  prod_threads.join_all();
  assert(work_queue->total_work == work_queue->work_done);
  free(work_queue);
  //contigs.insert(contigs.end(), master_contigs.begin(), master_contigs.end());
  //startKmer2unitig.clear();
}

/*
void find_unitigs_mt_worker(CQF_mt& cqf, const Params& params, vector_mt<string>& contigs, unordered_map_mt<string, int>& startKmer2unitig, WorkQueue* work_queue){
  uint32_t contig_id=0;
  while(!work_queue->master_done){
    if(!work_queue->get_next_work(contig_id)){
      std::this_thread::sleep_for(std::chrono::milliseconds(500));//sleep for 0.5 seconds
      continue;
    }
    string seq=contigs.at(contig_id);
    string kmer, kmer_RC;
    kmer = seq;
    kmer_RC = RC_DNA(kmer);
    string start_kmer = seq;
    int cid1, cid2;
    if(startKmer2unitig.find_mt(kmer, cid1)){
      if(startKmer2unitig.find_mt(kmer_RC, cid2)){
        if(cid1==contig_id && cid2==contig_id){//?
          get_unitig_forward(cqf, params, contigs, startKmer2unitig, work_queue, contig_id);
        }else if(cid1== -contig_id && cid2 == -contig_id){//?
          get_unitig_backward(cqf, params, contigs, startKmer2unitig, work_queue, contig_id);
        }
        continue;
      }
      if(cid1==contig_id){
        get_unitig_forward(cqf, params, contigs, startKmer2unitig, work_queue, contig_id);
        continue;
      }
    }else if(startKmer2unitig.find_mt(kmer_RC, cid2)){
      if(cid2==-contig_id){
        get_unitig_backward(cqf, params, contigs, startKmer2unitig, work_queue, contig_id);
      }
    }
  }
}
*/
/*
void find_unitigs_mt_worker(CQF_mt& cqf, const Params& options, vector_mt<Contig>& contigs, unordered_map_mt<string, int>& startKmer2unitig, WorkQueue* work_queue){
  uint32_t contig_id=0;
  while(!work_queue->master_done){
    if(!work_queue->get_next_work(contig_id)){
      std::this_thread::sleep_for(std::chrono::milliseconds(500));//sleep for 0.5 seconds
      continue;
    }
    get_unitig_forward(cqf, options, contigs, startKmer2unitig, work_queue, contig_id);
    work_queue->report_work_done();
  }
}
*/
/*
void find_unitigs_mt_worker(CQF_mt& cqf, const Params& options, vector_mt<Contig>& contigs, unordered_set_mt& startKmers, WorkQueue* work_queue){
  uint32_t contig_id=0;
  while(!work_queue->master_done){
    if(!work_queue->get_next_work(contig_id)){
      std::this_thread::sleep_for(std::chrono::milliseconds(500));//sleep for 0.5 seconds
      continue;
    }
    get_unitig_forward(cqf, options, contigs, startKmers, work_queue, contig_id);
    work_queue->report_work_done();
  }
}
*/
/*
void find_unitigs_mt_worker(CQF_mt& cqf, seqFile_batch& seqFiles, const Params& options, concurrent_vector<Contig>& contigs, unordered_set_mt& startKmers, WorkQueue* work_queue){
  concurrent_vector<Contig>::iterator contigIter;
  chunk dataChunk;
  while(!work_queue->master_done){
    if(!work_queue->get_next_work(contigIter)){
      if(seqFiles.getDataChunk(dataChunk)){
        processDataChunk(cqf, options, contigs, startKmers, work_queue, dataChunk);
      }else{
        std::this_thread::sleep_for(std::chrono::milliseconds(500));//sleep for 0.5 seconds
      }
    }else{
      get_unitig_forward(cqf, options, contigs, startKmers, work_queue, contigIter);
      work_queue->report_work_done();
    }
  }
}
*/

void find_unitigs_mt_worker(CQF_mt& cqf, seqFile_batch& seqFiles, const Params& options, concurrent_vector<Contig>& contigs, hash_map_mt& startKmer2unitig, WorkQueue* work_queue){
  concurrent_vector<Contig>::iterator contigIter;
  chunk dataChunk;
  while(!work_queue->master_done || (work_queue->work_done < work_queue->total_work)){
    if(!work_queue->get_next_work(contigIter)){
      if(seqFiles.getDataChunk(dataChunk)){
        processDataChunk(cqf, options, contigs, startKmer2unitig, work_queue, dataChunk);
      }else{
        std::this_thread::sleep_for(std::chrono::milliseconds(500));//sleep for 0.5 seconds
      }
    }else{
      get_unitig_forward(cqf, options, contigs, startKmer2unitig, work_queue, contigIter);
      work_queue->report_work_done();
    }
  }
}
/*
void get_unitig_forward(CQF_mt& cqf, const Params& params, vector_mt<string>& contigs, unordered_map_mt<string, int>& startKmer2unitig, WorkQueue* work_queue, int contig_id){
  array<bool, 4> candidates_before({false, false, false, false}), candidates_after({false, false, false, false});
  array<string, 4> kmer_befores, kmer_afters, kmer_befores_RC, kmer_afters_RC;
  int candidates_before_num, candidates_after_num;
  int nodes_before_num, nodes_after_num;
  uint64_t kmer_hash, kmer_RC_hash, current_kmer_hash, current_kmer_RC_hash;
  string kmer, kmer_RC, current_kmer;
  uint64_t kmer_count, kmer_RC_count;
  int idx;
  
  string contig_seq = contigs.at_mt(contig_id);
  current_kmer = contig_seq.substr(contig_seq.length()-params.K);

  while(true){
    kmer_afters = kmers_after(current_kmer);
    kmer_befores = kmers_before(kmer_afters[0]);
    candidates_before = candidates_after = {{false, false, false, false}};
    candidates_before_num = candidates_after_num = 0;
    nodes_before_num = nodes_after_num = 0;
    
    for(int x = 0; x<4; x++){
      kmer = kmer_afters[x];
      kmer_RC = RC_DNA(kmer);
      kmer_afters_RC[x] = kmer_RC;

      NTPC64(kmer.c_str(), params.K, kmer_hash, kmer_RC_hash);
      //kmer_hash = MurmurHash2(kmer);
      //kmer_RC_hash = MurmurHash2(kmer_RC);
      if(cqf.count_key_value_set_traveled(kmer_hash%cqf.qf->metadata->range, kmer_count) && cqf.count_key_value_set_traveled(kmer_RC_hash%cqf.qf->metadata->range, kmer_RC_count)){
        if(startKmer2unitig.find_mt(kmer, idx)){
          nodes_after_num++;
        }else if(kmer_count >= params.CONTIG_min_cov && kmer_RC_count >= params.CONTIG_min_cov){
          candidates_after[x] = true;
          candidates_after_num ++;
        }
      }else if(kmer_count >= params.CONTIG_min_cov && kmer_RC_count >= params.CONTIG_min_cov){
        candidates_after[x] = true;
        candidates_after_num ++;
      }
    }
    
    for(int x = 0; x<4; x++){
      if(kmer_befores[x] == current_kmer){
        continue;
      }
      kmer = kmer_befores[x];
      kmer_RC = RC_DNA(kmer);
      kmer_befores_RC[x] = kmer_RC;

      NTPC64(kmer.c_str(), params.K, kmer_hash, kmer_RC_hash);
      //kmer_hash = MurmurHash2(kmer);
      //kmer_RC_hash = MurmurHash2(kmer_RC);
      if(cqf.count_key_value_set_traveled(kmer_hash%cqf.qf->metadata->range, kmer_count) && cqf.count_key_value_set_traveled(kmer_RC_hash%cqf.qf->metadata->range, kmer_RC_count)){
        if(startKmer2unitig.find_mt(kmer_RC, idx)){
          nodes_before_num++;
        }else if(kmer_count >= params.CONTIG_min_cov && kmer_RC_count >= params.CONTIG_min_cov){
          candidates_before[x] = true;
          candidates_before_num++;
        }
      }else if(kmer_count >= params.CONTIG_min_cov && kmer_RC_count >= params.CONTIG_min_cov){
        candidates_before[x] = true;
        candidates_before_num++;
      }
    }
    if((nodes_before_num + candidates_before_num) || (nodes_after_num+candidates_after_num)>1){
      if(startKmer2unitig.find_mt(contig_seq.substr(0, params.K), idx)){
        //A contig has been constructed by another program from RC way
        if(abs(idx) != contig_id){
          contigs.set_mt(contig_id, ""); 
          break;
        }
      }else{
        throw logic_error("Unexpectedly not found start kmer!");
      }
      startKmer2unitig.insert_mt(RC_DNA(current_kmer), -contig_id);
      
      for(int x= 0; x < 4; x++){
        if(candidates_after[x]){
          if(!startKmer2unitig.find_mt(kmer_afters[x], idx)){
            int new_unitig_idx = contigs.push_back_mt(kmer_afters[x]);
            startKmer2unitig.insert_mt(kmer_afters[x], new_unitig_idx);
            work_queue->add_work(1); 
          }
        }
      }
      for(int x = 0; x < 4; x++){
        if(candidates_before[x]){
          if(!startKmer2unitig.find_mt(kmer_befores_RC[x], idx)){
            int new_unitig_idx = contigs.push_back_mt(kmer_befores_RC[x]);
            startKmer2unitig[kmer_befores_RC[x]] = -new_unitig_idx;
            work_queue->add_work(1);
          }
        }
      }
    }else if(candidates_after_num==1){
      for(int x = 0; x<4; x++){
        if(candidates_after[x]){
          current_kmer = kmer_afters[x];
          contig_seq += current_kmer.back();
          break;
        }
      }
      continue;
    }else{
      startKmer2unitig.insert_mt(RC_DNA(current_kmer), contig_id);
    }
    break;
  }
}
*/
/*
void get_unitig_forward(CQF_mt& cqf, const Params& options, vector_mt<Contig>& contigs, unordered_map_mt<string, int>& startKmer2unitig, WorkQueue* work_queue, const int& contig_id){
  array<bool, 4> candidates_before({false, false, false, false}), candidates_after({false, false, false, false});
  //array<string, 4> kmer_befores, kmer_afters, kmer_befores_RC, kmer_afters_RC;
  array<uint64_t, 4> kmer_abundance_befores, kmer_abundance_afters;
  auto abundance_min = options.kmer_abundance_min;
  auto K = options.K;

  int candidates_before_num, candidates_after_num;
  int nodes_before_num, nodes_after_num;
  uint64_t kmer_hash, kmer_RC_hash, current_kmer_hash, current_kmer_RC_hash;
  string kmer, kmer_RC, current_kmer, current_kmer_RC, current_kmer_fix;
  uint64_t kmer_count, kmer_RC_count;
  int idx;
  
  string contig_seq = contigs.at_mt(contig_id).seq;
  current_kmer = contig_seq.substr(contig_seq.length()-K);
  current_kmer_RC = RC_DNA(current_kmer);

  std::vector<int> abundances;
  abundances.push_back(int(contigs.at_mt(contig_id).median_abundance));

  NTPC64(current_kmer.c_str(), K, current_kmer_hash, current_kmer_RC_hash);
  int node_after_x, node_before_x;//useful only when there is only one node after and without candidates after.
  while(true){
    //NTPC64(current_kmer.c_str(), K, current_kmer_hash, current_kmer_RC_hash);
    current_kmer_fix = current_kmer.substr(1);
    //kmer_afters = kmers_after(current_kmer);
    //kmer_befores = kmers_before(kmer_afters[0]);
    candidates_before = candidates_after = {{false, false, false, false}};
    candidates_before_num = candidates_after_num = 0;
    nodes_before_num = nodes_after_num = 0;
    
    //kmers with current_kmer_fix as prefix
    for(int x = 0; x<4; x++){
      //kmer = kmer_afters[x];
      //kmer_RC = RC_DNA(kmer);
      //kmer_afters_RC[x] = kmer_RC;
      kmer_hash = current_kmer_hash; kmer_RC_hash = current_kmer_RC_hash;
      kmer_hash = NTPC64(current_kmer[0], DNA_bases[x], K, kmer_hash, kmer_RC_hash);
      //kmer_hash = MurmurHash2(kmer);
      //kmer_RC_hash = MurmurHash2(kmer_RC);
      if(cqf.count_key_value_set_traveled(kmer_hash%cqf.qf->metadata->range, kmer_count)){
        if(startKmer2unitig.find_mt(current_kmer_fix+DNA_bases[x], idx)){
          nodes_after_num++;
          kmer_abundance_afters[x] = kmer_count;
          node_after_x = x;
        }else if(kmer_count >= abundance_min){
          kmer_abundance_afters[x] = kmer_count;
          candidates_after[x] = true; //possible because of hash collisions
          candidates_after_num ++;
        }
      }else if(kmer_count >= abundance_min){
        kmer_abundance_afters[x] = kmer_count;
        candidates_after[x] = true;
        candidates_after_num ++;
      }
    }
    
    //kmers with RC(current_kmer_fix) as prefix
    NTPC64(current_kmer[0], 'A', options.K, current_kmer_hash, current_kmer_RC_hash);
    kmer = current_kmer_RC;
    for(int x = 0; x<4; x++){
      if(DNA_bases[x] == current_kmer_RC[K-1]){
        continue;
      }
      //kmer = kmer_befores[x];
      //kmer_RC = RC_DNA(kmer);
      //kmer_befores_RC[x] = kmer_RC;
      //kmer_hash = current_kmer_hash; kmer_RC_hash = current_kmer_RC_hash;
      //kmer_hash = NTPC64(current_kmer_RC[0], DNA_bases[x], K, kmer_hash, kmer_RC_hash);
      //kmer = current_kmer_RC;
      kmer[K-1] = DNA_bases[x];
      kmer_hash = current_kmer_hash; kmer_RC_hash = current_kmer_RC_hash;
      kmer_hash = NTPC64('T', DNA_bases[x], options.K, kmer_RC_hash, kmer_hash);
      //kmer_hash = NTPC64(kmer.c_str(), K, kmer_hash, kmer_RC_hash);
      //kmer_hash = MurmurHash2(kmer);
      //kmer_RC_hash = MurmurHash2(kmer_RC);
      if(cqf.count_key_value_set_traveled(kmer_hash%cqf.qf->metadata->range, kmer_count)){
        if(startKmer2unitig.find_mt(kmer, idx)){
          nodes_before_num++; node_before_x = x;
        }else if(kmer_count >= abundance_min){
          kmer_abundance_befores[x] = kmer_count;
          candidates_before[x] = true;
          candidates_before_num++;
        }
      }else if(kmer_count >= abundance_min){
        kmer_abundance_befores[x] = kmer_count;
        candidates_before[x] = true;
        candidates_before_num++;
      }
    }

    if((nodes_before_num + candidates_before_num) || (nodes_after_num+candidates_after_num)>1){ //no-linear extension
      // if(startKmer2unitig.find_mt(contig_seq.substr(0, params.K), idx)){
      //   //A contig has been constructed by another program from RC way
      //   if(abs(idx) != contig_id){
      //     contigs.set_mt(contig_id, ""); 
      //     break;
      //   }
      // }else{
      //   throw logic_error("Unexpectedly not found start kmer!");
      // }
      startKmer2unitig.insert_mt(current_kmer_RC, -contig_id);
      contigs.set_mt(contig_id, Contig(contig_seq, median(abundances)));

      for(int x= 0; x < 4; x++){
        if(candidates_after[x]){
          kmer = current_kmer_fix+DNA_bases[x];
          if(!startKmer2unitig.find_mt(kmer, idx)){
            int new_unitig_idx = contigs.push_back_mt(Contig(kmer, kmer_abundance_afters[x]));
            startKmer2unitig.insert_mt(kmer, new_unitig_idx);
            work_queue->add_work(1); 
          }
        }
      }
      kmer = current_kmer_RC;
      for(int x = 0; x < 4; x++){
        if(candidates_before[x]){
          kmer[K-1] = DNA_bases[x];
          if(!startKmer2unitig.find_mt(kmer, idx)){
            int new_unitig_idx = contigs.push_back_mt(Contig(kmer, kmer_abundance_befores[x]));
            //startKmer2unitig[kmer_befores_RC[x]] = -new_unitig_idx;
            startKmer2unitig.insert_mt(kmer, new_unitig_idx);
            work_queue->add_work(1);
          }
        }
      }
      break;
    }else if(candidates_after_num==1){ //only one candidate k-mer after
      for(int x = 0; x<4; x++){
        if(candidates_after[x]){
          current_kmer = current_kmer_fix+DNA_bases[x];
          current_kmer_RC = RC_DNAbase(DNA_bases[x])+current_kmer_RC.substr(0, K-1);
          contig_seq += DNA_bases[x];
          abundances.push_back(kmer_abundance_afters[x]);
          NTPC64('T', 'A', options.K, current_kmer_RC_hash, current_kmer_hash);
          NTPC64('T', DNA_bases[x], options.K, current_kmer_hash, current_kmer_RC_hash);
          break;
        }
      }
      continue;
    }else if (nodes_after_num==1){
      current_kmer = current_kmer_fix + DNA_bases[node_after_x];
      current_kmer_RC = RC_DNAbase(DNA_bases[node_after_x]) + current_kmer_RC.substr(0, K-1);
      contig_seq += DNA_bases[node_after_x];
      abundances.push_back(kmer_abundance_afters[node_after_x]);
      NTPC64('T', 'A', options.K, current_kmer_RC_hash, current_kmer_hash);
      NTPC64('T', DNA_bases[node_after_x], options.K, current_kmer_hash, current_kmer_RC_hash);
    }else{ //stop
      startKmer2unitig.insert_mt(current_kmer_RC, -contig_id);
      contigs.set_mt(contig_id, Contig(contig_seq, median(abundances)));
      break;
    }
  }
}
void get_unitig_forward(CQF_mt& cqf, const Params& options, vector_mt<Contig>& contigs, unordered_set_mt& startKmers, WorkQueue* work_queue, const int& contig_id){
  array<bool, 4> candidates_before({false, false, false, false}), candidates_after({false, false, false, false});
  //array<string, 4> kmer_befores, kmer_afters, kmer_befores_RC, kmer_afters_RC;
  array<uint64_t, 4> kmer_abundance_befores, kmer_abundance_afters;
  auto abundance_min = options.kmer_abundance_min;
  auto K = options.K;

  int candidates_before_num, candidates_after_num;
  int nodes_before_num, nodes_after_num;
  uint64_t kmer_hash, kmer_RC_hash, current_kmer_hash, current_kmer_RC_hash;
  string kmer, kmer_RC, current_kmer, current_kmer_RC, current_kmer_fix;
  uint64_t kmer_count, kmer_RC_count;
  int idx;
  
  string contig_seq = contigs.at_mt(contig_id).seq;
  current_kmer = contig_seq.substr(contig_seq.length()-K);
  current_kmer_RC = RC_DNA(current_kmer);

  std::vector<int> abundances;
  abundances.push_back(int(contigs.at_mt(contig_id).median_abundance));

  NTPC64(current_kmer.c_str(), K, current_kmer_hash, current_kmer_RC_hash);
  int node_after_x, node_before_x;//useful only when there is only one node after and without candidates after.
  while(true){
    //NTPC64(current_kmer.c_str(), K, current_kmer_hash, current_kmer_RC_hash);
    current_kmer_fix = current_kmer.substr(1);
    //kmer_afters = kmers_after(current_kmer);
    //kmer_befores = kmers_before(kmer_afters[0]);
    candidates_before = candidates_after = {{false, false, false, false}};
    candidates_before_num = candidates_after_num = 0;
    nodes_before_num = nodes_after_num = 0;
    
    //kmers with current_kmer_fix as prefix
    for(int x = 0; x<4; x++){
      //kmer = kmer_afters[x];
      //kmer_RC = RC_DNA(kmer);
      //kmer_afters_RC[x] = kmer_RC;
      kmer_hash = current_kmer_hash; kmer_RC_hash = current_kmer_RC_hash;
      kmer_hash = NTPC64(current_kmer[0], DNA_bases[x], K, kmer_hash, kmer_RC_hash);
      //kmer_hash = MurmurHash2(kmer);
      //kmer_RC_hash = MurmurHash2(kmer_RC);
      if(cqf.count_key_value_set_traveled(kmer_hash%cqf.qf->metadata->range, kmer_count)){
        if(startKmers.find(current_kmer_fix+DNA_bases[x]) != startKmers.end()){
          nodes_after_num++;
          kmer_abundance_afters[x] = kmer_count;
          node_after_x = x;
        }else if(kmer_count >= abundance_min){
          kmer_abundance_afters[x] = kmer_count;
          candidates_after[x] = true; //possible because of hash collisions
          candidates_after_num ++;
        }
      }else if(kmer_count >= abundance_min){
        kmer_abundance_afters[x] = kmer_count;
        candidates_after[x] = true;
        candidates_after_num ++;
      }
    }
    
    //kmers with RC(current_kmer_fix) as prefix
    NTPC64(current_kmer[0], 'A', options.K, current_kmer_hash, current_kmer_RC_hash);
    kmer = current_kmer_RC;
    for(int x = 0; x<4; x++){
      if(DNA_bases[x] == current_kmer_RC[K-1]){
        continue;
      }
      //kmer = kmer_befores[x];
      //kmer_RC = RC_DNA(kmer);
      //kmer_befores_RC[x] = kmer_RC;
      //kmer_hash = current_kmer_hash; kmer_RC_hash = current_kmer_RC_hash;
      //kmer_hash = NTPC64(current_kmer_RC[0], DNA_bases[x], K, kmer_hash, kmer_RC_hash);
      //kmer = current_kmer_RC;
      kmer[K-1] = DNA_bases[x];
      kmer_hash = current_kmer_hash; kmer_RC_hash = current_kmer_RC_hash;
      kmer_hash = NTPC64('T', DNA_bases[x], options.K, kmer_RC_hash, kmer_hash);
      //kmer_hash = NTPC64(kmer.c_str(), K, kmer_hash, kmer_RC_hash);
      //kmer_hash = MurmurHash2(kmer);
      //kmer_RC_hash = MurmurHash2(kmer_RC);
      if(cqf.count_key_value_set_traveled(kmer_hash%cqf.qf->metadata->range, kmer_count)){
        if(startKmers.find(kmer) != startKmers.end()){
          nodes_before_num++; node_before_x = x;
        }else if(kmer_count >= abundance_min){
          kmer_abundance_befores[x] = kmer_count;
          candidates_before[x] = true;
          candidates_before_num++;
        }
      }else if(kmer_count >= abundance_min){
        kmer_abundance_befores[x] = kmer_count;
        candidates_before[x] = true;
        candidates_before_num++;
      }
    }

    if((nodes_before_num + candidates_before_num) || (nodes_after_num+candidates_after_num)>1){ //no-linear extension
      // if(startKmer2unitig.find_mt(contig_seq.substr(0, params.K), idx)){
      //   //A contig has been constructed by another program from RC way
      //   if(abs(idx) != contig_id){
      //     contigs.set_mt(contig_id, ""); 
      //     break;
      //   }
      // }else{
      //   throw logic_error("Unexpectedly not found start kmer!");
      // }
      startKmers.insert(current_kmer_RC);
      contigs.set_mt(contig_id, Contig(contig_seq, median(abundances)));

      for(int x= 0; x < 4; x++){
        if(candidates_after[x]){
          kmer = current_kmer_fix+DNA_bases[x];
          if(startKmers.find(kmer)==startKmers.end()){
            int new_unitig_idx = contigs.push_back_mt(Contig(kmer, kmer_abundance_afters[x]));
            startKmers.insert(kmer);
            work_queue->add_work(1); 
          }
        }
      }
      kmer = current_kmer_RC;
      for(int x = 0; x < 4; x++){
        if(candidates_before[x]){
          kmer[K-1] = DNA_bases[x];
          if(startKmers.find(kmer) == startKmers.end()){
            int new_unitig_idx = contigs.push_back_mt(Contig(kmer, kmer_abundance_befores[x]));
            //startKmer2unitig[kmer_befores_RC[x]] = -new_unitig_idx;
            startKmers.insert(kmer);
            work_queue->add_work(1);
          }
        }
      }
      break;
    }else if(candidates_after_num==1){ //only one candidate k-mer after
      for(int x = 0; x<4; x++){
        if(candidates_after[x]){
          current_kmer = current_kmer_fix+DNA_bases[x];
          current_kmer_RC = RC_DNAbase(DNA_bases[x])+current_kmer_RC.substr(0, K-1);
          contig_seq += DNA_bases[x];
          abundances.push_back(kmer_abundance_afters[x]);
          NTPC64('T', 'A', options.K, current_kmer_RC_hash, current_kmer_hash);
          NTPC64('T', DNA_bases[x], options.K, current_kmer_hash, current_kmer_RC_hash);
          break;
        }
      }
      continue;
    }else if (nodes_after_num==1){//be careful about circles
      current_kmer = current_kmer_fix + DNA_bases[node_after_x];
      current_kmer_RC = RC_DNAbase(DNA_bases[node_after_x]) + current_kmer_RC.substr(0, K-1);
      contig_seq += DNA_bases[node_after_x];
      abundances.push_back(kmer_abundance_afters[node_after_x]);
      NTPC64('T', 'A', options.K, current_kmer_RC_hash, current_kmer_hash);
      NTPC64('T', DNA_bases[node_after_x], options.K, current_kmer_hash, current_kmer_RC_hash);
    }else{ //stop
      startKmers.insert(current_kmer_RC);
      contigs.set_mt(contig_id, Contig(contig_seq, median(abundances)));
      break;
    }
  }
}

void get_unitig_forward(CQF_mt& cqf, const Params& options, vector_mt<Contig>& contigs, unordered_map_mt<string, int>& startKmer2unitig, WorkQueue* work_queue, Contig& contig){
  int contig_id = 0;//fake contig id

  array<bool, 4> candidates_before({false, false, false, false}), candidates_after({false, false, false, false});
  //array<string, 4> kmer_befores, kmer_afters, kmer_befores_RC, kmer_afters_RC;
  array<uint64_t, 4> kmer_abundance_befores, kmer_abundance_afters;
  auto abundance_min = options.kmer_abundance_min;
  auto K = options.K;

  int candidates_before_num, candidates_after_num;
  int nodes_before_num, nodes_after_num;
  uint64_t kmer_hash, kmer_RC_hash, current_kmer_hash, current_kmer_RC_hash;
  string kmer, kmer_RC, current_kmer, current_kmer_RC, current_kmer_fix;
  uint64_t kmer_count, kmer_RC_count;
  int idx;
  
  string contig_seq = contig.seq;
  current_kmer = contig_seq.substr(contig_seq.length()-K);
  current_kmer_RC = RC_DNA(current_kmer);

  std::vector<int> abundances;
  abundances.push_back(int(contig.median_abundance));

  NTPC64(current_kmer.c_str(), K, current_kmer_hash, current_kmer_RC_hash);
  int node_after_x, node_before_x;//useful only when there is only one node after and without candidates after.
  while(true){
    //NTPC64(current_kmer.c_str(), K, current_kmer_hash, current_kmer_RC_hash);
    current_kmer_fix = current_kmer.substr(1);
    //kmer_afters = kmers_after(current_kmer);
    //kmer_befores = kmers_before(kmer_afters[0]);
    candidates_before = candidates_after = {{false, false, false, false}};
    candidates_before_num = candidates_after_num = 0;
    nodes_before_num = nodes_after_num = 0;
    
    //kmers with current_kmer_fix as prefix
    for(int x = 0; x<4; x++){
      //kmer = kmer_afters[x];
      //kmer_RC = RC_DNA(kmer);
      //kmer_afters_RC[x] = kmer_RC;
      kmer_hash = current_kmer_hash; kmer_RC_hash = current_kmer_RC_hash;
      kmer_hash = NTPC64(current_kmer[0], DNA_bases[x], K, kmer_hash, kmer_RC_hash);
      //kmer_hash = MurmurHash2(kmer);
      //kmer_RC_hash = MurmurHash2(kmer_RC);
      if(cqf.count_key_value_set_traveled(kmer_hash%cqf.qf->metadata->range, kmer_count)){
        if(startKmer2unitig.find_mt(current_kmer_fix+DNA_bases[x], idx)){
          nodes_after_num++;
          kmer_abundance_afters[x] = kmer_count;
          node_after_x = x;
        }else if(kmer_count >= abundance_min){
          kmer_abundance_afters[x] = kmer_count;
          candidates_after[x] = true; //possible because of hash collisions
          candidates_after_num ++;
        }
      }else if(kmer_count >= abundance_min){
        kmer_abundance_afters[x] = kmer_count;
        candidates_after[x] = true;
        candidates_after_num ++;
      }
    }
    
    //kmers with RC(current_kmer_fix) as prefix
    NTPC64(current_kmer[0], 'A', options.K, current_kmer_hash, current_kmer_RC_hash);
    kmer = current_kmer_RC;
    for(int x = 0; x<4; x++){
      if(DNA_bases[x] == current_kmer_RC[K-1]){
        continue;
      }
      //kmer = kmer_befores[x];
      //kmer_RC = RC_DNA(kmer);
      //kmer_befores_RC[x] = kmer_RC;
      //kmer_hash = current_kmer_hash; kmer_RC_hash = current_kmer_RC_hash;
      //kmer_hash = NTPC64(current_kmer_RC[0], DNA_bases[x], K, kmer_hash, kmer_RC_hash);
      //kmer = current_kmer_RC;
      kmer[K-1] = DNA_bases[x];
      kmer_hash = current_kmer_hash; kmer_RC_hash = current_kmer_RC_hash;
      kmer_hash = NTPC64('T', DNA_bases[x], options.K, kmer_RC_hash, kmer_hash);
      //kmer_hash = NTPC64(kmer.c_str(), K, kmer_hash, kmer_RC_hash);
      //kmer_hash = MurmurHash2(kmer);
      //kmer_RC_hash = MurmurHash2(kmer_RC);
      if(cqf.count_key_value_set_traveled(kmer_hash%cqf.qf->metadata->range, kmer_count)){
        if(startKmer2unitig.find_mt(kmer, idx)){
          nodes_before_num++; node_before_x = x;
        }else if(kmer_count >= abundance_min){
          kmer_abundance_befores[x] = kmer_count;
          candidates_before[x] = true;
          candidates_before_num++;
        }
      }else if(kmer_count >= abundance_min){
        kmer_abundance_befores[x] = kmer_count;
        candidates_before[x] = true;
        candidates_before_num++;
      }
    }

    if((nodes_before_num + candidates_before_num) || (nodes_after_num+candidates_after_num)>1){ //no-linear extension
      // if(startKmer2unitig.find_mt(contig_seq.substr(0, params.K), idx)){
      //   //A contig has been constructed by another program from RC way
      //   if(abs(idx) != contig_id){
      //     contigs.set_mt(contig_id, ""); 
      //     break;
      //   }
      // }else{
      //   throw logic_error("Unexpectedly not found start kmer!");
      // }
      startKmer2unitig.insert_mt(current_kmer_RC, -contig_id);
      //contigs.set_mt(contig_id, Contig(contig_seq, median(abundances)));
      contig.seq = contig_seq;
      contig.median_abundance = median(abundances);

      for(int x= 0; x < 4; x++){
        if(candidates_after[x]){
          kmer = current_kmer_fix+DNA_bases[x];
          if(!startKmer2unitig.find_mt(kmer, idx)){
            int new_unitig_idx = contigs.push_back_mt(Contig(kmer, kmer_abundance_afters[x]));
            startKmer2unitig.insert_mt(kmer, new_unitig_idx);
            work_queue->add_work(1); 
          }
        }
      }
      kmer = current_kmer_RC;
      for(int x = 0; x < 4; x++){
        if(candidates_before[x]){
          kmer[K-1] = DNA_bases[x];
          if(!startKmer2unitig.find_mt(kmer, idx)){
            int new_unitig_idx = contigs.push_back_mt(Contig(kmer, kmer_abundance_befores[x]));
            //startKmer2unitig[kmer_befores_RC[x]] = -new_unitig_idx;
            startKmer2unitig.insert_mt(kmer, new_unitig_idx);
            work_queue->add_work(1);
          }
        }
      }
      break;
    }else if(candidates_after_num==1){ //only one candidate k-mer after
      for(int x = 0; x<4; x++){
        if(candidates_after[x]){
          current_kmer = current_kmer_fix+DNA_bases[x];
          current_kmer_RC = RC_DNAbase(DNA_bases[x])+current_kmer_RC.substr(0, K-1);
          contig_seq += DNA_bases[x];
          abundances.push_back(kmer_abundance_afters[x]);
          NTPC64('T', 'A', options.K, current_kmer_RC_hash, current_kmer_hash);
          NTPC64('T', DNA_bases[x], options.K, current_kmer_hash, current_kmer_RC_hash);
          break;
        }
      }
      continue;
    }else if (nodes_after_num==1){
      current_kmer = current_kmer_fix + DNA_bases[node_after_x];
      current_kmer_RC = RC_DNAbase(DNA_bases[node_after_x]) + current_kmer_RC.substr(0, K-1);
      contig_seq += DNA_bases[node_after_x];
      abundances.push_back(kmer_abundance_afters[node_after_x]);
      NTPC64('T', 'A', options.K, current_kmer_RC_hash, current_kmer_hash);
      NTPC64('T', DNA_bases[node_after_x], options.K, current_kmer_hash, current_kmer_RC_hash);
    }else{ //stop
      startKmer2unitig.insert_mt(current_kmer_RC, -contig_id);
      //contigs.set_mt(contig_id, Contig(contig_seq, median(abundances)));
      contig.seq = contig_seq;
      contig.median_abundance = median(abundances);
      break;
    }
  }
}
*/

/*
void get_unitig_forward(CQF_mt& cqf, const Params& options, concurrent_vector<Contig>& contigs, unordered_set_mt& startKmers, WorkQueue* work_queue, concurrent_vector<Contig>::iterator& contigIter){
  array<bool, 4> candidates_before({false, false, false, false}), candidates_after({false, false, false, false});
  //array<string, 4> kmer_befores, kmer_afters, kmer_befores_RC, kmer_afters_RC;
  array<uint64_t, 4> kmer_abundance_befores, kmer_abundance_afters;
  auto abundance_min = options.kmer_abundance_min;
  auto K = options.K;

  int candidates_before_num, candidates_after_num;
  int nodes_before_num, nodes_after_num;
  uint64_t kmer_hash, kmer_RC_hash, current_kmer_hash, current_kmer_RC_hash;
  DNAString kmer, kmer_RC, current_kmer, current_kmer_RC, current_kmer_fix;
  uint64_t kmer_count, kmer_RC_count;
  int idx;
  
  DNAString& contig_seq = contigIter->seq;
  //string contig_seq = contigIter->seq.to_str();
  current_kmer = contig_seq.substr(contig_seq.length()-K);
  current_kmer_RC = current_kmer.get_RC();

  std::vector<int> abundances;
  abundances.push_back(int(contigIter->median_abundance));

  //NTPC64(current_kmer.c_str(), K, current_kmer_hash, current_kmer_RC_hash);
  NTPC64(current_kmer, K, current_kmer_hash, current_kmer_RC_hash);
  int node_after_x, node_before_x;//useful only when there is only one node after and without candidates after.
  while(true){
    //NTPC64(current_kmer.c_str(), K, current_kmer_hash, current_kmer_RC_hash);
    current_kmer_fix = current_kmer.substr(1);
    //kmer_afters = kmers_after(current_kmer);
    //kmer_befores = kmers_before(kmer_afters[0]);
    candidates_before = candidates_after = {{false, false, false, false}};
    candidates_before_num = candidates_after_num = 0;
    nodes_before_num = nodes_after_num = 0;
    
    //kmers with current_kmer_fix as prefix
    for(int x = 0; x<4; x++){
      //kmer = kmer_afters[x];
      //kmer_RC = RC_DNA(kmer);
      //kmer_afters_RC[x] = kmer_RC;
      kmer_hash = current_kmer_hash; kmer_RC_hash = current_kmer_RC_hash;
      kmer_hash = NTPC64(current_kmer[0], DNA_bases[x], K, kmer_hash, kmer_RC_hash);
      //kmer_hash = MurmurHash2(kmer);
      //kmer_RC_hash = MurmurHash2(kmer_RC);
      bool isTraveled=cqf.count_key_value_set_traveled(kmer_hash%cqf.qf->metadata->range, kmer_count);
      if(kmer_count>= abundance_min){//if not traveled, then go
        if(isTraveled && startKmers.find(current_kmer_fix+DNA_bases[x]) != startKmers.end()){
          nodes_after_num++;
          kmer_abundance_afters[x] = kmer_count;
          node_after_x = x;
        //}else if(kmer_count >= abundance_min){
        }else{
          kmer_abundance_afters[x] = kmer_count;
          candidates_after[x] = true; //possible because of hash collisions
          candidates_after_num ++;
        }
      }
    }
    
    //kmers with RC(current_kmer_fix) as prefix
    NTPC64(current_kmer[0], 'A', options.K, current_kmer_hash, current_kmer_RC_hash);
    kmer = current_kmer_RC;
    for(int x = 0; x<4; x++){
      if(DNA_bases[x] == current_kmer_RC[K-1]){
        continue;
      }
      //kmer = kmer_befores[x];
      //kmer_RC = RC_DNA(kmer);
      //kmer_befores_RC[x] = kmer_RC;
      //kmer_hash = current_kmer_hash; kmer_RC_hash = current_kmer_RC_hash;
      //kmer_hash = NTPC64(current_kmer_RC[0], DNA_bases[x], K, kmer_hash, kmer_RC_hash);
      //kmer = current_kmer_RC;
      kmer.replace(K-1, DNA_bases[x]);
      kmer_hash = current_kmer_hash; kmer_RC_hash = current_kmer_RC_hash;
      kmer_hash = NTPC64('T', DNA_bases[x], options.K, kmer_RC_hash, kmer_hash);
      //kmer_hash = NTPC64(kmer.c_str(), K, kmer_hash, kmer_RC_hash);
      //kmer_hash = MurmurHash2(kmer);
      //kmer_RC_hash = MurmurHash2(kmer_RC);
      bool isTraveled = cqf.count_key_value_set_traveled(kmer_hash%cqf.qf->metadata->range, kmer_count);
      if(kmer_count >= abundance_min){
        if(isTraveled && startKmers.find(kmer) != startKmers.end()){
          nodes_before_num++; node_before_x = x;
        //}else if(kmer_count >= abundance_min){
        }else{
          kmer_abundance_befores[x] = kmer_count;
          candidates_before[x] = true;
          candidates_before_num++;
        }
      }
    }

    if((nodes_before_num + candidates_before_num) || (nodes_after_num+candidates_after_num)>1){ //no-linear extension
      // if(startKmer2unitig.find_mt(contig_seq.substr(0, params.K), idx)){
      //   //A contig has been constructed by another program from RC way
      //   if(abs(idx) != contig_id){
      //     contigs.set_mt(contig_id, ""); 
      //     break;
      //   }
      // }else{
      //   throw logic_error("Unexpectedly not found start kmer!");
      // }
      startKmers.insert(current_kmer_RC);
      //contigs.set_mt(contig_id, Contig(contig_seq, median(abundances)));
      //contigIter->seq = contig_seq;
      contigIter->median_abundance = median(abundances);

      for(int x= 0; x < 4; x++){
        if(candidates_after[x]){
          kmer = current_kmer_fix+DNA_bases[x];
          if(startKmers.find(kmer)==startKmers.end()){
            startKmers.insert(kmer);
            auto it = contigs.push_back(Contig(kmer, kmer_abundance_afters[x]));
            work_queue->add_work(it); 
          }
        }
      }
      kmer = current_kmer_RC;
      for(int x = 0; x < 4; x++){
        if(candidates_before[x]){
          kmer.replace(K-1, DNA_bases[x]);
          if(startKmers.find(kmer) == startKmers.end()){
            startKmers.insert(kmer);
            auto it = contigs.push_back(Contig(kmer, kmer_abundance_befores[x]));
            //startKmer2unitig[kmer_befores_RC[x]] = -new_unitig_idx;
            work_queue->add_work(it);
          }
        }
      }
      break;
    }else if(candidates_after_num==1){ //only one candidate k-mer after
      for(int x = 0; x<4; x++){
        if(candidates_after[x]){
          current_kmer = current_kmer_fix+DNA_bases[x];
          current_kmer_RC = RC_DNAbase(DNA_bases[x])+current_kmer_RC.substr(0, K-1);
          contig_seq += DNA_bases[x];
          abundances.push_back(kmer_abundance_afters[x]);
          NTPC64('T', 'A', options.K, current_kmer_RC_hash, current_kmer_hash);
          NTPC64('T', DNA_bases[x], options.K, current_kmer_hash, current_kmer_RC_hash);
          break;
        }
      }
      continue;
    }else if (nodes_after_num==1){//be careful about circles
      current_kmer = current_kmer_fix + DNA_bases[node_after_x];
      current_kmer_RC = RC_DNAbase(DNA_bases[node_after_x]) + current_kmer_RC.substr(0, K-1);
      contig_seq += DNA_bases[node_after_x];
      abundances.push_back(kmer_abundance_afters[node_after_x]);
      NTPC64('T', 'A', options.K, current_kmer_RC_hash, current_kmer_hash);
      NTPC64('T', DNA_bases[node_after_x], options.K, current_kmer_hash, current_kmer_RC_hash);
    }else{ //stop
      startKmers.insert(current_kmer_RC);
      //contigIter->seq = contig_seq;
      contigIter->median_abundance = median(abundances);
      //contigs.set_mt(contig_id, Contig(contig_seq, median(abundances)));
      break;
    }
  }
}
*/

//True if it was inserted. 
//if dnaStr exists, insert only when idx is smaller than existing one; otherwise simply insert the new record. 
bool insert_or_replace(hash_map_mt& startKmer2unitig, const DNAString& dnastr, const size_t idx){
  hash_map_mt::accessor access;
  if(startKmer2unitig.insert(access, dnastr) || access->second >= idx){// insert: True if new pair was inserted; false if key was already in the map.
    access->second = idx; access.release();
    return true;
  } access.release();
  return false;
}

//TODO
void get_unitig_forward(CQF_mt& cqf, const Params& options, concurrent_vector<Contig>& contigs, hash_map_mt& startKmer2unitig, WorkQueue* work_queue, concurrent_vector<Contig>::iterator& contigIter){
  array<bool, 4> candidates_before({false, false, false, false}), candidates_after({false, false, false, false});
  //array<string, 4> kmer_befores, kmer_afters, kmer_befores_RC, kmer_afters_RC;
  array<uint64_t, 4> kmer_abundance_befores, kmer_abundance_afters;
  auto abundance_min = options.kmer_abundance_min;
  auto K = options.K;

  int candidates_before_num, candidates_after_num;
  int nodes_before_num, nodes_after_num;
  uint64_t kmer_hash, kmer_RC_hash, current_kmer_hash, current_kmer_RC_hash;
  DNAString kmer, kmer_RC, current_kmer, current_kmer_RC, current_kmer_fix, first_kmer;
  uint64_t kmer_count, kmer_RC_count;
  int idx;
  
  DNAString& contig_seq = contigIter->seq;
  first_kmer = contig_seq.substr(0, options.K);
  //string contig_seq = contigIter->seq.to_str();
  current_kmer = contig_seq.substr(contig_seq.length()-K);
  current_kmer_RC = current_kmer.get_RC();

  std::vector<int> abundances(contigIter->seq.length()-K+1, int(contigIter->median_abundance));

  hash_map_mt::const_accessor const_access;
  //NTPC64(current_kmer.c_str(), K, current_kmer_hash, current_kmer_RC_hash);
  NTPC64(current_kmer, K, current_kmer_hash, current_kmer_RC_hash);
  int node_after_x, node_before_x;//useful only when there is only one node after and without candidates after.
  while(true){
    //NTPC64(current_kmer.c_str(), K, current_kmer_hash, current_kmer_RC_hash);
    current_kmer_fix = current_kmer.substr(1);
    //kmer_afters = kmers_after(current_kmer);
    //kmer_befores = kmers_before(kmer_afters[0]);
    candidates_before = candidates_after = {{false, false, false, false}};
    candidates_before_num = candidates_after_num = 0;
    nodes_before_num = nodes_after_num = 0;
    
    //kmers with current_kmer_fix as prefix
    for(int x = 0; x<4; x++){
      //kmer = kmer_afters[x];
      //kmer_RC = RC_DNA(kmer);
      //kmer_afters_RC[x] = kmer_RC;
      kmer_hash = current_kmer_hash; kmer_RC_hash = current_kmer_RC_hash;
      kmer_hash = NTPC64(current_kmer[0], DNA_bases[x], K, kmer_hash, kmer_RC_hash);
      //kmer_hash = MurmurHash2(kmer);
      //kmer_RC_hash = MurmurHash2(kmer_RC);

      bool isTraveled=cqf.count_key_value_set_traveled(kmer_hash%cqf.qf->metadata->range, kmer_count);
      if(kmer_count>= abundance_min){//if not traveled, then go
        if(isTraveled && startKmer2unitig.find(const_access, current_kmer_fix+DNA_bases[x])){
          nodes_after_num++;
          kmer_abundance_afters[x] = kmer_count;
          node_after_x = x;
          const_access.release();
        //}else if(kmer_count >= abundance_min){
        }else{
          kmer_abundance_afters[x] = kmer_count;
          candidates_after[x] = true; //possible because of hash collisions
          candidates_after_num ++;
        }
      }
    }
    
    //kmers with RC(current_kmer_fix) as prefix
    NTPC64(current_kmer[0], 'A', options.K, current_kmer_hash, current_kmer_RC_hash);
    kmer = current_kmer_RC;
    for(int x = 0; x<4; x++){
      if(DNA_bases[x] == current_kmer_RC[K-1]){
        continue;
      }
      //kmer = kmer_befores[x];
      //kmer_RC = RC_DNA(kmer);
      //kmer_befores_RC[x] = kmer_RC;
      //kmer_hash = current_kmer_hash; kmer_RC_hash = current_kmer_RC_hash;
      //kmer_hash = NTPC64(current_kmer_RC[0], DNA_bases[x], K, kmer_hash, kmer_RC_hash);
      //kmer = current_kmer_RC;
      kmer.replace(K-1, DNA_bases[x]);
      kmer_hash = current_kmer_hash; kmer_RC_hash = current_kmer_RC_hash;
      kmer_hash = NTPC64('T', DNA_bases[x], options.K, kmer_RC_hash, kmer_hash);
      //kmer_hash = NTPC64(kmer.c_str(), K, kmer_hash, kmer_RC_hash);
      //kmer_hash = MurmurHash2(kmer);
      //kmer_RC_hash = MurmurHash2(kmer_RC);
      bool isTraveled = cqf.count_key_value_set_traveled(kmer_hash%cqf.qf->metadata->range, kmer_count);
      if(kmer_count >= abundance_min){
        if(isTraveled && startKmer2unitig.find(const_access, kmer)){
          nodes_before_num++; node_before_x = x;
        //}else if(kmer_count >= abundance_min){
          const_access.release();
        }else{
          kmer_abundance_befores[x] = kmer_count;
          candidates_before[x] = true;
          candidates_before_num++;
        }
      }
    }

    if((nodes_before_num + candidates_before_num) || (nodes_after_num+candidates_after_num)>1){ //no-linear extension
      // if(startKmer2unitig.find_mt(contig_seq.substr(0, params.K), idx)){
      //   //A contig has been constructed by another program from RC way
      //   if(abs(idx) != contig_id){
      //     contigs.set_mt(contig_id, ""); 
      //     break;
      //   }
      // }else{
      //   throw logic_error("Unexpectedly not found start kmer!");
      // }
      if(!insert_or_replace(startKmer2unitig, DNAString(current_kmer_RC), (contigIter-contigs.begin()))){
        contigIter->clear();
        break;
      }

      //startKmers.insert(current_kmer_RC);
      contigIter->median_abundance = median(abundances);

      for(int x= 0; x < 4; x++){
        if(candidates_after[x]){
          kmer = current_kmer_fix+DNA_bases[x];

          hash_map_mt::accessor access;
          if(startKmer2unitig.insert(access, kmer)){
            auto it = contigs.push_back(Contig(kmer, kmer_abundance_afters[x]));
            access->second = it - contigs.begin();
            work_queue->add_work(it);
          }
          access.release();
        }
      }
      kmer = current_kmer_RC;
      for(int x = 0; x < 4; x++){
        if(candidates_before[x]){
          kmer.replace(K-1, DNA_bases[x]);
          hash_map_mt::accessor access;
          if(startKmer2unitig.insert(access, kmer)){
            auto it = contigs.push_back(Contig(kmer, kmer_abundance_befores[x]));
            access->second = it - contigs.begin();
            work_queue->add_work(it);
          }
          access.release();
        }
      }
      break;
    }else if(candidates_after_num==1){ //only one candidate k-mer after #Warning: circles
      int x = 0;
      for(x = 0; x<4; x++){
        if(candidates_after[x]){
          break;
        }
      }
      current_kmer = current_kmer_fix+DNA_bases[x];
          
      if(current_kmer == first_kmer){//it is an pure circle
        int contig_id = contigIter-contigs.begin();
        if(!insert_or_replace(startKmer2unitig, DNAString(first_kmer), contig_id) || !insert_or_replace(startKmer2unitig, DNAString(current_kmer_RC), contig_id)){
          contigIter->clear();
        }else{
          contigIter->median_abundance = median(abundances);
        }
        break;
      }else{
        current_kmer_RC = RC_DNAbase(DNA_bases[x])+current_kmer_RC.substr(0, K-1);
        contig_seq += DNA_bases[x];
        abundances.push_back(kmer_abundance_afters[x]);
        NTPC64('T', 'A', options.K, current_kmer_RC_hash, current_kmer_hash);
        NTPC64('T', DNA_bases[x], options.K, current_kmer_hash, current_kmer_RC_hash);
        continue;
      }
    }else if (nodes_after_num==1){//be careful about circles
      // current_kmer = current_kmer_fix + DNA_bases[node_after_x];
      // current_kmer_RC = RC_DNAbase(DNA_bases[node_after_x]) + current_kmer_RC.substr(0, K-1);
      // contig_seq += DNA_bases[node_after_x];
      // abundances.push_back(kmer_abundance_afters[node_after_x]);
      // NTPC64('T', 'A', options.K, current_kmer_RC_hash, current_kmer_hash);
      // NTPC64('T', DNA_bases[node_after_x], options.K, current_kmer_hash, current_kmer_RC_hash);
      if(!insert_or_replace(startKmer2unitig, DNAString(current_kmer_RC), (contigIter-contigs.begin()))){
        contigIter->clear();
      }else{
        contigIter->median_abundance = median(abundances);
      }
      break;
    }else{ //stop
      if(!insert_or_replace(startKmer2unitig, DNAString(current_kmer_RC), (contigIter-contigs.begin()))){
        contigIter->clear();
      }else{
        contigIter->median_abundance = median(abundances);
      }
      //startKmers.insert(current_kmer_RC);
      //contigIter->seq = contig_seq;
      
      //contigs.set_mt(contig_id, Contig(contig_seq, median(abundances)));
      break;
    }
  }
}

/*
void get_unitig_forward(CQF_mt& cqf, const Params& options, concurrent_vector<Contig>& contigs, unordered_map_mt<string, int>& startKmer2unitig, WorkQueue* work_queue, Contig& contig){
  int contig_id = 0;//fake contig id

  array<bool, 4> candidates_before({false, false, false, false}), candidates_after({false, false, false, false});
  //array<string, 4> kmer_befores, kmer_afters, kmer_befores_RC, kmer_afters_RC;
  array<uint64_t, 4> kmer_abundance_befores, kmer_abundance_afters;
  auto abundance_min = options.kmer_abundance_min;
  auto K = options.K;

  int candidates_before_num, candidates_after_num;
  int nodes_before_num, nodes_after_num;
  uint64_t kmer_hash, kmer_RC_hash, current_kmer_hash, current_kmer_RC_hash;
  string kmer, kmer_RC, current_kmer, current_kmer_RC, current_kmer_fix;
  uint64_t kmer_count, kmer_RC_count;
  int idx;
  
  string contig_seq = contig.seq;
  current_kmer = contig_seq.substr(contig_seq.length()-K);
  current_kmer_RC = RC_DNA(current_kmer);

  std::vector<int> abundances;
  abundances.push_back(int(contig.median_abundance));

  NTPC64(current_kmer.c_str(), K, current_kmer_hash, current_kmer_RC_hash);
  int node_after_x, node_before_x;//useful only when there is only one node after and without candidates after.
  while(true){
    //NTPC64(current_kmer.c_str(), K, current_kmer_hash, current_kmer_RC_hash);
    current_kmer_fix = current_kmer.substr(1);
    //kmer_afters = kmers_after(current_kmer);
    //kmer_befores = kmers_before(kmer_afters[0]);
    candidates_before = candidates_after = {{false, false, false, false}};
    candidates_before_num = candidates_after_num = 0;
    nodes_before_num = nodes_after_num = 0;
    
    //kmers with current_kmer_fix as prefix
    for(int x = 0; x<4; x++){
      //kmer = kmer_afters[x];
      //kmer_RC = RC_DNA(kmer);
      //kmer_afters_RC[x] = kmer_RC;
      kmer_hash = current_kmer_hash; kmer_RC_hash = current_kmer_RC_hash;
      kmer_hash = NTPC64(current_kmer[0], DNA_bases[x], K, kmer_hash, kmer_RC_hash);
      //kmer_hash = MurmurHash2(kmer);
      //kmer_RC_hash = MurmurHash2(kmer_RC);
      if(cqf.count_key_value_set_traveled(kmer_hash%cqf.qf->metadata->range, kmer_count)){
        if(startKmer2unitig.find_mt(current_kmer_fix+DNA_bases[x], idx)){
          nodes_after_num++;
          kmer_abundance_afters[x] = kmer_count;
          node_after_x = x;
        }else if(kmer_count >= abundance_min){
          kmer_abundance_afters[x] = kmer_count;
          candidates_after[x] = true; //possible because of hash collisions
          candidates_after_num ++;
        }
      }else if(kmer_count >= abundance_min){
        kmer_abundance_afters[x] = kmer_count;
        candidates_after[x] = true;
        candidates_after_num ++;
      }
    }
    
    //kmers with RC(current_kmer_fix) as prefix
    NTPC64(current_kmer[0], 'A', options.K, current_kmer_hash, current_kmer_RC_hash);
    kmer = current_kmer_RC;
    for(int x = 0; x<4; x++){
      if(DNA_bases[x] == current_kmer_RC[K-1]){
        continue;
      }
      //kmer = kmer_befores[x];
      //kmer_RC = RC_DNA(kmer);
      //kmer_befores_RC[x] = kmer_RC;
      //kmer_hash = current_kmer_hash; kmer_RC_hash = current_kmer_RC_hash;
      //kmer_hash = NTPC64(current_kmer_RC[0], DNA_bases[x], K, kmer_hash, kmer_RC_hash);
      //kmer = current_kmer_RC;
      kmer[K-1] = DNA_bases[x];
      kmer_hash = current_kmer_hash; kmer_RC_hash = current_kmer_RC_hash;
      kmer_hash = NTPC64('T', DNA_bases[x], options.K, kmer_RC_hash, kmer_hash);
      //kmer_hash = NTPC64(kmer.c_str(), K, kmer_hash, kmer_RC_hash);
      //kmer_hash = MurmurHash2(kmer);
      //kmer_RC_hash = MurmurHash2(kmer_RC);
      if(cqf.count_key_value_set_traveled(kmer_hash%cqf.qf->metadata->range, kmer_count)){
        if(startKmer2unitig.find_mt(kmer, idx)){
          nodes_before_num++; node_before_x = x;
        }else if(kmer_count >= abundance_min){
          kmer_abundance_befores[x] = kmer_count;
          candidates_before[x] = true;
          candidates_before_num++;
        }
      }else if(kmer_count >= abundance_min){
        kmer_abundance_befores[x] = kmer_count;
        candidates_before[x] = true;
        candidates_before_num++;
      }
    }

    if((nodes_before_num + candidates_before_num) || (nodes_after_num+candidates_after_num)>1){ //no-linear extension
      // if(startKmer2unitig.find_mt(contig_seq.substr(0, params.K), idx)){
      //   //A contig has been constructed by another program from RC way
      //   if(abs(idx) != contig_id){
      //     contigs.set_mt(contig_id, ""); 
      //     break;
      //   }
      // }else{
      //   throw logic_error("Unexpectedly not found start kmer!");
      // }
      startKmer2unitig.insert_mt(current_kmer_RC, -contig_id);
      //contigs.set_mt(contig_id, Contig(contig_seq, median(abundances)));
      contig.seq = contig_seq;
      contig.median_abundance = median(abundances);

      for(int x= 0; x < 4; x++){
        if(candidates_after[x]){
          kmer = current_kmer_fix+DNA_bases[x];
          if(!startKmer2unitig.find_mt(kmer, idx)){
            int new_unitig_idx = contigs.push_back_mt(Contig(kmer, kmer_abundance_afters[x]));
            startKmer2unitig.insert_mt(kmer, new_unitig_idx);
            work_queue->add_work(1); 
          }
        }
      }
      kmer = current_kmer_RC;
      for(int x = 0; x < 4; x++){
        if(candidates_before[x]){
          kmer[K-1] = DNA_bases[x];
          if(!startKmer2unitig.find_mt(kmer, idx)){
            int new_unitig_idx = contigs.push_back_mt(Contig(kmer, kmer_abundance_befores[x]));
            //startKmer2unitig[kmer_befores_RC[x]] = -new_unitig_idx;
            startKmer2unitig.insert_mt(kmer, new_unitig_idx);
            work_queue->add_work(1);
          }
        }
      }
      break;
    }else if(candidates_after_num==1){ //only one candidate k-mer after
      for(int x = 0; x<4; x++){
        if(candidates_after[x]){
          current_kmer = current_kmer_fix+DNA_bases[x];
          current_kmer_RC = RC_DNAbase(DNA_bases[x])+current_kmer_RC.substr(0, K-1);
          contig_seq += DNA_bases[x];
          abundances.push_back(kmer_abundance_afters[x]);
          NTPC64('T', 'A', options.K, current_kmer_RC_hash, current_kmer_hash);
          NTPC64('T', DNA_bases[x], options.K, current_kmer_hash, current_kmer_RC_hash);
          break;
        }
      }
      continue;
    }else if (nodes_after_num==1){
      current_kmer = current_kmer_fix + DNA_bases[node_after_x];
      current_kmer_RC = RC_DNAbase(DNA_bases[node_after_x]) + current_kmer_RC.substr(0, K-1);
      contig_seq += DNA_bases[node_after_x];
      abundances.push_back(kmer_abundance_afters[node_after_x]);
      NTPC64('T', 'A', options.K, current_kmer_RC_hash, current_kmer_hash);
      NTPC64('T', DNA_bases[node_after_x], options.K, current_kmer_hash, current_kmer_RC_hash);
    }else{ //stop
      startKmer2unitig.insert_mt(current_kmer_RC, -contig_id);
      //contigs.set_mt(contig_id, Contig(contig_seq, median(abundances)));
      contig.seq = contig_seq;
      contig.median_abundance = median(abundances);
      break;
    }
  }
}
*/

/*
void get_unitig_backward(CQF_mt& cqf, const Params& params, vector_mt<string>& contigs, unordered_map_mt<string, int>& startKmer2unitig, WorkQueue* work_queue, int contig_id){
  array<bool, 4> candidates_before, candidates_after;
  array<string, 4> kmer_befores, kmer_afters, kmer_befores_RC, kmer_afters_RC;
  int candidates_before_num, candidates_after_num;
  int nodes_before_num, nodes_after_num;
  uint64_t kmer_hash, kmer_RC_hash;
  string kmer, kmer_RC, current_kmer;
  uint64_t kmer_count, kmer_RC_count;
  int idx;

  string contig_seq = contigs.at_mt(contig_id); 
  current_kmer = contig_seq.substr(0, params.K);
  while(true){
    candidates_before = candidates_after = {{false, false, false, false}};
    kmer_befores = kmers_before(current_kmer);
    kmer_afters = kmers_after(kmer_befores[0]);
    candidates_before_num = candidates_after_num = 0;
    nodes_before_num = nodes_after_num = 0;

    for(int x = 0; x<4; x++){
      if(kmer_afters[x]==current_kmer){
        continue;
      }
      kmer = kmer_afters[x];
      kmer_RC = RC_DNA(kmer);
      kmer_afters_RC[x] = kmer_RC;

      NTPC64(kmer.c_str(), params.K, kmer_hash, kmer_RC_hash);
      //kmer_hash = MurmurHash2(kmer);
      //kmer_RC_hash = MurmurHash2(kmer_RC);
      if(cqf.count_key_value_set_traveled(kmer_hash%cqf.qf->metadata->range, kmer_count) && cqf.count_key_value_set_traveled(kmer_RC_hash%cqf.qf->metadata->range, kmer_RC_count)){
        if(startKmer2unitig.find_mt(kmer, idx)){
          nodes_after_num++;
        }else if(kmer_count >= params.CONTIG_min_cov && kmer_RC_count >= params.CONTIG_min_cov){
          candidates_after[x] = true;
          candidates_after_num ++;
        }
      }else if(kmer_count >= params.CONTIG_min_cov && kmer_RC_count >= params.CONTIG_min_cov){
        candidates_after[x] = true;
        candidates_after_num ++;        
      }
    }  
    for(int x = 0; x<4; x++){
      kmer = kmer_befores[x];
      kmer_RC = RC_DNA(kmer);
      kmer_befores_RC[x] = kmer_RC;

      NTPC64(kmer.c_str(), params.K, kmer_hash, kmer_RC_hash);
      //kmer_hash = MurmurHash2(kmer);
      //kmer_RC_hash  = MurmurHash2(kmer_RC);
      if(cqf.count_key_value_set_traveled(kmer_hash%cqf.qf->metadata->range, kmer_count) && cqf.count_key_value_set_traveled(kmer_RC_hash%cqf.qf->metadata->range, kmer_RC_count)){
        if(startKmer2unitig.find_mt(kmer_RC, idx)){
          nodes_before_num++; 
        }else if (kmer_count >= params.CONTIG_min_cov && kmer_RC_count >= params.CONTIG_min_cov){
          candidates_before[x] = true;
          candidates_before_num++;
        }
      }else if(kmer_count >= params.CONTIG_min_cov && kmer_RC_count >= params.CONTIG_min_cov){
        candidates_before[x] = true;
        candidates_before_num++;
      }
    }
    if((nodes_before_num + candidates_before_num)>1 || (nodes_after_num+candidates_after_num)){
      if(startKmer2unitig.find_mt(RC_DNA(contig_seq.substr(contig_seq.length()-params.K)), idx)){
        //A contig has been constructed by another program from RC way
        if(abs(idx) != contig_id){
          contigs.set_mt(contig_id, ""); 
          break;
        }
      }else{
        throw logic_error("Unexpectedly not found start kmer!");
      }
      startKmer2unitig.insert_mt(current_kmer, contig_id);
      
      for(int x= 0; x < 4; x++){
        if(candidates_after[x]){
          if(!startKmer2unitig.find_mt(kmer_afters[x], idx)){
            int new_unitig_idx = contigs.push_back_mt(kmer_afters[x]);
            startKmer2unitig.insert_mt(kmer_afters[x], new_unitig_idx);
            work_queue->add_work(1); 
          }
        }
      }
      for(int x = 0; x < 4; x++){
        if(candidates_before[x]){
          if(!startKmer2unitig.find_mt(kmer_befores_RC[x], idx)){
            int new_unitig_idx = contigs.push_back_mt(kmer_befores_RC[x]);
            startKmer2unitig[kmer_befores_RC[x]] = -new_unitig_idx;
            work_queue->add_work(1);
          }
        }
      }     
    }else if(candidates_before_num==1){
      for(int x = 0; x<4; x++){
        if(candidates_before[x]){
          current_kmer = kmer_befores[x];
          contig_seq = current_kmer.front()+contig_seq;
          break;
        }
      }
      continue;
    }else{
      startKmer2unitig.insert_mt(current_kmer, contig_id);
    }
    break;
  }
}
*/

/*
void get_unitig_backward(CQF_mt& cqf, const Params& options, vector_mt<Contig>& contigs, unordered_map_mt<string, int>& startKmer2unitig, WorkQueue* work_queue, int contig_id){
  auto abundance_min = options.kmer_abundance_min;
  auto K = options.K;

  array<bool, 4> candidates_before, candidates_after;
  //array<string, 4> kmer_befores, kmer_afters, kmer_befores_RC, kmer_afters_RC;
  int candidates_before_num, candidates_after_num;
  int nodes_before_num, nodes_after_num;
  uint64_t kmer_hash, kmer_RC_hash;
  string kmer, kmer_RC, current_kmer;
  uint64_t kmer_count, kmer_RC_count;                                                       
  int idx;

  string contig_seq = contigs.at_mt(contig_id); 
  current_kmer = contig_seq.substr(0, K);
  current_kmer_RC = RC_DNA(current_kmer);

  while(true){


    candidates_before = candidates_after = {{false, false, false, false}};
    //kmer_befores = kmers_before(current_kmer);
    //kmer_afters = kmers_after(kmer_befores[0]);
    candidates_before_num = candidates_after_num = 0;
    nodes_before_num = nodes_after_num = 0;

    for(int x = 0; x<4; x++){
      if(kmer_afters[x]==current_kmer){
        continue;
      }
      kmer = kmer_afters[x];
      kmer_RC = RC_DNA(kmer);
      kmer_afters_RC[x] = kmer_RC;

      NTPC64(kmer.c_str(), K, kmer_hash, kmer_RC_hash);
      if(kmer_hash > kmer_RC_hash){
        kmer_hash = kmer_RC_hash;
      }
      //kmer_hash = MurmurHash2(kmer);
      //kmer_RC_hash = MurmurHash2(kmer_RC);
      if(cqf.count_key_value_set_traveled(kmer_hash%cqf.qf->metadata->range, kmer_count)){
        if(startKmer2unitig.find_mt(kmer, idx)){
          nodes_after_num++;
        }else if(kmer_count >= abundance_min){
          candidates_after[x] = true;
          candidates_after_num ++;
        }
      }else if(kmer_count >= abundance_min){
        candidates_after[x] = true;
        candidates_after_num ++;        
      }
    }  
    for(int x = 0; x<4; x++){
      kmer = kmer_befores[x];
      kmer_RC = RC_DNA(kmer);
      kmer_befores_RC[x] = kmer_RC;

      NTPC64(kmer.c_str(), K, kmer_hash, kmer_RC_hash);
      if(kmer_hash > kmer_RC_hash){
        kmer_hash = kmer_RC_hash;
      }
      //kmer_hash = MurmurHash2(kmer);
      //kmer_RC_hash  = MurmurHash2(kmer_RC);
      if(cqf.count_key_value_set_traveled(kmer_hash%cqf.qf->metadata->range, kmer_count)){
        if(startKmer2unitig.find_mt(kmer_RC, idx)){
          nodes_before_num++; 
        }else if (kmer_count >= abundance_min){
          candidates_before[x] = true;
          candidates_before_num++;
        }
      }else if(kmer_count >= abundance_min){
        candidates_before[x] = true;
        candidates_before_num++;
      }
    }
    if((nodes_before_num + candidates_before_num)>1 || (nodes_after_num+candidates_after_num)){
      // if(startKmer2unitig.find_mt(RC_DNA(contig_seq.substr(contig_seq.length()-params.K)), idx)){
      //   //A contig has been constructed by another program from RC way
      //   if(abs(idx) != contig_id){
      //     contigs.set_mt(contig_id, ""); 
      //     break;
      //   }
      // }else{
      //   throw logic_error("Unexpectedly not found start kmer!");
      // }
      startKmer2unitig.insert_mt(current_kmer, contig_id);
      contigs.set_mt(contig_id, contig_seq);
      
      for(int x= 0; x < 4; x++){
        if(candidates_after[x]){
          if(!startKmer2unitig.find_mt(kmer_afters[x], idx)){
            int new_unitig_idx = contigs.push_back_mt(kmer_afters[x]);
            startKmer2unitig.insert_mt(kmer_afters[x], new_unitig_idx);
            work_queue->add_work(1); 
          }
        }
      }
      for(int x = 0; x < 4; x++){
        if(candidates_before[x]){
          if(!startKmer2unitig.find_mt(kmer_befores_RC[x], idx)){
            int new_unitig_idx = contigs.push_back_mt(kmer_befores_RC[x]);
            startKmer2unitig[kmer_befores_RC[x]] = -new_unitig_idx;
            work_queue->add_work(1);
          }
        }
      }
      break;     
    }else if(candidates_before_num==1){
      for(int x = 0; x<4; x++){
        if(candidates_before[x]){
          current_kmer = kmer_befores[x];
          contig_seq = current_kmer.front()+contig_seq;
          break;
        }
      }
      continue;
    }else{
      startKmer2unitig.insert_mt(current_kmer, contig_id);
      contigs.set_mt(contig_id, contig_seq);
      break;
    }
  }
}
*/

#if 0
void find_unitigs_hash_mt(CQF_mt& cqf, const Sequences& seqs, const Params& params, vector<string>& contigs){
  string unitig, kmer, kmer_RC;
  uint64_t count, kmer_hash, kmer_RC_hash;

  string seq_RC;

  vector_mt<string> unitigs;
  vector<Unitig_node> unitig_nodes;
  unordered_map_mt<string, int> startKmer2unitig(1ULL<<30);//, endKmer2unitig;

  int counter = 0;
  for(auto seq:seqs){
    if(seq.length()<params.K){
      continue;
    }
    //seq_RC = RC_DNA(seq);
    int seq_len = seq.length();
    int step = seq_len/3;
    for(int x = 0; x<=seq_len-params.K; x += step){
      unitigs.clear();
      unitig_nodes.clear();
      startKmer2unitig.clear(); endKmer2unitig.clear();   

      kmer = seq.substr(x, params.K);
      kmer_RC = RC_DNA(kmer);//seq_RC.substr(seq_len-params.K-x, params.K);
      kmer_hash = MurmurHash2(kmer);
      kmer_RC_hash = MurmurHash2(kmer_RC);
      
      //if(kmer=="CAGAACTCAAGGAAACACTTATGTTTACAGGT"){
      //  cout<<"Good!"<<endl;
      //}else{
      //  continue;
      //}

      if(cqf.is_traveled(kmer_hash)){
        continue;
      }else if(cqf.count(kmer_hash) < params.CONTIG_min_cov || cqf.count(kmer_RC_hash) < params.CONTIG_min_cov){
        cqf.set_traveled(kmer_hash);
        cqf.set_traveled(kmer_RC_hash);
        continue;
      }
      counter++;
      
      /*
      if(counter==14){
        fasta_write(contigs, params.OUT_CONTIG_min_len, dir+"counter13.fa");
      }
      if(contigs.size()>26046){
        cout<<counter<<endl;
        contig_summary(contigs);
        fasta_write(contigs, params.OUT_CONTIG_min_len, dir+"tmp.fa");
      }
      */

      cqf.set_traveled(kmer_hash);
      cqf.set_traveled(kmer_RC_hash);   
    
      startKmer2unitig[kmer] = 1;
      endKmer2unitig[kmer_RC] = -1;
      unitigs.resize(2);
      unitigs[1] = kmer;
      unitig_nodes.resize(2);

      get_unitig_forward1(cqf, params, unitigs, unitig_nodes, startKmer2unitig, endKmer2unitig, 1, kmer);
      get_unitig_backward1(cqf, params, unitigs, unitig_nodes, startKmer2unitig, endKmer2unitig, 1, kmer);
      int current_unitig_idx = 1, unitigs_before_num, unitigs_after_num;
      while(current_unitig_idx < unitigs.size()-1){
        current_unitig_idx++;
        unitigs_before_num = unitig_nodes[current_unitig_idx].unitigs_before.size();
        unitigs_after_num = unitig_nodes[current_unitig_idx].unitigs_after.size();
        
        if(unitigs_before_num>0 && unitigs_after_num==0){
          //check for existence of the same unitigs
          if(startKmer2unitig[unitigs[current_unitig_idx]]!=current_unitig_idx){
            for(auto ele:unitig_nodes[current_unitig_idx].unitigs_before){
              if(ele>0){
                unitig_nodes[ele].unitigs_after.erase(current_unitig_idx);
              }else{
                unitig_nodes[-ele].unitigs_before.erase(-current_unitig_idx); 
              }
            }
            unitig_nodes[current_unitig_idx].unitigs_before.clear();
            unitigs[current_unitig_idx] = "";
            continue;
          }
          get_unitig_forward1(cqf, params, unitigs, unitig_nodes, startKmer2unitig, endKmer2unitig, current_unitig_idx, unitigs[current_unitig_idx]);
        }else if(unitigs_before_num==0 && unitigs_after_num>0){
          //check for existence of the same unitigs
          if(endKmer2unitig[unitigs[current_unitig_idx]]!=current_unitig_idx){
            for(auto ele:unitig_nodes[current_unitig_idx].unitigs_after){
              if(ele>0){
                unitig_nodes[ele].unitigs_before.erase(current_unitig_idx);
              }else{
                unitig_nodes[-ele].unitigs_after.erase(-current_unitig_idx);
              }
            }
            unitig_nodes[current_unitig_idx].unitigs_after.clear();
            unitigs[current_unitig_idx] = "";
            continue;
          }
          get_unitig_backward1(cqf, params, unitigs, unitig_nodes, startKmer2unitig, endKmer2unitig, current_unitig_idx, unitigs[current_unitig_idx]);
        }else if(unitigs_before_num==0 && unitigs_after_num ==0){
          if(startKmer2unitig.find(unitigs[current_unitig_idx]) != startKmer2unitig.end()){
            //check for existence of the same unitigs
            if(startKmer2unitig[unitigs[current_unitig_idx]]!=current_unitig_idx){
              for(auto ele:unitig_nodes[current_unitig_idx].unitigs_before){
                if(ele>0){
                  unitig_nodes[ele].unitigs_after.erase(current_unitig_idx);
                }else{
                  unitig_nodes[-ele].unitigs_before.erase(-current_unitig_idx);  
                }
              }
              unitig_nodes[current_unitig_idx].unitigs_before.clear();
              unitigs[current_unitig_idx] = "";
              continue;
            }
            get_unitig_forward1(cqf, params, unitigs, unitig_nodes, startKmer2unitig, endKmer2unitig, current_unitig_idx, unitigs[current_unitig_idx]);
          }
          if(endKmer2unitig.find(unitigs[current_unitig_idx]) != endKmer2unitig.end()){
            //check for existence of the same unitigs
            if(endKmer2unitig[unitigs[current_unitig_idx]]!=current_unitig_idx){
              for(auto ele:unitig_nodes[current_unitig_idx].unitigs_after){
                if(ele>0){
                  unitig_nodes[ele].unitigs_before.erase(current_unitig_idx);
                }else{
                  unitig_nodes[-ele].unitigs_after.erase(-current_unitig_idx);
                }
              }
              unitig_nodes[current_unitig_idx].unitigs_after.clear();
              unitigs[current_unitig_idx] = "";
              continue;
            }
            get_unitig_backward1(cqf, params, unitigs, unitig_nodes, startKmer2unitig, endKmer2unitig, current_unitig_idx, unitigs[current_unitig_idx]);
          }
        }else{
          throw std::runtime_error("Unexpected cases: unextended unitig with non-zero before and after unitigs!");
        }
      }

      size_t unitig_num = unitigs.size();
      vector<bool> traveled(unitig_num, false); traveled[0]=true;
      
      /*
      cout<<"*"<<unitigs.size()<<" unitigs in total."<<endl;
      for(int x = 0; x< unitigs.size(); x++){
        cout<<"unitig "<<x<<endl;
        cout<<"\tseq: "<<unitigs[x]<<endl;
        cout<<"\tunitigs_before: ";
        for(auto ele:unitig_nodes[x].unitigs_before){
          cout<<ele<<" ";
        }
        cout<<endl;
        cout<<"\tunitigs_after: ";
        for(auto ele:unitig_nodes[x].unitigs_after){
          cout<<ele<<" ";
        }
        cout<<endl;
      }
      */
      remove_tips(params, unitigs, unitig_nodes, traveled);
      deduce_contigs(params, unitigs, unitig_nodes, contigs, traveled);
    }
  }
  /*
  cout<<"*"<<contigs.size()<<" contigs in total."<<endl;
  for(int x  = 0; x<contigs.size(); x++){
    cout<<"contig "<<x<<endl;
    cout<<"seq: "<<contigs[x]<<endl;
  }
  */
}


void ntHash(){
  string read = "ATGC";
  int K = 4;
  uint64_t hash, hash_RC;
  NTPC64(read.c_str(), K, hash, hash_RC);
  cout<<read<<endl;
  cout<<hash<<"\t"<<hash_RC<<endl;

  read = "atgc";
  NTPC64(read.c_str(), K, hash, hash_RC);
  cout<<read<<endl;
  cout<<hash<<"\t"<<hash_RC<<endl;
  
  read = "atgN";
  NTPC64(read.c_str(), K, hash, hash_RC);
  cout<<read<<endl;
  cout<<hash<<"\t"<<hash_RC<<endl;

  read = "atgn";
  NTPC64(read.c_str(), K, hash, hash_RC);
  cout<<read<<endl;
  cout<<hash<<"\t"<<hash_RC<<endl;
 
  read = "qftg";
  NTPC64(read.c_str(), K, hash, hash_RC);
  cout<<read<<endl;
  cout<<hash<<"\t"<<hash_RC<<endl;

  if(argc<3){
    cout<<"No enough parameters."<<endl;
    cout<<argv[0]<<" <K> <cqf2load>"<<endl;
  }

}
#endif
