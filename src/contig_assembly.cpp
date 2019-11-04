#define GRAPH_TRAVERSE

#include "cqf/CQF_mt.h"
#include "base/Utility.h"
#include "base/unordered_map_mt.h"
#include "base/vector_mt.h"
#include "base/Params.h"
//#include "base/nthash.h"
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>

//using namespace boost::program_options;
//namespace po = boost::program_options;
//using boost::program_options::variables_map;

boost::program_options::variables_map get_opts(int argc, char* argv[]){
  namespace po = boost::program_options;
  po::options_description desc(string(argv[0])+"  <options>\nOptions:");
  desc.add_options() 
    ("help,h", "print help messages") 
    (",k", po::value<int>()->required(), "kmer length") 
    ("input,i", po::value<string>()->required(), "a file containing list of input file name(s), should be absolute address or file names when in the running directory.")
    ("format,f", po::value<char>()->default_value('f'), "format of the input: g(gzip); b(bzip2); f(plain fastq)")
    ("cqf,f", po::value<string>()->required(), "the counting quotient filter built with the same 'k'")
    ("abundance_min, s", po::value<int>()->default_value(2), "minimum coverage of a solid k-mer to start the assembly of a contig") 
    (",t", po::value<int>()->default_value(16), "number of threads")
    ("output,o", po::value<string>(), "output contig file name (fasta)");
    
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  
  if (argc==1 || vm.count("help")) {
    cerr << desc << "\n";
    exit(0);
  }
  
  po::notify(vm);

  return vm;
}

struct WorkQueue;
//void find_unitigs_mt_master(CQF_mt& cqf, const vector<string>& seqFiles, const Params& params, vector_mt<string>& contigs);
void find_unitigs_mt_master(CQF_mt& cqf, seqFile_batch& seqFiles, const boost::program_options::variables_map& options, vector_mt<string>& unitigs);

//void find_unitigs_mt_worker(CQF_mt& cqf, const Params& params, vector_mt<string>& contigs, unordered_map_mt<string, int>& startKmer2unitig, WorkQueue* work_queue);
void find_unitigs_mt_worker(CQF_mt& cqf, const boost::program_options::variables_map& options, vector_mt<string>& contigs, unordered_map_mt<string, int>& startKmer2unitig, WorkQueue* work_queue);

//void get_unitig_forward(CQF_mt& cqf, const Params& params, vector_mt<string>& contigs, unordered_map_mt<string, int>& startKmer2unitig, WorkQueue* work_queue, int contig_id);
void get_unitig_forward(CQF_mt& cqf, const boost::program_options::variables_map& options, vector_mt<string>& contigs, unordered_map_mt<string, int>& startKmer2unitig, WorkQueue* work_queue, int contig_id);

//void get_unitig_backward(CQF_mt& cqf, const Params& params, vector_mt<string>& contigs, unordered_map_mt<string, int>& startKmer2unitig, WorkQueue* work_queue, int contig_id);
void get_unitig_backward(CQF_mt& cqf, const boost::program_options::variables_map& options, vector_mt<string>& contigs, unordered_map_mt<string, int>& startKmer2unitig, WorkQueue* work_queue, int contig_id);

int main(int argc, char* argv[]){
  //string fileList = "/research/kevinyip11/hshi/NCBI/F.vesca/files.list1";
  //string fileList = "/research/kevinyip11/hshi/Tools/art_src_MountRainier_Linux/44_data/chr6_cox_hap2/MSv1_C50_L250/files.list";
  //string cqfFile = "/research/kevinyip11/hshi/Tools/squeakr/F.vesca.k55.s33.nthash.ser";
  //string prefix = "/research/kevinyip11/hshi/Tools/art_src_MountRainier_Linux/44_data/chr6_cox_hap2/MSv1_C50_L250/single_end";
  //string dir="/research/kevinyip11/hshi/Tools/art_src_MountRainier_Linux/44_data/chr6_cox_hap2/MSv1_C50_L250/";
  //tring dir="/research/kevinyip11/hshi/NCBI/F.vesca/";
  boost::program_options::variables_map options = get_opts(argc, argv);

  vector<string> seqFileNames;
  ifstream fin;
  fin.open(options["input"].as<string>(), ios::in);
  string line;
  while(getline(fin, line)){
    if(line.empty())
      continue;
    seqFileNames.push_back(line);
  }
  FILE_TYPE ftype=FILE_TYPE::FASTQ;
  FILE_MODE fmode;
  switch(options["format"].as<char>()){
    case 'g':
      fmode = FILE_MODE::GZIP;
      break;
    case 'b':
      fmode = FILE_MODE::BZIP2;
      break;
    case 'f':
      fmode = FILE_MODE::TEXT;
      break;
    default:
      break;
  }
  seqFile_batch seqFiles(seqFileNames, ftype, fmode);   

  DisplayCurrentDateTime(); 
  cout<<"[CQF] loading cqf from disk"<<endl;
  CQF_mt cqf_mt;
  cqf_mt.load(options["cqf"].as<string>());
  cout<<"[CQF] cqf loaded!"<<endl;

  DisplayCurrentDateTime();
  cout<<"[Unitig] constructing unitigs"<<endl;
  vector_mt<string> contigs;
  contigs.resize(1);
  find_unitigs_mt_master(cqf_mt, seqFiles, options, contigs);
  
  DisplayCurrentDateTime();  
  contig_summary(contigs);
  return 0;
}

//use 1-based counter, e.g. the first work has work_id 1.
struct WorkQueue{
  volatile atomic<uint32_t> next_work;
  volatile atomic<uint32_t> total_work;
  volatile atomic<uint32_t> work_done;
  volatile atomic<bool> master_done;
  boost::mutex mut;

  WorkQueue(){
    next_work = 0;
    total_work = 0;
    master_done = false;
  }
  WorkQueue(uint32_t n, uint32_t t){
    next_work = n;
    total_work = t;
    master_done = false;
  }
  
  bool get_next_work(uint32_t& work_id){
    boost::unique_lock<boost::mutex> lock(mut);
    if(next_work < total_work){
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
    lock.unlock();
  }
};

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

void find_unitigs_mt_master(CQF_mt& cqf, seqFile_batch& seqFiles, const boost::program_options::variables_map& options, vector_mt<string>& contigs){
  int K = options["-k"].as<int>();
  int thread_num = options["-t"].as<int>();

  WorkQueue* work_queue = new WorkQueue();
  unordered_map_mt<string, int> startKmer2unitig(1000);

  int abundance_min = options["abundance_min"].as<int>();

  boost::thread_group prod_threads;
  for(int t = 0; t<thread_num; t++){
    prod_threads.add_thread(new boost::thread(find_unitigs_mt_worker, boost::ref(cqf), boost::ref(options), boost::ref(contigs), boost::ref(startKmer2unitig), work_queue));
  }
  
  string kmer, kmer_RC;
  uint64_t kmer_hash, kmer_RC_hash;
  uint64_t kmer_count, kmer_RC_count;
  size_t contig_id;
  
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
      if(seq.length()<K){
        continue;
      }
      //seq_RC = RC_DNA(seq);
      int seq_len = seq.length();
      int step = seq_len/3;
      for(int x = 0; x<=seq_len-K; x += step){
        kmer = seq.substr(x, K);
        kmer_RC = RC_DNA(kmer);//seq_RC.substr(seq_len-params.K-x, params.K);
        if(kmer.find_first_of("nN")!=string::npos){
          continue;
        }
        //kmer_hash = MurmurHash2(kmer);
        //kmer_RC_hash = MurmurHash2(kmer_RC);
        NTPC64(kmer.c_str(), K, kmer_hash, kmer_RC_hash);
        
        if(cqf.count_key_value_set_traveled(kmer_hash%cqf.qf->metadata->range, kmer_count) && cqf.count_key_value_set_traveled(kmer_RC_hash%cqf.qf->metadata->range, kmer_RC_count)){
          continue;
        }else if(kmer_count < abundance_min || kmer_RC_count < abundance_min){
          continue;
        }

        contig_id = contigs.push_back_mt(seq);
        work_queue->add_skip_work(1);   
        startKmer2unitig.insert_mt(kmer, contig_id); //may not be the start k-mer
        startKmer2unitig.insert_mt(kmer_RC, -contig_id); //may not be the end k-mer

        get_unitig_forward(cqf, options, contigs, startKmer2unitig, work_queue, contig_id);
        get_unitig_backward(cqf, options, contigs, startKmer2unitig, work_queue, contig_id);

        //wait for the end of this run.
        while(work_queue->work_done != work_queue->total_work){
          std::this_thread::sleep_for(std::chrono::milliseconds(1000));//sleep for 1 seconds  
        }
      }
      //skip two lines
      dataChunk.skipLines(2);
    }
  }
  work_queue->master_done = true;
  prod_threads.join_all();
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

void find_unitigs_mt_worker(CQF_mt& cqf, const boost::program_options::variables_map& options, vector_mt<string>& contigs, unordered_map_mt<string, int>& startKmer2unitig, WorkQueue* work_queue){
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

void get_unitig_forward(CQF_mt& cqf, const boost::program_options::variables_map& options, vector_mt<string>& contigs, unordered_map_mt<string, int>& startKmer2unitig, WorkQueue* work_queue, int contig_id){
  array<bool, 4> candidates_before({false, false, false, false}), candidates_after({false, false, false, false});
  array<string, 4> kmer_befores, kmer_afters, kmer_befores_RC, kmer_afters_RC;
  auto abundance_min = options["abundance_min"].as<int>();
  auto K = options["k"].as<int>();

  int candidates_before_num, candidates_after_num;
  int nodes_before_num, nodes_after_num;
  uint64_t kmer_hash, kmer_RC_hash, current_kmer_hash, current_kmer_RC_hash;
  string kmer, kmer_RC, current_kmer;
  uint64_t kmer_count, kmer_RC_count;
  int idx;
  
  string contig_seq = contigs.at_mt(contig_id);
  current_kmer = contig_seq.substr(contig_seq.length()-K);

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

      NTPC64(kmer.c_str(), K, kmer_hash, kmer_RC_hash);
      //kmer_hash = MurmurHash2(kmer);
      //kmer_RC_hash = MurmurHash2(kmer_RC);
      if(cqf.count_key_value_set_traveled(kmer_hash%cqf.qf->metadata->range, kmer_count) && cqf.count_key_value_set_traveled(kmer_RC_hash%cqf.qf->metadata->range, kmer_RC_count)){
        if(startKmer2unitig.find_mt(kmer, idx)){
          nodes_after_num++;
        }else if(kmer_count >= abundance_min && kmer_RC_count >= abundance_min){
          candidates_after[x] = true; //possible because of hash collisions
          candidates_after_num ++;
        }
      }else if(kmer_count >= abundance_min && kmer_RC_count >= abundance_min){
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

      NTPC64(kmer.c_str(), K, kmer_hash, kmer_RC_hash);
      //kmer_hash = MurmurHash2(kmer);
      //kmer_RC_hash = MurmurHash2(kmer_RC);
      if(cqf.count_key_value_set_traveled(kmer_hash%cqf.qf->metadata->range, kmer_count) && cqf.count_key_value_set_traveled(kmer_RC_hash%cqf.qf->metadata->range, kmer_RC_count)){
        if(startKmer2unitig.find_mt(kmer_RC, idx)){
          nodes_before_num++;
        }else if(kmer_count >= abundance_min && kmer_RC_count >= abundance_min){
          candidates_before[x] = true;
          candidates_before_num++;
        }
      }else if(kmer_count >= abundance_min && kmer_RC_count >= abundance_min){
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
      startKmer2unitig.insert_mt(RC_DNA(current_kmer), -contig_id);
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
            //startKmer2unitig[kmer_befores_RC[x]] = -new_unitig_idx;
            startKmer2unitig.insert_mt(kmer_befores_RC[x], -new_unitig_idx);
            work_queue->add_work(1);
          }
        }
      }
      break;
    }else if(candidates_after_num==1){ //only one candidate k-mer after
      for(int x = 0; x<4; x++){
        if(candidates_after[x]){
          current_kmer = kmer_afters[x];
          contig_seq += current_kmer.back();
          break;
        }
      }
      continue;
    }else{ //stop
      startKmer2unitig.insert_mt(RC_DNA(current_kmer), -contig_id);
      contigs.set_mt(contig_id, contig_seq);
      break;
    }
  }
}

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

void get_unitig_backward(CQF_mt& cqf, const boost::program_options::variables_map& options, vector_mt<string>& contigs, unordered_map_mt<string, int>& startKmer2unitig, WorkQueue* work_queue, int contig_id){
  auto abundance_min = options["abundance_min"].as<int>();
  auto K = options["k"].as<int>();

  array<bool, 4> candidates_before, candidates_after;
  array<string, 4> kmer_befores, kmer_afters, kmer_befores_RC, kmer_afters_RC;
  int candidates_before_num, candidates_after_num;
  int nodes_before_num, nodes_after_num;
  uint64_t kmer_hash, kmer_RC_hash;
  string kmer, kmer_RC, current_kmer;
  uint64_t kmer_count, kmer_RC_count;
  int idx;

  string contig_seq = contigs.at_mt(contig_id); 
  current_kmer = contig_seq.substr(0, K);
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

      NTPC64(kmer.c_str(), K, kmer_hash, kmer_RC_hash);
      //kmer_hash = MurmurHash2(kmer);
      //kmer_RC_hash = MurmurHash2(kmer_RC);
      if(cqf.count_key_value_set_traveled(kmer_hash%cqf.qf->metadata->range, kmer_count) && cqf.count_key_value_set_traveled(kmer_RC_hash%cqf.qf->metadata->range, kmer_RC_count)){
        if(startKmer2unitig.find_mt(kmer, idx)){
          nodes_after_num++;
        }else if(kmer_count >= abundance_min && kmer_RC_count >= abundance_min){
          candidates_after[x] = true;
          candidates_after_num ++;
        }
      }else if(kmer_count >= abundance_min && kmer_RC_count >= abundance_min){
        candidates_after[x] = true;
        candidates_after_num ++;        
      }
    }  
    for(int x = 0; x<4; x++){
      kmer = kmer_befores[x];
      kmer_RC = RC_DNA(kmer);
      kmer_befores_RC[x] = kmer_RC;

      NTPC64(kmer.c_str(), K, kmer_hash, kmer_RC_hash);
      //kmer_hash = MurmurHash2(kmer);
      //kmer_RC_hash  = MurmurHash2(kmer_RC);
      if(cqf.count_key_value_set_traveled(kmer_hash%cqf.qf->metadata->range, kmer_count) && cqf.count_key_value_set_traveled(kmer_RC_hash%cqf.qf->metadata->range, kmer_RC_count)){
        if(startKmer2unitig.find_mt(kmer_RC, idx)){
          nodes_before_num++; 
        }else if (kmer_count >= abundance_min && kmer_RC_count >= abundance_min){
          candidates_before[x] = true;
          candidates_before_num++;
        }
      }else if(kmer_count >= abundance_min && kmer_RC_count >= abundance_min){
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
