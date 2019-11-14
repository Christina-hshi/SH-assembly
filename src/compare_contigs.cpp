#include "base/global.h"

boost::program_options::variables_map get_opts(int argc, char* argv[]){
  namespace po = boost::program_options;
  po::options_description desc(string(argv[0])+"  <options>\nOptions:");
  desc.add_options() 
    ("help,h", "print help messages") 
    (",k", po::value<int>()->required(), "kmer length, kmers are used as seed to index sequences.") 
    ("input1,1", po::value<string>()->required(), "first DNA sequence file in fasta format")
    ("input2,2", po::value<string>()->required(), "second DNA sequence file in fasta format")
    ("output,o", po::value<string>()->default_value("unitigs.fa"), "output shared DNA sequences");
    
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  
  if (argc==1 || vm.count("help")) {
    cerr << desc << "\n";
    exit(0);
  }
  
  po::notify(vm);
  return vm;
}

int main(int argc, char* argv[]){
  namespace po = boost::program_options;
  po::variables_map options = get_opts(argc, argv);

  int k = options["-k"].as<int>();
  string input1 = options["input1"].as<string>();
  string input2 = options["input2"].as<string>();
  string output = options["output"].as<string>();

  ifstream fin1, fin2;
  ofstream fout;
  fin1.open(input1, ios::in);
  fin2.open(input2, ios::in);
  fout.open(output, ios::out);

  //load all sequences 
  std::vector<string> seqs1, seqs2;
  string line, seq;
  getline(fin1, line);
  while(getline(fin1, line)){
    seq = line;
    while(getline(fin1, line)){
      if(line[0]=='>') break;
      seq += line;
    }
    seqs1.push_back(seq);
  }
  getline(fin2, line);
  while(getline(fin2, line)){
    seq = line;
    while(getline(fin2, line)){
      if(line[0]=='>') break;
      seq += line;
    }
    seqs2.push_back(seq);
  }

  fin1.close();
  fin2.close();
  fout.close();

  return 0;
}