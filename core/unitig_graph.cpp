/*
 * ============================================================================
 *  Filename:  unitig_graph.cpp
 *  
 *  Codes by
 *    "SH-assembly"
 *        Authors: Christina SHI <hshi@cse.cuhk.edu.hk>
 *                 Kevin Yip <kevinyip@cse.cuhk.edu.hk>
 *        
 * ============================================================================
 */

#include "unitig_graph.h"

string UnitigNode::to_str() const{
    string tmp="Nodes before: ";
    for(auto ele : beforeNodes){
      tmp += to_string(ele)+"|";
    }
    tmp.resize(tmp.length()-1);
    tmp += "\nNodes after: ";
    for(auto ele : afterNodes){
      tmp += to_string(ele)+"|";
    }
    tmp.resize(tmp.length()-1);
    return tmp+"\n";
  }

bool load_unitig_graph(const string& file, concurrent_vector<Contig>& unitigs){
	unitigs.resize(1);
	ifstream fin;
	fin.open(file, std::ios::in);
	if(fin.is_open()){
		string line, seq;
		int ave_abundance=0;
		while(getline(fin, line)){
			if(line.empty()){
	      continue;
	    }else if(line[0] != '>'){
	      continue;
	    }
	    int pos = line.find("km:f:");
	    if(pos == std::string::npos){
	    	std::cerr<<"[Load unitig graph] the mean abundance of k-mers in unitigs should be provided."<<endl;
	    	return false;
	    }else{
	    	pos += 5;
	    	int end = line.find(" ", pos);
	    	ave_abundance = std::stoi(line.substr(pos, end-pos));
	    }
	    getline(fin, seq);
	    unitigs.push_back(Contig(seq, ave_abundance));
		}
		fin.close();
		return true;
	}
	return false;
}

bool load_unitig_graph(const Params& options, const string& file, concurrent_vector<Contig>& unitigs, hash_map_mt& startKmer2unitig, concurrent_vector<UnitigNode>& unitigNodes){
	unitigs.resize(1);
	unitigNodes.resize(1);
	ifstream fin;
	fin.open(file, std::ios::in);
	int contig_idx = 0;
	if(fin.is_open()){
		string line, seq;
		int ave_abundance=0;
		while(getline(fin, line)){
			if(line.empty()){
	      continue;
	    }else if(line[0] != '>'){
	      continue;
	    }
	    int pos = line.find("km:f:");
	    if(pos == std::string::npos){
	    	std::cerr<<"[Load unitig graph] the mean abundance of k-mers in unitigs should be provided."<<endl;
	    	return false;
	    }else{
	    	pos += 5;
	    	int end = line.find(" ", pos);
	    	ave_abundance = std::stoi(line.substr(pos, end-pos));
	    }
	    getline(fin, seq);
	    contig_idx++;
	    unitigs.push_back(Contig(seq, ave_abundance));
	    //TODO: load connection between unitigs

	    //trace the start kmer of contig
	    DNAString kmer, last_kmer_RC;
	    kmer = seq.substr(0, options.K);
	    hash_map_mt::accessor access;
      if(startKmer2unitig.insert(access, kmer)){
      	access->second = contig_idx; 
      	access.release();
      }
      last_kmer_RC = seq.substr(seq.length() - options.K);
      last_kmer_RC.RC();
      if(!(kmer == last_kmer_RC)){
	      if(startKmer2unitig.insert(access, last_kmer_RC)){
	      	access->second = -contig_idx; 
	      	access.release();
	      }
	    }
		}
		fin.close();
		return true;
	}
	return false;
}