/*
 * ============================================================================
 *  Filename:  unitig_graph.h
 *  
 *  Codes by
 *    "SH-assembly"
 *        Authors: Christina SHI <hshi@cse.cuhk.edu.hk>
 *                 Kevin Yip <kevinyip@cse.cuhk.edu.hk>
 *        
 * ============================================================================
 */

#pragma once

#include "base/global.h"
#include "base/Utility.h"
#include "base/DNA_string.h"
#include "base/Params.h"
//#include "base/Hash.h"

using tbb::concurrent_vector;
using tbb::concurrent_vector;
using tbb::concurrent_queue;

struct Edge{
  int toNode;
  // tbb::atomic<int> coverage;
  // Edge(int node, int c=0):toNode(node), coverage(c){}
  Edge(int node):toNode(node){}
};

struct UnitigNode{
  //vector<Edge> beforeNodes;
  //vector<Edge> afterNodes;
  vector<int> beforeNodes;
  vector<int> afterNodes;
  UnitigNode(){}
  // string to_str() const{
  //   string tmp="Nodes before: ";
  //   for(auto ele : beforeNodes){
  //     tmp += to_string(ele.toNode)+"("+to_string(ele.coverage)+") ";
  //   }
  //   tmp += "\nNodes after: ";
  //   for(auto ele : afterNodes){
  //     tmp += to_string(ele.toNode)+"("+to_string(ele.coverage)+") ";
  //   }
  //   return tmp+"\n";
  // }
  string to_str() const;
};

namespace std{
  template<>
  struct hash<DNAString>
  {
     size_t operator () (const DNAString& x) const
     {
        //cout<<"hash being called"<<endl;
        return MurmurHash2(x.data(), x.data_size());
     }
  };

  template<>
  struct equal_to<DNAString>
  {
    bool operator()(const DNAString &lhs, const DNAString &rhs) const 
    {
        //cout<<"equal being called"<<endl;
        return lhs == rhs;
    }
  };
}

namespace tbb{
  template<>
  struct tbb_hash<DNAString>
  { 
    size_t operator () (const DNAString& x) const
     {
        return MurmurHash2(x.data(), x.data_size());
     }
  };
}

struct DNAStringCompare {
    static size_t hash( const DNAString& x ) {
      return MurmurHash2(x.data(), x.data_size());
    }
    //! True if strings are equal
    static bool equal( const DNAString& x, const DNAString& y ) {
        return x==y;
    }
};

typedef tbb::concurrent_unordered_set<DNAString> unordered_set_mt;
//typedef tbb::concurrent_unordered_map<DNAString, int> unordered_map_mt;
typedef tbb::concurrent_hash_map<DNAString, int, DNAStringCompare> hash_map_mt;

//unitigs are stored from 1
bool load_unitig_graph(const string& file, concurrent_vector<Contig>& unitigs);
bool load_unitig_graph(const Params& options, const string& file, concurrent_vector<Contig>& unitigs, hash_map_mt& startKmer2unitig, concurrent_vector<UnitigNode>& unitigNodes);
