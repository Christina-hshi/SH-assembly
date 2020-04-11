/*
 * ============================================================================
 *  
 *        Authors: Christina SHI <hshi@cse.cuhk.edu.hk>
 *                 Kevin Yip <kevinyip@cse.cuhk.edu.hk>
 *
 * ============================================================================
 */

#pragma once

#include "global.h"
#include "Utility.h"
#include "DNA_string.h"
#include "Params.h"
#include "cqf/CQF_mt.h"

using tbb::concurrent_vector;
using tbb::concurrent_vector;
using tbb::concurrent_queue;

struct Edge{
  int toNode;
  int coverage;
  Edge(int node, int c=0):toNode(node), coverage(c){

  }
};

struct UnitigNode{
  vector<Edge> beforeNodes;
  vector<Edge> afterNodes;
  UnitigNode(){}
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
