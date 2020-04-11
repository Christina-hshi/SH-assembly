/*
 * ============================================================================
 *  
 *        Authors: Christina SHI <hshi@cse.cuhk.edu.hk>
 *                 Kevin Yip <kevinyip@cse.cuhk.edu.hk>
 *
 * ============================================================================
 */

#ifndef _UNORDERED_MAP_MT_H
#define _UNORDERED_MAP_MT_H

#include<unordered_map>
#include<boost/thread.hpp>
#include<limits>

using namespace::std;
template < class Key,                                  // unordered_map::key_type
         class T,                                      // unordered_map::mapped_type
         class Hash = hash<Key>,                       // unordered_map::hasher
         class Pred = equal_to<Key>,                   // unordered_map::key_equal
         class Alloc = allocator< pair<const Key,T> >  // unordered_map::allocator_type
         > 
class unordered_map_mt : public unordered_map<Key, T, Hash, Pred, Alloc>{
  typedef unordered_map<Key, T, Hash, Pred, Alloc> _Base;

private:
  boost::shared_mutex* bucket_mutex;
public:
  typedef typename _Base::value_type value_type;
  typedef typename _Base::size_type  size_type;
  typedef typename _Base::hasher     hasher;
  typedef typename _Base::key_equal  key_equal;
  typedef typename _Base::allocator_type allocator_type;
  
  unordered_map_mt(){ bucket_mutex = NULL;}

  explicit unordered_map_mt( size_type n, 
      const hasher& hf = hasher(),
      const key_equal& eql = key_equal(),
      const allocator_type& alloc = allocator_type())
  :_Base(n, hf, eql, alloc){
    bucket_mutex = new boost::shared_mutex[n];
    //set the max_load_factor to as large as possible to avoid rehashing
    this->max_load_factor(numeric_limits<float>::max());//prevent from rehashing
  }
   
  void insert_mt( const Key& key, const T& value){
    size_t bucket_idx = this->bucket(key);
    boost::unique_lock<boost::shared_mutex> lock(bucket_mutex[bucket_idx]);
    this->insert(std::make_pair(key, value));
  }
 
  void erase_mt( const Key& key){
    size_t bucket_idx = this->bucket(key);
    boost::unique_lock<boost::shared_mutex> lock(bucket_mutex[bucket_idx]);
    this->erase(key);
  }

  bool find_mt( const Key& key, T& value){
    size_t bucket_idx = this->bucket(key);
    boost::shared_lock<boost::shared_mutex> lock(bucket_mutex[bucket_idx]);
    auto it = this->find(key);
    if(this->find(key) != this->end()){
      value = it->second;
      return true;
    }
    return false;
  }
 
  void insert_or_assign_mt(const Key& key, T& value){
    size_t bucket_idx = this->bucket(key);
    boost::unique_lock<boost::shared_mutex> lock(bucket_mutex[bucket_idx]);
    this->insert_or_assign(key, value);
  }

  ~unordered_map_mt(){
    if(bucket_mutex!=NULL){
      delete []bucket_mutex;
    }
  }
};

#endif
