/*
 * ============================================================================
 *  
 *        Authors: Christina SHI <hshi@cse.cuhk.edu.hk>
 *                 Kevin Yip <kevinyip@cse.cuhk.edu.hk>
 *
 * ============================================================================
 */

#ifndef _VECTOR_MT_H
#define _VECTOR_MT_H
#include<vector>
#include<mutex>
#include<boost/thread.hpp>
#include<atomic>

using namespace::std;
template < class T, class Alloc = allocator<T> >
class vector_mt : public vector<T, Alloc>{
private:
  boost::shared_mutex mutex_mt;
public:
  vector_mt(const Alloc& __a = Alloc()):vector<T, Alloc>(__a){
  }
  
  vector_mt(const vector_mt& v) = delete;
  
  //return the current size
  size_t push_back_mt(const T& x){
    size_t s;
    boost::unique_lock<boost::shared_mutex> lock(mutex_mt);
    s = this->size();
    this->push_back(x);
    lock.unlock();
    return s;
  }
  
  T at_mt(size_t n){
    boost::shared_lock<boost::shared_mutex> lock(mutex_mt);
    return this->at(n); 
  }
 
  void set_mt(size_t n, T value){
    boost::shared_lock<boost::shared_mutex> lock(mutex_mt);
    this->at(n) = value; 
  }

  size_t size_mt(){
    boost::shared_lock<boost::shared_mutex> lock(mutex_mt);
    return this->size();
  }
};

#endif
