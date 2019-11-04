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

#define mt_log multithread_io::log

class multithread_io{
  public:
    static std::mutex mtx;
    static void log(string m){
      std::unique_lock<std::mutex> locker(mtx);
      cout<<m<<std::flush;
      locker.unlock();
    }
};
