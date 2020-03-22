/** For storing global settings
*/

#pragma once

#include<string>
#include<iostream>
#include<sstream>
#include<fstream>
#include<vector>
#include<map>
#include<unordered_map>
#include<set>
#include<list>
#include<stack>
#include<fstream>
#include<sstream>
#include<stdexcept>
#include<ctime>
#include<cstdlib>
//#include<cmath>
#include<cstring>
#include<thread>
#include<chrono>
#include<mutex>
#include<condition_variable>
#include<algorithm>
#include<limits>
#include<numeric>
#include<assert.h>
//#include<malloc.h>
//#include<malloc/malloc.h>
//#include<stdlib.h>
#include<queue>
//#include<atomic>
#include<sys/resource.h>
#include<sys/time.h>

#include<boost/config.hpp>
#include<boost/tokenizer.hpp>
#include<boost/bimap.hpp>
#include<boost/program_options.hpp>
//#include "Config.h"
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/atomic.hpp>

#include <tbb/concurrent_unordered_set.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_queue.h>
#include <tbb/concurrent_hash_map.h>
//#include "tbb/atomic.h"

//modules used as default in this project
using std::cin;
using std::cout;
using std::endl;
using std::string;
using std::ostream;
using std::ifstream;
using std::ofstream;
using std::fstream;
using std::istringstream;
using std::ostringstream;
using std::stringstream;
using std::vector;
using std::map;
using std::unordered_map;
using std::set;
using std::list;
using std::pair;
using std::bitset;
using std::array;
using std::stack;
using std::to_string;

using std::streampos;

typedef unsigned int u_int;
typedef unsigned long int ul_int;
typedef long long int ll_int;

typedef unsigned int SEQ_ID_T; //sequence id type: to adapt to different number of sequences
typedef int CHR_LEN_T; //chromosome len type(used as offset and s_pos): negative value is needed.
typedef u_int* sketch_t; //sketch type

typedef vector<string> vec_str;
//#define EN_DEBUG_MSG

typedef boost::tokenizer<boost::char_separator<char>> tokenizer;

enum SEQ_LIB_TYPE{MP, PE, SE};
enum CQF_MODE{EMPTY, FILEMAP, MEMORY};
enum FILE_TYPE{FASTA, FASTQ};
enum FILE_MODE{TEXT, GZIP, BZIP2};  


const char DNA_bases[4] = {'A', 'C', 'G','T'};

const std::string currentDateTime();


#if 0
std::ostream& operator<<( std::ostream& dest, __int128_t value ){
  std::ostream::sentry s( dest );
  if ( s ) {
    __uint128_t tmp = value < 0 ? -value : value;
    char buffer[ 128 ];
    char* d = std::end( buffer );
    do
    {
      -- d;
      *d = "0123456789"[ tmp % 10 ];
      tmp /= 10;
    } while ( tmp != 0 );
    if ( value < 0 ) {
      -- d;
      *d = '-';
    }
    int len = std::end( buffer ) - d;
    if ( dest.rdbuf()->sputn( d, len ) != len ) {
      dest.setstate( std::ios_base::badbit );
    }
  }
  return dest;
}

std::ostream& operator<<( std::ostream& dest, __uint128_t value ){
  std::ostream::sentry s( dest );
  if ( s ) {
    __uint128_t tmp = value;
    char buffer[ 128 ];
    char* d = std::end( buffer );
    do
    {
      -- d;
      *d = "0123456789"[ tmp % 10 ];
      tmp /= 10;
    } while ( tmp != 0 );
    if ( value < 0 ) {
      -- d;
      *d = '-';
    }
    int len = std::end( buffer ) - d;
    if ( dest.rdbuf()->sputn( d, len ) != len ) {
      dest.setstate( std::ios_base::badbit );
    }
  }
  return dest;
}
#endif
