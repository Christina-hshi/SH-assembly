/*
 * ============================================================================
 *  Filename:  chunk.h
 *  
 *  Mixture of codes by
 *    "CQF-deNoise"
 *        Authors: Christina SHI <hshi@cse.cuhk.edu.hk>
 *                 Kevin Yip <kevinyip@cse.cuhk.edu.hk>
 *        
 *    "A General-Purpose Counting Filter: Making Every Bit Count"    
 *        Authors:  Prashant Pandey <ppandey@cs.stonybrook.edu>
 *                  Rob Johnson <robj@vmware.com>   
 *
 * ============================================================================
 */

#ifndef _CHUNK_H_
#define _CHUNK_H_

#include <stdio.h>

namespace kmercounting {
	class chunk {
		public:
			chunk();
			chunk(char *reads, uint32_t size);
			inline char *get_reads();
			inline uint32_t get_size();
			inline void set(char *reads, uint32_t size);
			bool readLine(string& str);
			bool skipLines(uint32_t num=1);//skip lines;
			inline void reset_pointer();

		private:
			char *_reads;
			uint32_t _size;
			char *_end;//end of the chunk
			char *_pointer;//reading from pointer
	};

	chunk::chunk()
	{
		_reads = NULL;
		_size = 0;
		_pointer = NULL; _end=NULL;
	}

	chunk::chunk(char *reads, uint32_t size)
	{
		_reads = reads;
		_size = size;
		_end = reads + size;
		_pointer = reads; 
	}

	inline char *chunk::get_reads()
	{
		return _reads;
	}

	inline uint32_t chunk::get_size()
	{
		return _size;
	}

	inline void chunk::set(char *reads, uint32_t size){
		_reads = reads;
		_size = size;
		_end = reads + size;
		_pointer = reads;
	}
	bool chunk::readLine(string& str){
		if(_pointer == _end){
			return false;
		}
		auto tmp = static_cast<char*>(memchr(_pointer, '\n', this->_end - this->_pointer));
		str.assign(_pointer, tmp-this->_pointer);
		_pointer = tmp + 1;
		return true;
	}
	bool chunk::skipLines(uint32_t num){
		while(num--){
			if(_pointer == _end){
				return false;
			}
			auto tmp = static_cast<char*>(memchr(_pointer, '\n', this->_end - this->_pointer));
			_pointer = tmp+1;
		}
		return true;
	}
	inline void chunk::reset_pointer(){
		_pointer = _reads;
	}
}

#endif
