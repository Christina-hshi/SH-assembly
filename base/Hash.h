#pragma once
/**
* For hashing the DNA sequence to generate fingerprint for each sequence, which are used to approximate sequence similarity
*/

#include "global.h"
//#include "Utility.h"

uint32_t MurmurHash2 ( const void * key, int len, u_int seed=0);

uint32_t MurmurHash2 ( string seq, u_int seed=0);

/*
uint32_t MurmurHash2_kmer( const string seq, int k, u_int seed);

Apply MurmurHash2 to all kmers in seq and reverse complement of seq.
And store "sketch_size" smallest hash values in sketch.
#duplicate hashValues are discarded.

void MinHash_MurmurHash2_kmer(const string seq, int k, u_int seed, u_int* sketch, int sketch_size);
*/
