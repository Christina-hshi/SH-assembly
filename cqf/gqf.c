/*
 * ============================================================================
 *  Filename:  gqf.h
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


#include <stdlib.h>
#if 0
# include <assert.h>
#else
# define assert(x)
#endif
#include <string.h>
#include <inttypes.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "gqf.h"

/******************************************************************
 * Code for managing the metadata bits and slots w/o interpreting *
 * the content of the slots.
 ******************************************************************/

//#define LOG_WAIT_TIME

/* Must be >= 6.  6 seems fastest. */
//#define BLOCK_OFFSET_BITS (6)
#define BLOCK_OFFSET_BITS (8)
#define BLOCK_OFFSET_MAX ((1ULL<<(BLOCK_OFFSET_BITS))-1ULL)

//#define SLOTS_PER_BLOCK (1ULL << BLOCK_OFFSET_BITS)
//#define SLOTS_PER_BLOCK (1ULL << 6)
#define METADATA_WORDS_PER_BLOCK ((SLOTS_PER_BLOCK + 63) / 64)

//#define NUM_SLOTS_TO_LOCK (1ULL<<16)
#define CLUSTER_SIZE (1ULL<<14)

#define METADATA_WORD(qf,field,slot_index) (get_block((qf), (slot_index) / \
					SLOTS_PER_BLOCK)->field[((slot_index)  % SLOTS_PER_BLOCK) / 64])

#define BITMASK(nbits) ((nbits) == 64 ? 0xffffffffffffffff : (1ULL << (nbits)) \
												- 1ULL)
#define MAX_VALUE(nbits) ((1ULL << (nbits)) - 1)
#define BILLION 1000000000L

typedef struct __attribute__ ((__packed__)) qfblock {
	/* Code works with uint16_t, uint32_t, etc, but uint8_t seems just as fast as
   * anything else */
	uint8_t offset; 
	uint64_t occupieds[METADATA_WORDS_PER_BLOCK];
	uint64_t runends[METADATA_WORDS_PER_BLOCK];
#ifdef GRAPH_TRAVERSE
	uint64_t traveled[METADATA_WORDS_PER_BLOCK];
#endif

#if BITS_PER_SLOT == 8
	uint8_t  slots[SLOTS_PER_BLOCK];
#elif BITS_PER_SLOT == 16
	uint16_t  slots[SLOTS_PER_BLOCK];
#elif BITS_PER_SLOT == 32
	uint32_t  slots[SLOTS_PER_BLOCK];
#elif BITS_PER_SLOT == 64
	uint64_t  slots[SLOTS_PER_BLOCK];
#elif BITS_PER_SLOT != 0
	uint8_t   slots[SLOTS_PER_BLOCK * BITS_PER_SLOT / 8];
#else
	uint8_t   slots[];
#endif
} qfblock;

static __inline__ unsigned long long rdtsc(void)
{
	unsigned hi, lo;
	__asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
	return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}

#ifdef LOG_WAIT_TIME
/*static inline*/ bool qf_spin_lock(QF *cf, volatile int *lock, uint64_t idx, bool
																flag_spin)
{
	struct timespec start, end;
	bool ret;

	clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);
	if (!flag_spin) {
		ret = !__sync_lock_test_and_set(lock, 1);
		clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end);
		cf->mem->wait_times[idx].locks_acquired_single_attempt++;
		cf->mem->wait_times[idx].total_time_single += BILLION * (end.tv_sec -
																												start.tv_sec) +
			end.tv_nsec - start.tv_nsec;
	} else {
		if (!__sync_lock_test_and_set(lock, 1)) {
			clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end);
			cf->mem->wait_times[idx].locks_acquired_single_attempt++;
			cf->mem->wait_times[idx].total_time_single += BILLION * (end.tv_sec -
																													start.tv_sec) +
			end.tv_nsec - start.tv_nsec;
		} else {
			while (__sync_lock_test_and_set(lock, 1))
				while (*lock);
			clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end);
			cf->mem->wait_times[idx].total_time_spinning += BILLION * (end.tv_sec -
																														start.tv_sec) +
				end.tv_nsec - start.tv_nsec;
		}
		ret = true;
	}
	cf->mem->wait_times[idx].locks_taken++;

	return ret;

	/*start = rdtsc();*/
	/*if (!__sync_lock_test_and_set(lock, 1)) {*/
		/*clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end);*/
		/*cf->mem->wait_times[idx].locks_acquired_single_attempt++;*/
		/*cf->mem->wait_times[idx].total_time_single += BILLION * (end.tv_sec -
		 * start.tv_sec) + end.tv_nsec - start.tv_nsec;*/
	/*} else {*/
		/*while (__sync_lock_test_and_set(lock, 1))*/
			/*while (*lock);*/
		/*clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end);*/
		/*cf->mem->wait_times[idx].total_time_spinning += BILLION * (end.tv_sec -
		 * start.tv_sec) + end.tv_nsec - start.tv_nsec;*/
	/*}*/

	/*end = rdtsc();*/
	/*cf->mem->wait_times[idx].locks_taken++;*/
	/*return;*/
}
#else
/**
 * Try to acquire a lock once and return even if the lock is busy.
 * If spin flag is set, then spin until the lock is available.
 */
/*static inline*/ bool qf_spin_lock(volatile int *lock, bool flag_spin)
{
	if (!flag_spin) {
		return !__sync_lock_test_and_set(lock, 1);
	} else {
		while (__sync_lock_test_and_set(lock, 1))
			while (*lock);
		return true;
	}

	return false;
}
#endif

/*static inline*/ void qf_spin_unlock(volatile int *lock)
{
	__sync_lock_release(lock);
	return;
}

static bool qf_lock(QF *cf, uint64_t hash_bucket_index, bool spin, bool flag)
{
	uint64_t hash_bucket_lock_offset  = hash_bucket_index % NUM_SLOTS_TO_LOCK;
	if (flag) {
#ifdef LOG_WAIT_TIME
		if (!qf_spin_lock(cf, &cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK],
											hash_bucket_index/NUM_SLOTS_TO_LOCK, spin))
			return false;
		if (NUM_SLOTS_TO_LOCK - hash_bucket_lock_offset <= CLUSTER_SIZE) {
			if (!qf_spin_lock(cf, &cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1],
												hash_bucket_index/NUM_SLOTS_TO_LOCK+1, spin)) {
				qf_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
				return false;
			}
		}
#else
		if (!qf_spin_lock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK], spin))
			return false;
		if (NUM_SLOTS_TO_LOCK - hash_bucket_lock_offset <= CLUSTER_SIZE) {
			if (!qf_spin_lock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1],
												spin)) {
				qf_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
				return false;
			}
		}
#endif
	} else {
		/* take the lock for two lock-blocks; the lock-block in which the
		 * hash_bucket_index falls and the next lock-block */

#ifdef LOG_WAIT_TIME
		if (hash_bucket_index >= NUM_SLOTS_TO_LOCK && hash_bucket_lock_offset <=
				CLUSTER_SIZE) {
			if (!qf_spin_lock(cf, &cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK-1], hash_bucket_index/NUM_SLOTS_TO_LOCK-1, spin))
				return false;
		}
		if (!qf_spin_lock(cf, &cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK], hash_bucket_index/NUM_SLOTS_TO_LOCK, spin)) {
			if (hash_bucket_index >= NUM_SLOTS_TO_LOCK && hash_bucket_lock_offset <=
					CLUSTER_SIZE)
				qf_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK-1]);
			return false;
		}
		if (!qf_spin_lock(cf, &cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1], hash_bucket_index/NUM_SLOTS_TO_LOCK+1,
											spin)) {
			qf_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
			if (hash_bucket_index >= NUM_SLOTS_TO_LOCK && hash_bucket_lock_offset <=
					CLUSTER_SIZE)
				qf_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK-1]);
			return false;
		}
#else
		if (hash_bucket_index >= NUM_SLOTS_TO_LOCK && hash_bucket_lock_offset <=
				CLUSTER_SIZE) {
			if (!qf_spin_lock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK-1], spin))
				return false;
		}
		if (!qf_spin_lock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK], spin)) {
			if (hash_bucket_index >= NUM_SLOTS_TO_LOCK && hash_bucket_lock_offset <=
					CLUSTER_SIZE)
				qf_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK-1]);
			return false;
		}
		if (!qf_spin_lock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1],
											spin)) {
			qf_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
			if (hash_bucket_index >= NUM_SLOTS_TO_LOCK && hash_bucket_lock_offset <=
					CLUSTER_SIZE)
				qf_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK-1]);
			return false;
		}
#endif
	}
	return true;
}

static void qf_unlock(QF *cf, uint64_t hash_bucket_index, bool flag)
{
	uint64_t hash_bucket_lock_offset  = hash_bucket_index % NUM_SLOTS_TO_LOCK;
	if (flag) {
		if (NUM_SLOTS_TO_LOCK - hash_bucket_lock_offset <= CLUSTER_SIZE) {
			qf_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1]);
		}
		qf_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
	} else {
		qf_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK+1]);
		qf_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK]);
		if (hash_bucket_index >= NUM_SLOTS_TO_LOCK && hash_bucket_lock_offset <=
				CLUSTER_SIZE)
			qf_spin_unlock(&cf->mem->locks[hash_bucket_index/NUM_SLOTS_TO_LOCK-1]);
	}
}

static void modify_metadata(QF *cf, uint64_t *metadata, int cnt)
{
#ifdef LOG_WAIT_TIME
	qf_spin_lock(cf, &cf->mem->metadata_lock,cf->metadata->num_locks, true);
#else
	qf_spin_lock(&cf->mem->metadata_lock, true);
#endif
	*metadata = *metadata + cnt;
	qf_spin_unlock(&cf->mem->metadata_lock);
	return;
}

static inline int popcnt(uint64_t val)
{
	asm("popcnt %[val], %[val]"
			: [val] "+r" (val)
			:
			: "cc");
	return val;
}

static inline int64_t bitscanreverse(uint64_t val)
{
	if (val == 0) {
		return -1;
	} else {
		asm("bsr %[val], %[val]"
				: [val] "+r" (val)
				:
				: "cc");
		return val;
	}
}

static inline int popcntv(const uint64_t val, int ignore)
{
	if (ignore % 64)
		return popcnt (val & ~BITMASK(ignore % 64));
	else
		return popcnt(val);
}

// Returns the number of 1s up to (and including) the pos'th bit
// Bits are numbered from 0
static inline int bitrank(uint64_t val, int pos) {
	val = val & ((2ULL << pos) - 1);
	asm("popcnt %[val], %[val]"
			: [val] "+r" (val)
			:
			: "cc");
	return val;
}

/**
 * Returns the position of the k-th 1 in the 64-bit word x.
 * k is 0-based, so k=0 returns the position of the first 1.
 *
 * Uses the broadword selection algorithm by Vigna [1], improved by Gog
 * and Petri [2] and Vigna [3].
 *
 * [1] Sebastiano Vigna. Broadword Implementation of Rank/Select
 *    Queries. WEA, 2008
 *
 * [2] Simon Gog, Matthias Petri. Optimized succinct data
 * structures for massive data. Softw. Pract. Exper., 2014
 *
 * [3] Sebastiano Vigna. MG4J 5.2.1. http://mg4j.di.unimi.it/
 * The following code is taken from
 * https://github.com/facebook/folly/blob/b28186247104f8b90cfbe094d289c91f9e413317/folly/experimental/Select64.h
 */
const uint8_t kSelectInByte[2048] = {
	8, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0,
	1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0,
	2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0,
	1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0,
	3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 7, 0,
	1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0,
	2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0,
	1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
	4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0,
	1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 8, 8, 8, 1,
	8, 2, 2, 1, 8, 3, 3, 1, 3, 2, 2, 1, 8, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2,
	2, 1, 8, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1, 3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2, 1,
	4, 3, 3, 1, 3, 2, 2, 1, 8, 6, 6, 1, 6, 2, 2, 1, 6, 3, 3, 1, 3, 2, 2, 1, 6, 4,
	4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 6, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1,
	3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 8, 7, 7, 1, 7, 2,
	2, 1, 7, 3, 3, 1, 3, 2, 2, 1, 7, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1,
	7, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1, 3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2, 1, 4, 3,
	3, 1, 3, 2, 2, 1, 7, 6, 6, 1, 6, 2, 2, 1, 6, 3, 3, 1, 3, 2, 2, 1, 6, 4, 4, 1,
	4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 6, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1, 3, 2,
	2, 1, 5, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 8, 8, 8, 8, 8, 8, 8, 2,
	8, 8, 8, 3, 8, 3, 3, 2, 8, 8, 8, 4, 8, 4, 4, 2, 8, 4, 4, 3, 4, 3, 3, 2, 8, 8,
	8, 5, 8, 5, 5, 2, 8, 5, 5, 3, 5, 3, 3, 2, 8, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3,
	4, 3, 3, 2, 8, 8, 8, 6, 8, 6, 6, 2, 8, 6, 6, 3, 6, 3, 3, 2, 8, 6, 6, 4, 6, 4,
	4, 2, 6, 4, 4, 3, 4, 3, 3, 2, 8, 6, 6, 5, 6, 5, 5, 2, 6, 5, 5, 3, 5, 3, 3, 2,
	6, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3, 3, 2, 8, 8, 8, 7, 8, 7, 7, 2, 8, 7,
	7, 3, 7, 3, 3, 2, 8, 7, 7, 4, 7, 4, 4, 2, 7, 4, 4, 3, 4, 3, 3, 2, 8, 7, 7, 5,
	7, 5, 5, 2, 7, 5, 5, 3, 5, 3, 3, 2, 7, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3,
	3, 2, 8, 7, 7, 6, 7, 6, 6, 2, 7, 6, 6, 3, 6, 3, 3, 2, 7, 6, 6, 4, 6, 4, 4, 2,
	6, 4, 4, 3, 4, 3, 3, 2, 7, 6, 6, 5, 6, 5, 5, 2, 6, 5, 5, 3, 5, 3, 3, 2, 6, 5,
	5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3, 3, 2, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 3, 8, 8, 8, 8, 8, 8, 8, 4, 8, 8, 8, 4, 8, 4, 4, 3, 8, 8, 8, 8, 8, 8,
	8, 5, 8, 8, 8, 5, 8, 5, 5, 3, 8, 8, 8, 5, 8, 5, 5, 4, 8, 5, 5, 4, 5, 4, 4, 3,
	8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 6, 8, 6, 6, 3, 8, 8, 8, 6, 8, 6, 6, 4, 8, 6,
	6, 4, 6, 4, 4, 3, 8, 8, 8, 6, 8, 6, 6, 5, 8, 6, 6, 5, 6, 5, 5, 3, 8, 6, 6, 5,
	6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7,
	7, 3, 8, 8, 8, 7, 8, 7, 7, 4, 8, 7, 7, 4, 7, 4, 4, 3, 8, 8, 8, 7, 8, 7, 7, 5,
	8, 7, 7, 5, 7, 5, 5, 3, 8, 7, 7, 5, 7, 5, 5, 4, 7, 5, 5, 4, 5, 4, 4, 3, 8, 8,
	8, 7, 8, 7, 7, 6, 8, 7, 7, 6, 7, 6, 6, 3, 8, 7, 7, 6, 7, 6, 6, 4, 7, 6, 6, 4,
	6, 4, 4, 3, 8, 7, 7, 6, 7, 6, 6, 5, 7, 6, 6, 5, 6, 5, 5, 3, 7, 6, 6, 5, 6, 5,
	5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 5, 8, 8, 8, 8, 8, 8, 8, 5, 8, 8, 8, 5, 8, 5, 5, 4, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 6, 8, 6,
	6, 4, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 6, 8, 6, 6, 5, 8, 8, 8, 6, 8, 6, 6, 5,
	8, 6, 6, 5, 6, 5, 5, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8,
	8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 4, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7,
	8, 7, 7, 5, 8, 8, 8, 7, 8, 7, 7, 5, 8, 7, 7, 5, 7, 5, 5, 4, 8, 8, 8, 8, 8, 8,
	8, 7, 8, 8, 8, 7, 8, 7, 7, 6, 8, 8, 8, 7, 8, 7, 7, 6, 8, 7, 7, 6, 7, 6, 6, 4,
	8, 8, 8, 7, 8, 7, 7, 6, 8, 7, 7, 6, 7, 6, 6, 5, 8, 7, 7, 6, 7, 6, 6, 5, 7, 6,
	6, 5, 6, 5, 5, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 5, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 6,
	8, 6, 6, 5, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7,
	8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 5, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 6, 8, 8, 8, 8,
	8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 6, 8, 8, 8, 7, 8, 7, 7, 6, 8, 7, 7, 6, 7, 6,
	6, 5, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 6, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7
};

static inline uint64_t _select64(uint64_t x, int k)
{
	if (k >= popcnt(x)) { return 64; }
	
	const uint64_t kOnesStep4  = 0x1111111111111111ULL;
	const uint64_t kOnesStep8  = 0x0101010101010101ULL;
	const uint64_t kMSBsStep8  = 0x80ULL * kOnesStep8;

	uint64_t s = x;
	s = s - ((s & 0xA * kOnesStep4) >> 1);
	s = (s & 0x3 * kOnesStep4) + ((s >> 2) & 0x3 * kOnesStep4);
	s = (s + (s >> 4)) & 0xF * kOnesStep8;
	uint64_t byteSums = s * kOnesStep8;

	uint64_t kStep8 = k * kOnesStep8;
	uint64_t geqKStep8 = (((kStep8 | kMSBsStep8) - byteSums) & kMSBsStep8);
	uint64_t place = popcnt(geqKStep8) * 8;
	uint64_t byteRank = k - (((byteSums << 8) >> place) & (uint64_t)(0xFF));
	return place + kSelectInByte[((x >> place) & 0xFF) | (byteRank << 8)];
}

// Returns the position of the rank'th 1.  (rank = 0 returns the 1st 1)
// Returns 64 if there are fewer than rank+1 1s.
static inline uint64_t bitselect(uint64_t val, int rank) {
#ifdef __SSE4_2_
	uint64_t i = 1ULL << rank;
	asm("pdep %[val], %[mask], %[val]"
			: [val] "+r" (val)
			: [mask] "r" (i));
	asm("tzcnt %[bit], %[index]"
			: [index] "=r" (i)
			: [bit] "g" (val)
			: "cc");
	return i;
#endif
	return _select64(val, rank);
}

static inline uint64_t bitselectv(const uint64_t val, int ignore, int rank)
{
	return bitselect(val & ~BITMASK(ignore % 64), rank);
}

#if BITS_PER_SLOT > 0
/*static inline*/ qfblock * get_block(const QF *qf, uint64_t block_index)
{
	return &qf->blocks[block_index];
}
#else
/*static inline*/ qfblock * get_block(const QF *qf, uint64_t block_index)
{
	return (qfblock *)(((char *)qf->blocks) + block_index * (sizeof(qfblock) +
						SLOTS_PER_BLOCK * qf->metadata->bits_per_slot / 8));
}
#endif

static inline int is_runend(const QF *qf, uint64_t index)
{
	return (METADATA_WORD(qf, runends, index) >> ((index % SLOTS_PER_BLOCK) %
																								64)) & 1ULL;
}

static inline int is_occupied(const QF *qf, uint64_t index)
{
	return (METADATA_WORD(qf, occupieds, index) >> ((index % SLOTS_PER_BLOCK) %
																									64)) & 1ULL;
}

#if BITS_PER_SLOT == 8 || BITS_PER_SLOT == 16 || BITS_PER_SLOT == 32 || BITS_PER_SLOT == 64

static inline uint64_t get_slot(const QF *qf, uint64_t index)
{
	assert(index < qf->metadata->xnslots);
	return get_block(qf, index / SLOTS_PER_BLOCK)->slots[index % SLOTS_PER_BLOCK];
}

static inline void set_slot(const QF *qf, uint64_t index, uint64_t value)
{
	assert(index < qf->metadata->xnslots);
	get_block(qf, index / SLOTS_PER_BLOCK)->slots[index % SLOTS_PER_BLOCK] =
		value & BITMASK(qf->metadata->bits_per_slot);
}

#elif BITS_PER_SLOT > 0

/* Little-endian code ....  Big-endian is TODO */

static inline uint64_t get_slot(const QF *qf, uint64_t index)
{
	/* Should use __uint128_t to support up to 64-bit remainders, but gcc seems
	 * to generate buggy code.  :/  */
	assert(index < qf->metadata->xnslots);
	uint64_t *p = (uint64_t *)&get_block(qf, index /
																			 SLOTS_PER_BLOCK)->slots[(index %
																																SLOTS_PER_BLOCK)
																			 * BITS_PER_SLOT / 8];
	return (uint64_t)(((*p) >> (((index % SLOTS_PER_BLOCK) * BITS_PER_SLOT) %
															8)) & BITMASK(BITS_PER_SLOT));
}

static inline void set_slot(const QF *qf, uint64_t index, uint64_t value)
{
	/* Should use __uint128_t to support up to 64-bit remainders, but gcc seems
	 * to generate buggy code.  :/  */
	assert(index < qf->metadata->xnslots);
	uint64_t *p = (uint64_t *)&get_block(qf, index /
																			 SLOTS_PER_BLOCK)->slots[(index %
																																SLOTS_PER_BLOCK)
																			 * BITS_PER_SLOT / 8];
	uint64_t t = *p;
	uint64_t mask = BITMASK(BITS_PER_SLOT);
	uint64_t v = value;
	int shift = ((index % SLOTS_PER_BLOCK) * BITS_PER_SLOT) % 8;
	mask <<= shift;
	v <<= shift;
	t &= ~mask;
	t |= v;
	*p = t;
}

#else

/* Little-endian code ....  Big-endian is TODO */

static inline uint64_t get_slot(const QF *qf, uint64_t index)
{
	assert(index < qf->metadata->xnslots);
	/* Should use __uint128_t to support up to 64-bit remainders, but gcc seems
	 * to generate buggy code.  :/  */
	uint64_t *p = (uint64_t *)&get_block(qf, index /
																			 SLOTS_PER_BLOCK)->slots[(index %
																																SLOTS_PER_BLOCK)
																			 * qf->metadata->bits_per_slot / 8];
	return (uint64_t)(((*p) >> (((index % SLOTS_PER_BLOCK) *
															 qf->metadata->bits_per_slot) % 8)) &
										BITMASK(qf->metadata->bits_per_slot));
}

static inline void set_slot(const QF *qf, uint64_t index, uint64_t value)
{
	assert(index < qf->metadata->xnslots);
	/* Should use __uint128_t to support up to 64-bit remainders, but gcc seems
	 * to generate buggy code.  :/  */
	uint64_t *p = (uint64_t *)&get_block(qf, index /
																			 SLOTS_PER_BLOCK)->slots[(index %
																																SLOTS_PER_BLOCK)
																			 * qf->metadata->bits_per_slot / 8];
	uint64_t t = *p;
	uint64_t mask = BITMASK(qf->metadata->bits_per_slot);
	uint64_t v = value;
	int shift = ((index % SLOTS_PER_BLOCK) * qf->metadata->bits_per_slot) % 8;
	mask <<= shift;
	v <<= shift;
	t &= ~mask;
	t |= v;
	*p = t;
}

#endif

static inline uint64_t run_end(const QF *qf, uint64_t hash_bucket_index);

/*static inline*/ uint64_t block_offset(const QF *qf, uint64_t blockidx)
{
	/* If we have extended counters and a 16-bit (or larger) offset
		 field, then we can safely ignore the possibility of overflowing
		 that field. */
	if (sizeof(qf->blocks[0].offset) > 1 || 
			get_block(qf, blockidx)->offset < BITMASK(8*sizeof(qf->blocks[0].offset)))
		return get_block(qf, blockidx)->offset;

	return run_end(qf, SLOTS_PER_BLOCK * blockidx - 1) - SLOTS_PER_BLOCK *
		blockidx + 1;
}

//CHRISTINA
//uint64_t b_offset(const QF* qf, uint64_t blockidx){
//	return get_block(qf, blockidx)->offset;
//}

uint64_t run_end_strict(const QF *qf, uint64_t hash_bucket_index);
uint64_t block_offset_strict(const QF* qf, uint64_t blockidx){	
	return run_end(qf, SLOTS_PER_BLOCK * blockidx - 1) - SLOTS_PER_BLOCK * blockidx + 1;
}
uint64_t run_end_strict(const QF *qf, uint64_t hash_bucket_index)
{
	uint64_t bucket_block_index       = hash_bucket_index / SLOTS_PER_BLOCK;
	uint64_t bucket_intrablock_offset = hash_bucket_index % SLOTS_PER_BLOCK;
	uint64_t bucket_blocks_offset = block_offset_strict(qf, bucket_block_index);

	uint64_t bucket_intrablock_rank   = bitrank(get_block(qf,
																				bucket_block_index)->occupieds[0],
																				bucket_intrablock_offset);

	if (bucket_intrablock_rank == 0) {
		if (bucket_blocks_offset <= bucket_intrablock_offset)
			return hash_bucket_index;
		else
			return SLOTS_PER_BLOCK * bucket_block_index + bucket_blocks_offset - 1;
	}

	uint64_t runend_block_index  = bucket_block_index + bucket_blocks_offset /
		SLOTS_PER_BLOCK;
	uint64_t runend_ignore_bits  = bucket_blocks_offset % SLOTS_PER_BLOCK;
	uint64_t runend_rank         = bucket_intrablock_rank - 1;
	uint64_t runend_block_offset = bitselectv(get_block(qf,
																						runend_block_index)->runends[0],
																						runend_ignore_bits, runend_rank);
	if (runend_block_offset == SLOTS_PER_BLOCK) {
		if (bucket_blocks_offset == 0 && bucket_intrablock_rank == 0) {
			/* The block begins in empty space, and this bucket is in that region of
			 * empty space */
			return hash_bucket_index;
		} else {
			do {
				runend_rank        -= popcntv(get_block(qf,
																								runend_block_index)->runends[0],
																			runend_ignore_bits);
				runend_block_index++;
				runend_ignore_bits  = 0;
				runend_block_offset = bitselectv(get_block(qf,
																									 runend_block_index)->runends[0],
																				 runend_ignore_bits, runend_rank);
			} while (runend_block_offset == SLOTS_PER_BLOCK);
		}
	}

	uint64_t runend_index = SLOTS_PER_BLOCK * runend_block_index +
		runend_block_offset;
	if (runend_index < hash_bucket_index)
		return hash_bucket_index;
	else
		return runend_index;
}
//CHRISTINA END


static inline uint64_t run_end(const QF *qf, uint64_t hash_bucket_index)
{
	uint64_t bucket_block_index       = hash_bucket_index / SLOTS_PER_BLOCK;
	uint64_t bucket_intrablock_offset = hash_bucket_index % SLOTS_PER_BLOCK;
	uint64_t bucket_blocks_offset = block_offset(qf, bucket_block_index);

	uint64_t bucket_intrablock_rank   = bitrank(get_block(qf,
																				bucket_block_index)->occupieds[0],
																				bucket_intrablock_offset);

	if (bucket_intrablock_rank == 0) {
		if (bucket_blocks_offset <= bucket_intrablock_offset)
			return hash_bucket_index;
		else
			return SLOTS_PER_BLOCK * bucket_block_index + bucket_blocks_offset - 1;
	}

	uint64_t runend_block_index  = bucket_block_index + bucket_blocks_offset /
		SLOTS_PER_BLOCK;
	uint64_t runend_ignore_bits  = bucket_blocks_offset % SLOTS_PER_BLOCK;
	uint64_t runend_rank         = bucket_intrablock_rank - 1;
	uint64_t runend_block_offset = bitselectv(get_block(qf,
																						runend_block_index)->runends[0],
																						runend_ignore_bits, runend_rank);
	if (runend_block_offset == SLOTS_PER_BLOCK) {
		if (bucket_blocks_offset == 0 && bucket_intrablock_rank == 0) {
			/* The block begins in empty space, and this bucket is in that region of
			 * empty space */
			return hash_bucket_index;
		} else {
			do {
				runend_rank        -= popcntv(get_block(qf,
																								runend_block_index)->runends[0],
																			runend_ignore_bits);
				runend_block_index++;
				runend_ignore_bits  = 0;
				runend_block_offset = bitselectv(get_block(qf,
																									 runend_block_index)->runends[0],
																				 runend_ignore_bits, runend_rank);
			} while (runend_block_offset == SLOTS_PER_BLOCK);
		}
	}

	uint64_t runend_index = SLOTS_PER_BLOCK * runend_block_index +
		runend_block_offset;
	if (runend_index < hash_bucket_index)
		return hash_bucket_index;
	else
		return runend_index;
}

static inline int offset_lower_bound(const QF *qf, uint64_t slot_index)
{
	const qfblock * b = get_block(qf, slot_index / SLOTS_PER_BLOCK);
	const uint64_t slot_offset = slot_index % SLOTS_PER_BLOCK;
	const uint64_t boffset = b->offset;
	const uint64_t occupieds = b->occupieds[0] & BITMASK(slot_offset+1);
	assert(SLOTS_PER_BLOCK == 64);
	if (boffset <= slot_offset) {
		const uint64_t runends = (b->runends[0] & BITMASK(slot_offset)) >> boffset;
		return popcnt(occupieds) - popcnt(runends);
	}
	return boffset - slot_offset + popcnt(occupieds);
}

static inline int is_empty(const QF *qf, uint64_t slot_index)
{
	return offset_lower_bound(qf, slot_index) == 0;
}

static inline int might_be_empty(const QF *qf, uint64_t slot_index)
{
	return !is_occupied(qf, slot_index)
		&& !is_runend(qf, slot_index);
}

static inline int probably_is_empty(const QF *qf, uint64_t slot_index)
{
	return get_slot(qf, slot_index) == 0
		&& !is_occupied(qf, slot_index)
		&& !is_runend(qf, slot_index);
}

/*static inline*/ uint64_t find_first_empty_slot(const QF *qf, uint64_t from)
{
	do {
		int t = offset_lower_bound(qf, from);
		assert(t>=0);
		if (t == 0)
			break;
		from = from + t;
	} while(1);
	return from;
}

//CHRISTINA
uint64_t find_first_nonempty_slot(const QF* qf, uint64_t from){
	if(is_occupied(qf, from)){
		return from;
	}
	uint64_t block_index = from / SLOTS_PER_BLOCK;
	uint64_t rank = bitrank(get_block(qf, block_index)->occupieds[0],
			from % SLOTS_PER_BLOCK);
	uint64_t next_run = bitselect(get_block(qf,
				block_index)->occupieds[0],
			rank);
	if (next_run == 64) {
		rank = 0;
		while (next_run == 64 && block_index < qf->metadata->nblocks) {
			block_index++;
			next_run = bitselect(get_block(qf, block_index)->occupieds[0],
					rank);
		}
	}
	next_run = block_index*SLOTS_PER_BLOCK+next_run;
	if(block_index > qf->metadata->nblocks || next_run >= qf->metadata->xnslots){
		return qf->metadata->xnslots;
	}
	return next_run;
}

uint32_t increaseOffsetBy2(uint32_t offset){
	offset += 2;
	if(offset > BLOCK_OFFSET_MAX){
		return BLOCK_OFFSET_MAX;
	}else
		return offset;
}
uint32_t increaseOffsetBy1(uint32_t offset){
	offset++;
	if(offset > BLOCK_OFFSET_MAX){
		return BLOCK_OFFSET_MAX;
	}else
		return offset;
}
//CHRISTINA END


static inline uint64_t shift_into_b(const uint64_t a, const uint64_t b,
																		const int bstart, const int bend,
																		const int amount)
{
	const uint64_t a_component = bstart == 0 ? (a >> (64 - amount)) : 0;
	const uint64_t b_shifted_mask = BITMASK(bend - bstart) << bstart;
	const uint64_t b_shifted = ((b_shifted_mask & b) << amount) & b_shifted_mask;
	const uint64_t b_mask = ~b_shifted_mask;
	return a_component | b_shifted | (b & b_mask);
}

#if BITS_PER_SLOT == 8 || BITS_PER_SLOT == 16 || BITS_PER_SLOT == 32 || BITS_PER_SLOT == 64

static inline void shift_remainders(QF *qf, uint64_t start_index, uint64_t
																		empty_index)
{
	uint64_t start_block  = start_index / SLOTS_PER_BLOCK;
	uint64_t start_offset = start_index % SLOTS_PER_BLOCK;
	uint64_t empty_block  = empty_index / SLOTS_PER_BLOCK;
	uint64_t empty_offset = empty_index % SLOTS_PER_BLOCK;

	assert (start_index <= empty_index && empty_index < qf->metadata->xnslots);

	while (start_block < empty_block) {
		memmove(&get_block(qf, empty_block)->slots[1], 
						&get_block(qf, empty_block)->slots[0],
						empty_offset * sizeof(qf->blocks[0].slots[0]));
		get_block(qf, empty_block)->slots[0] = get_block(qf,
																			empty_block-1)->slots[SLOTS_PER_BLOCK-1];
		empty_block--;
		empty_offset = SLOTS_PER_BLOCK-1;
	}

	memmove(&get_block(qf, empty_block)->slots[start_offset+1], 
					&get_block(qf, empty_block)->slots[start_offset],
					(empty_offset - start_offset) * sizeof(qf->blocks[0].slots[0]));
}

static inline void shift_remainders_by2(QF *qf, uint64_t start_index, uint64_t
																		empty_index)
{
	uint64_t start_block  = start_index / SLOTS_PER_BLOCK;
	uint64_t start_offset = start_index % SLOTS_PER_BLOCK;
	uint64_t empty_block  = empty_index / SLOTS_PER_BLOCK;
	uint64_t empty_offset = empty_index % SLOTS_PER_BLOCK;

	assert (start_index <= empty_index && empty_index < qf->metadata->xnslots);

	if(empty_offset==0){
		get_block(qf, empty_block)->slots[0] = get_block(qf, empty_block-1)->slots[SLOTS_PER_BLOCK-2];
		empty_block--;
		empty_offset = SLOTS_PER_BLOCK;
	}
	while (start_block < empty_block) {
		memmove(&get_block(qf, empty_block)->slots[2], 
						&get_block(qf, empty_block)->slots[0],
						(empty_offset-2) * sizeof(qf->blocks[0].slots[0]));
		get_block(qf, empty_block)->slots[0] = get_block(qf,
																			empty_block-1)->slots[SLOTS_PER_BLOCK-2];
		get_block(qf, empty_block)->slots[1] = get_block(qf,
																			empty_block-1)->slots[SLOTS_PER_BLOCK-1];
		empty_block--;
		empty_offset = SLOTS_PER_BLOCK;
	}
	
	if(empty_offset > start_offset){
		memmove(&get_block(qf, empty_block)->slots[start_offset+2], 
					&get_block(qf, empty_block)->slots[start_offset],
					(SLOTS_PER_BLOCK - empty_offset+1) * sizeof(qf->blocks[0].slots[0]));
	}else if(start_offset == SLOTS_PER_BLOCK-1){//when start_offset == SLOTS_PER_BLOCK-1, the first slot in next block should be 0 instead
		get_block(qf, empty_block+1)->slots[0] = 0;
	}
}
#else

#define REMAINDER_WORD(qf, i) ((uint64_t *)&(get_block(qf, (i)/qf->metadata->bits_per_slot)->slots[8 * ((i) % qf->metadata->bits_per_slot)]))

static inline void shift_remainders(QF *qf, const uint64_t start_index, const
																		uint64_t empty_index)
{
	uint64_t last_word = (empty_index + 1) * qf->metadata->bits_per_slot / 64;
	const uint64_t first_word = start_index * qf->metadata->bits_per_slot / 64;
	int bend = ((empty_index + 1) * qf->metadata->bits_per_slot) % 64;
	const int bstart = (start_index * qf->metadata->bits_per_slot) % 64;

	while (last_word != first_word) {
		*REMAINDER_WORD(qf, last_word) = shift_into_b(*REMAINDER_WORD(qf, last_word-1),
																									*REMAINDER_WORD(qf, last_word),
																									0, bend, qf->metadata->bits_per_slot);
		last_word--;
		bend = 64;
	}
	*REMAINDER_WORD(qf, last_word) = shift_into_b(0, *REMAINDER_WORD(qf,
																																	 last_word),
																								bstart, bend,
																								qf->metadata->bits_per_slot);
}

static inline void shift_remainders_by2(QF *qf, const uint64_t start_index, const
																		uint64_t empty_index)
{
	uint64_t last_word = (empty_index + 1) * qf->metadata->bits_per_slot / 64;
	const uint64_t first_word = start_index * qf->metadata->bits_per_slot / 64;
	int bend = ((empty_index + 1) * qf->metadata->bits_per_slot) % 64;
	const int bstart = (start_index * qf->metadata->bits_per_slot) % 64;

	while (last_word != first_word) {
		*REMAINDER_WORD(qf, last_word) = shift_into_b(*REMAINDER_WORD(qf, last_word-1),
																									*REMAINDER_WORD(qf, last_word),
																									0, bend, qf->metadata->bits_per_slot*2);
		last_word--;
		bend = 64;
	}
	*REMAINDER_WORD(qf, last_word) = shift_into_b(0, *REMAINDER_WORD(qf, last_word), bstart, bend, qf->metadata->bits_per_slot*2);
	//the slot is not set to 0 when there is only one slot in last word.
}
#endif

/* static inline */ void qf_dump_block(const QF *qf, uint64_t i)
{
	uint64_t j;
	printf("#Block: %d\n", i);
	
	printf("offset: %-192d", get_block(qf, i)->offset);
	printf("\n");

	for (j = 0; j < SLOTS_PER_BLOCK; j++)
		printf("%02lx ", j);
	printf("\n");
	printf("occupieds\n");
	for (j = 0; j < SLOTS_PER_BLOCK; j++)
		printf(" %d ", (get_block(qf, i)->occupieds[j/64] & (1ULL << (j%64))) ? 1 : 0);
	printf("\n");
	
	printf("runends\n");
	for (j = 0; j < SLOTS_PER_BLOCK; j++)
		printf(" %d ", (get_block(qf, i)->runends[j/64] & (1ULL << (j%64))) ? 1 : 0);
	printf("\n");
	
	printf("slots\n");
#if BITS_PER_SLOT == 8 || BITS_PER_SLOT == 16 || BITS_PER_SLOT == 32
	for (j = 0; j < SLOTS_PER_BLOCK; j++)
		printf("%02x ", get_block(qf, i)->slots[j]);
#elif BITS_PER_SLOT == 64
	for (j = 0; j < SLOTS_PER_BLOCK; j++)
		printf("%02lx ", get_block(qf, i)->slots[j]);
#else
	for (j = 0; j < SLOTS_PER_BLOCK * qf->metadata->bits_per_slot / 8; j++)
		printf("%02x ", get_block(qf, i)->slots[j]);
#endif

	printf("\n");

	printf("\n");
}

void qf_dump(const QF *qf)
{
	uint64_t i;

	printf("%lu %lu %lu\n",
				 qf->metadata->nblocks,
				 qf->metadata->ndistinct_elts,
				 qf->metadata->nelts);

	for (i = 0; i < qf->metadata->nblocks; i++) {
		qf_dump_block(qf, i);
	}

}

static inline void find_next_n_empty_slots(QF *qf, uint64_t from, uint64_t n,
																					 uint64_t *indices)
{
	while (n) {
		indices[--n] = find_first_empty_slot(qf, from);
		from = indices[n] + 1;
	}
}

static inline void shift_slots(QF *qf, int64_t first, uint64_t last, uint64_t
															 distance)
{
	if(last<first)
		return;
	int64_t i;
	if (distance == 1)
		shift_remainders(qf, first, last+1);
	else
		for (i = last; i >= first; i--)
			set_slot(qf, i + distance, get_slot(qf, i));
}

static inline void shift_runends(QF *qf, int64_t first, uint64_t last,
																 uint64_t distance)
{
	if(last< first)
		return;
	assert(last < qf->metadata->xnslots && distance < 64);
	uint64_t first_word = first / 64;
	uint64_t bstart = first % 64;
	uint64_t last_word = (last + distance + 1) / 64;
	uint64_t bend = (last + distance + 1) % 64;

	if (last_word != first_word) {
		METADATA_WORD(qf, runends, 64*last_word) = shift_into_b(METADATA_WORD(qf, runends, 64*(last_word-1)), 
																														METADATA_WORD(qf, runends, 64*last_word),
																														0, bend, distance);
		bend = 64;
		last_word--;
		while (last_word != first_word) {
			METADATA_WORD(qf, runends, 64*last_word) = shift_into_b(METADATA_WORD(qf, runends, 64*(last_word-1)), 
																															METADATA_WORD(qf, runends, 64*last_word),
																															0, bend, distance);
			last_word--;
		}
		uint64_t offset=bend-bstart;
		if(offset <= distance){
			uint64_t over=distance-offset;
			if(over){
				//correct the next word
				METADATA_WORD(qf, runends, 64*(first_word+1)) &= ~((1ULL<<over)-1);
			}
			METADATA_WORD(qf, runends, 64*(first_word)) &= ((1ULL<<(64-offset))-1);
			return;
		}
	}
	METADATA_WORD(qf, runends, 64*last_word) = shift_into_b(0, METADATA_WORD(qf,runends,64*last_word),bstart, bend, distance);
}

static inline void insert_replace_slots_and_shift_remainders_and_runends_and_offsets(QF		*qf, 
																																										 int		 operation, 
																																										 uint64_t		 bucket_index,
																																										 uint64_t		 overwrite_index,
																																										 const uint64_t	*remainders, 
																																										 uint64_t		 total_remainders,
																																										 uint64_t		 noverwrites)
{
	uint64_t empties[67];
	uint64_t i;
	int64_t ninserts = total_remainders - noverwrites;
	uint64_t insert_index = overwrite_index + noverwrites;

	if (ninserts > 0) {
		/* First, shift things to create n empty spaces where we need them. */
		find_next_n_empty_slots(qf, insert_index, ninserts, empties);

		for (i = 0; i < ninserts - 1; i++)
			shift_slots(qf, empties[i+1] + 1, empties[i] - 1, i + 1);
		shift_slots(qf, insert_index, empties[ninserts - 1] - 1, ninserts);

		for (i = 0; i < ninserts - 1; i++)
			shift_runends(qf, empties[i+1] + 1, empties[i] - 1, i + 1);
		shift_runends(qf, insert_index, empties[ninserts - 1] - 1, ninserts);


		for (i = noverwrites; i < total_remainders - 1; i++)
			METADATA_WORD(qf, runends, overwrite_index + i) &= ~(1ULL <<
																													 (((overwrite_index
																															+ i) %
																														 SLOTS_PER_BLOCK)
																														% 64));

		switch (operation) {
			case 0: /* insert into empty bucket */
				assert (noverwrites == 0);
				METADATA_WORD(qf, runends, overwrite_index + total_remainders - 1) |=
					1ULL << (((overwrite_index + total_remainders - 1) %
										SLOTS_PER_BLOCK) % 64);
				break;
			case 1: /* append to bucket */
				METADATA_WORD(qf, runends, overwrite_index + noverwrites - 1)      &=
					~(1ULL << (((overwrite_index + noverwrites - 1) % SLOTS_PER_BLOCK) %
										 64));
				METADATA_WORD(qf, runends, overwrite_index + total_remainders - 1) |=
					1ULL << (((overwrite_index + total_remainders - 1) %
										SLOTS_PER_BLOCK) % 64);
				break;
			case 2: /* insert into bucket */
				METADATA_WORD(qf, runends, overwrite_index + total_remainders - 1) &=
					~(1ULL << (((overwrite_index + total_remainders - 1) %
											SLOTS_PER_BLOCK) % 64));
				break;
			default:
				fprintf(stderr, "Invalid operation %d\n", operation);
				abort();
		}

		uint64_t npreceding_empties = 0;
		for (i = bucket_index / SLOTS_PER_BLOCK + 1; i <= empties[0]/SLOTS_PER_BLOCK; i++) {
			while (npreceding_empties < ninserts &&
						 empties[ninserts - 1 - npreceding_empties]  / SLOTS_PER_BLOCK < i)
				npreceding_empties++;
			
			if (get_block(qf, i)->offset + ninserts - npreceding_empties < BITMASK(8*sizeof(qf->blocks[0].offset)))
				get_block(qf, i)->offset += ninserts - npreceding_empties;
			else
				get_block(qf, i)->offset = (uint8_t) BITMASK(8*sizeof(qf->blocks[0].offset));
		}
	}
						
	for (i = 0; i < total_remainders; i++)
		set_slot(qf, overwrite_index + i, remainders[i]);
	
	/*modify_metadata(qf, &qf->metadata->noccupied_slots, ninserts);*/
}

static inline void remove_replace_slots_and_shift_remainders_and_runends_and_offsets(QF		        *qf,
																																										 int		 operation,
																																										 uint64_t		 bucket_index,
																																										 uint64_t		 overwrite_index,
																																										 const uint64_t	*remainders,
																																										 uint64_t		 total_remainders,
																																										 uint64_t		 old_length)
{
	uint64_t i;

	// Update the slots
	for (i = 0; i < total_remainders; i++)
		set_slot(qf, overwrite_index + i, remainders[i]);

	// If this is the last thing in its run, then we may need to set a new runend bit
	if (is_runend(qf, overwrite_index + old_length - 1)) {
	  if (total_remainders > 0) { 
	    // If we're not deleting this entry entirely, then it will still the last entry in this run
	    METADATA_WORD(qf, runends, overwrite_index + total_remainders - 1) |= 1ULL << ((overwrite_index + total_remainders - 1) % 64);
	  } else if (overwrite_index > bucket_index &&
		     !is_runend(qf, overwrite_index - 1)) {
	    // If we're deleting this entry entirely, but it is not the first entry in this run,
	    // then set the preceding entry to be the runend
	    METADATA_WORD(qf, runends, overwrite_index - 1) |= 1ULL << ((overwrite_index - 1) % 64);
	  }
	}

	// shift slots back one run at a time
	uint64_t original_bucket = bucket_index;
	uint64_t current_bucket = bucket_index;
	uint64_t current_slot = overwrite_index + total_remainders;
	uint64_t current_distance = old_length - total_remainders;

	while (current_distance > 0) {
		if (is_runend(qf, current_slot + current_distance - 1)) {
			do {
				current_bucket++;
			} while (current_bucket < current_slot + current_distance &&
							 !is_occupied(qf, current_bucket));
		}

		if (current_bucket <= current_slot) {
			set_slot(qf, current_slot, get_slot(qf, current_slot + current_distance));
			if (is_runend(qf, current_slot) != 
					is_runend(qf, current_slot + current_distance))
				METADATA_WORD(qf, runends, current_slot) ^= 1ULL << (current_slot % 64);
			current_slot++;

		} else if (current_bucket <= current_slot + current_distance) {
			uint64_t i;
			for (i = current_slot; i < current_slot + current_distance; i++) {
				set_slot(qf, i, 0);
				METADATA_WORD(qf, runends, i) &= ~(1ULL << (i % 64));
			}

			current_distance = current_slot + current_distance - current_bucket;
			current_slot = current_bucket;
		} else {
			current_distance = 0;
		}
	}
	
	// reset the occupied bit of the hash bucket index if the hash is the
	// only item in the run and is removed completely.
	if (operation && !total_remainders)
		METADATA_WORD(qf, occupieds, bucket_index) &= ~(1ULL << (bucket_index % 64));
	
	// update the offset bits.
	// find the number of occupied slots in the original_bucket block.
	// Then find the runend slot corresponding to the last run in the
	// original_bucket block.
	// Update the offset of the block to which it belongs.
	uint64_t original_block = original_bucket / SLOTS_PER_BLOCK;
	while (1 && old_length > total_remainders) {	// we only update offsets if we shift/delete anything
		int32_t last_occupieds_bit = bitscanreverse(get_block(qf, original_block)->occupieds[0]);
		// there is nothing in the block
		if (last_occupieds_bit == -1) {
			if (get_block(qf, original_block + 1)->offset == 0)
				break;
			get_block(qf, original_block + 1)->offset = 0;
		} else {
			uint64_t last_occupieds_hash_index = SLOTS_PER_BLOCK * original_block + last_occupieds_bit;
			uint64_t runend_index = run_end(qf, last_occupieds_hash_index);
			// runend spans across the block
			// update the offset of the next block
			if (runend_index / SLOTS_PER_BLOCK == original_block) { // if the run ends in the same block
				if (get_block(qf, original_block + 1)->offset == 0)
					break;
				get_block(qf, original_block + 1)->offset = 0;
			} else if (runend_index / SLOTS_PER_BLOCK == original_block + 1) { // if the last run spans across one block
				if (get_block(qf, original_block + 1)->offset == (runend_index % SLOTS_PER_BLOCK) + 1)
					break;
				get_block(qf, original_block + 1)->offset = (runend_index % SLOTS_PER_BLOCK) + 1;
			} else { // if the last run spans across multiple blocks
				uint64_t i;
				for (i = original_block + 1; i < runend_index / SLOTS_PER_BLOCK - 1; i++)
					get_block(qf, i)->offset = SLOTS_PER_BLOCK;
				if (get_block(qf, runend_index / SLOTS_PER_BLOCK)->offset == (runend_index % SLOTS_PER_BLOCK) + 1)
					break;
				get_block(qf, runend_index / SLOTS_PER_BLOCK)->offset = (runend_index % SLOTS_PER_BLOCK) + 1;
			}
		}
		original_block++;
	}

	/*int num_slots_freed = old_length - total_remainders;*/
	/*modify_metadata(qf, &qf->metadata->noccupied_slots, -num_slots_freed);*/
	/*qf->metadata->noccupied_slots -= (old_length - total_remainders);*/
	if (!total_remainders) {
		/*modify_metadata(qf, &qf->metadata->ndistinct_elts, -1);*/
		/*qf->metadata->ndistinct_elts--;*/
	}
}

/*****************************************************************************
 * Code that uses the above to implement a QF with keys and inline counters. *
 *****************************************************************************/

/* 
	 Counter format:
	 0 xs:    <empty string>
	 1 x:     x
	 >=2 0s:    00c..cd  The highest bit in c is set, while the highest bit in d is not set.
	 >=2 xs:    xbc...cd For x != 0, b < x, the higest bit in c is set while the highest bit in d is not set.
	 */
/* static inline */uint64_t *encode_counter(QF *qf, uint64_t remainder, uint64_t
																			 counter, uint64_t *slots)
{
	uint64_t digit = remainder;
	uint64_t least_significant_bit_mask = 1ULL << (qf->metadata->bits_per_slot-1);
	uint64_t base = least_significant_bit_mask;
	uint64_t *p = slots;

	if (counter == 0)
		return p;

	if (counter == 1){
		*--p = remainder;
		return p;
	}

	counter--;
	digit = counter % base;
	counter /= base;
	*--p = digit;
	while(counter){
		digit = (counter % base) | least_significant_bit_mask;
		*--p = digit;
		counter /= base;
	}
	if(digit > remainder){
		*--p =0;
	}
	*--p = remainder;
	return p;
}

/* Returns the length of the encoding. 
REQUIRES: index points to first slot of a counter. */
static inline uint64_t decode_counter(const QF *qf, uint64_t index, uint64_t
																			*remainder, uint64_t *count)
{
	//uint64_t base;
	uint64_t rem;
	uint64_t cnt;
	uint64_t digit;
	uint64_t end;
	uint64_t least_significant_bit_mask = 1ULL << (qf->metadata->bits_per_slot-1);
	uint64_t base = least_significant_bit_mask;
	uint64_t counter_mask = base-1;

	*remainder = rem = get_slot(qf, index);

	if (is_runend(qf, index)) {
		*count = 1; 
		return index;
	}

	digit = get_slot(qf, index + 1);
	
	if(digit > rem){
		*count=1;
		return index;
	}else{
		cnt=0;
		end = index+1;
		if(digit==0){
			end++;
			digit = get_slot(qf, end);
		}
		while(digit & least_significant_bit_mask){
			cnt = cnt*base + (digit & counter_mask);
			end++;
			digit = get_slot(qf, end);
		}
		cnt = cnt*base + digit;
		*count = cnt+1;
		return end;
	}
}

/* return the next slot which corresponds to a 
 * different element 
 * */
static inline uint64_t next_slot(QF *qf, uint64_t current) 
{
	uint64_t rem = get_slot(qf, current);
	current++;

	while (get_slot(qf, current) == rem && current <= qf->metadata->nslots) {
		current++;
	}
	return current;
}

static inline bool insert1(QF *qf, __uint128_t hash, bool lock, bool spin)
{
	uint64_t hash_remainder           = hash & BITMASK(qf->metadata->bits_per_slot);
	uint64_t hash_bucket_index        = hash >> qf->metadata->bits_per_slot;
	uint64_t hash_bucket_block_offset = hash_bucket_index % SLOTS_PER_BLOCK;
	uint64_t least_significant_bit_mask = 1ULL << (qf->metadata->bits_per_slot-1);
	uint64_t counter_mask = least_significant_bit_mask - 1;
	uint64_t slot_mask = (1ULL<<(qf->metadata->bits_per_slot))-1;

	if (lock) {
		if (!qf_lock(qf, hash_bucket_index, spin, true))
			return false;
	}
	if (is_empty(qf, hash_bucket_index) /* might_be_empty(qf, hash_bucket_index) && runend_index == hash_bucket_index */) {
		METADATA_WORD(qf, runends, hash_bucket_index) |= 1ULL <<
			(hash_bucket_block_offset % 64);
		set_slot(qf, hash_bucket_index, hash_remainder);
		METADATA_WORD(qf, occupieds, hash_bucket_index) |= 1ULL <<
			(hash_bucket_block_offset % 64);

		/*modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);*/
		/*modify_metadata(qf, &qf->metadata->noccupied_slots, 1);*/
		/*modify_metadata(qf, &qf->metadata->nelts, 1);*/
	} else {
		uint64_t runend_index              = run_end(qf, hash_bucket_index);
		//int operation = 0; /* Insert into empty bucket */
		uint64_t insert_index = runend_index + 1;
		//uint64_t new_value = hash_remainder;

		/* printf("RUNSTART: %02lx RUNEND: %02lx\n", runstart_index, runend_index); */

		uint64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf,
																																	 hash_bucket_index
																																	 - 1) + 1;
		if(is_occupied(qf, hash_bucket_index)){
			uint64_t current_remainder = get_slot(qf, runstart_index);
			uint64_t tmp;
			//find where the remainder is
			while(current_remainder < hash_remainder && runstart_index < runend_index){
				tmp = get_slot(qf, runstart_index+1);
				if(tmp > current_remainder){
					//occuring only once
					current_remainder = tmp;
					runstart_index++;
					continue;
				}else{
					//occuring more than once
					if(tmp==0){
						runstart_index += 2;
						tmp = get_slot(qf, runstart_index);
					}else{
						runstart_index ++;
					}
					while(tmp & least_significant_bit_mask){
						runstart_index++;
						tmp = get_slot(qf, runstart_index);
					}
					runstart_index++;
					current_remainder = get_slot(qf, runstart_index);
				}
			}
			if(runstart_index == runend_index){
				if(current_remainder < hash_remainder){
					runstart_index++;
				}
			}
			if(runstart_index > runend_index){//largest new remainder
				insert_index = runend_index+1;
				uint64_t empty_slot_index;
				if(is_empty(qf, insert_index)){
					empty_slot_index = insert_index;
					METADATA_WORD(qf, runends, insert_index-1) &= ~(1ULL << (((insert_index-1) %SLOTS_PER_BLOCK) %64));
					METADATA_WORD(qf, runends, insert_index)   |= 1ULL << ((insert_index % SLOTS_PER_BLOCK) % 64);
				}else{
					empty_slot_index = find_first_empty_slot(qf, insert_index);
					shift_remainders(qf, insert_index, empty_slot_index);
					shift_runends(qf, insert_index-1, empty_slot_index-1, 1);
					METADATA_WORD(qf, runends, insert_index-1) &= ~(1ULL << (((insert_index-1) %SLOTS_PER_BLOCK) %64));

				}
				for (uint64_t i = hash_bucket_index / SLOTS_PER_BLOCK + 1; i <=
						 empty_slot_index/SLOTS_PER_BLOCK; i++) {
					if (get_block(qf, i)->offset < BITMASK(8*sizeof(qf->blocks[0].offset)))
						get_block(qf, i)->offset++;
					//assert(get_block(qf, i)->offset != 0);
				}
				set_slot(qf, insert_index, hash_remainder);
			}else{//remainder inserted within
				//extend counter of existing remainder
				if(current_remainder == hash_remainder){
					//locate the end of the counter
					uint64_t counter_end=runstart_index+1;
					if(runstart_index == runend_index){
						counter_end = runstart_index;
					}else{
						tmp=get_slot(qf, counter_end);
						if(tmp > hash_remainder){
							counter_end = runstart_index;
						}else{
							if(tmp == 0){
								counter_end++;
								tmp = get_slot(qf, counter_end);
							}
							while(tmp & least_significant_bit_mask){
								counter_end++;
								tmp = get_slot(qf, counter_end);
							}
						}
					}
					//update the counter
					if(runstart_index == counter_end){//counter being 1
						if(hash_remainder == 0){
							//shift to get two empty slots after insert_index
							insert_index = counter_end+1;
							uint64_t empty_slot_index1 = find_first_empty_slot(qf, insert_index);
							uint64_t empty_slot_index2 = find_first_empty_slot(qf, empty_slot_index1+1);
							if(empty_slot_index2 != empty_slot_index1+1){
								shift_remainders(qf, empty_slot_index1+1, empty_slot_index2);
								shift_runends(qf, empty_slot_index1+1, empty_slot_index2-1, 1);
								//update offset
								for(uint64_t i = empty_slot_index1/SLOTS_PER_BLOCK+1; i <= empty_slot_index2/SLOTS_PER_BLOCK; i++){
									if(get_block(qf, i)->offset < BITMASK(8*sizeof(qf->blocks[0].offset))){
										get_block(qf, i)->offset++;
									}
								}
							}	
							shift_remainders_by2(qf, insert_index, empty_slot_index1+1);
							set_slot(qf, insert_index, 0);
							set_slot(qf, insert_index+1, 1);
							//shift_runends(qf, insert_index, empty_slot_index1-1, 2);
							if(insert_index > runend_index){
								shift_runends(qf, insert_index-1, empty_slot_index1-1, 2);
								//METADATA_WORD(qf, runends, insert_index-1) &= ~(1ULL << (((insert_index-1)%SLOTS_PER_BLOCK)%64));
							}else{
								shift_runends(qf, insert_index, empty_slot_index1-1, 2);
								//METADATA_WORD(qf, runends, (insert_index+1)) &= ~(1ULL << (((insert_index+1)%SLOTS_PER_BLOCK)%64));
							}
							if((runend_index>insert_index?(runend_index-insert_index):(insert_index - runend_index)) < 2){
							//if(abs(runend_index-insert_index) < 2){
								METADATA_WORD(qf, runends, runend_index) &= ~(1ULL << ((runend_index%SLOTS_PER_BLOCK)%64));
							}

							//update offset
							for(uint64_t i = hash_bucket_index/SLOTS_PER_BLOCK+1; i <= empty_slot_index1/SLOTS_PER_BLOCK; i++){
								get_block(qf, i)->offset = increaseOffsetBy2(get_block(qf, i)->offset);
							}
							if(empty_slot_index1%SLOTS_PER_BLOCK == SLOTS_PER_BLOCK-1){
								uint64_t i = empty_slot_index1/SLOTS_PER_BLOCK+1;
								if(get_block(qf, i)->offset == 0){
									get_block(qf, i)->offset =1;
								}
							}
						}else{//non-zero remainder occurring once
							insert_index = counter_end + 1;
							uint64_t empty_slot_index = find_first_empty_slot(qf, insert_index);
							shift_remainders(qf, insert_index, empty_slot_index);

							set_slot(qf, insert_index, 1);
							if(insert_index > runend_index){
								shift_runends(qf, runend_index, empty_slot_index-1, 1);
								METADATA_WORD(qf, runends, runend_index) &= ~(1ULL << ((runend_index%SLOTS_PER_BLOCK)%64));
							}else{
								shift_runends(qf, insert_index, empty_slot_index-1, 1);
							}
							
							for (uint64_t i = hash_bucket_index / SLOTS_PER_BLOCK + 1; i <=
									 empty_slot_index/SLOTS_PER_BLOCK; i++) {
								if (get_block(qf, i)->offset < BITMASK(8*sizeof(qf->blocks[0].offset)))
									get_block(qf, i)->offset++;
								//assert(get_block(qf, i)->offset != 0);
							}
						}
					}else{//counter being >1
						uint64_t zeroStart = get_slot(qf, runstart_index+1)==0?1:0;//whether counter starts with zero
						//increase the counter from lowest slot to highest slot
						uint64_t digit;
						bool advance=true;
						digit = get_slot(qf, counter_end);
						if(digit == counter_mask){
							advance = true;
							set_slot(qf, counter_end, 0);
						}else{
							advance = false;
							set_slot(qf, counter_end, digit+1);
						}
						for(uint64_t i = counter_end-1; i > (runstart_index+zeroStart) && advance; i--){
							digit = get_slot(qf, i);
							if(digit == slot_mask){
								set_slot(qf, i, least_significant_bit_mask);
								advance= true;
							}else{
								set_slot(qf, i, digit+1);
								advance = false;
								break;
							}
						}
						if(advance){
							if(zeroStart){
								if(hash_remainder > least_significant_bit_mask){
									set_slot(qf, runstart_index+1, least_significant_bit_mask+1);
								}else{
									insert_index = runstart_index+2;
									uint64_t empty_slot_index = find_first_empty_slot(qf, insert_index);
									shift_remainders(qf, insert_index, empty_slot_index);

									set_slot(qf, insert_index, least_significant_bit_mask+1);
									shift_runends(qf, insert_index, empty_slot_index-1, 1);
									
									for (uint64_t i = hash_bucket_index / SLOTS_PER_BLOCK + 1; i <=
											 empty_slot_index/SLOTS_PER_BLOCK; i++) {
										if (get_block(qf, i)->offset < BITMASK(8*sizeof(qf->blocks[0].offset)))
											get_block(qf, i)->offset++;
										//assert(get_block(qf, i)->offset != 0);
									}
									if(insert_index == runend_index)
										METADATA_WORD(qf, runends, runend_index) &= ~(1ULL << ((runend_index%SLOTS_PER_BLOCK)%64));
								}
							}else{
								insert_index = runstart_index+1;		
								uint64_t empty_slot_index = find_first_empty_slot(qf, insert_index);
								shift_remainders(qf, insert_index, empty_slot_index);

								set_slot(qf, insert_index, least_significant_bit_mask+1);
								shift_runends(qf, insert_index, empty_slot_index-1, 1);
								
								for (uint64_t i = hash_bucket_index / SLOTS_PER_BLOCK + 1; i <=
										 empty_slot_index/SLOTS_PER_BLOCK; i++) {
									if (get_block(qf, i)->offset < BITMASK(8*sizeof(qf->blocks[0].offset)))
										get_block(qf, i)->offset++;
									//assert(get_block(qf, i)->offset != 0);
								}
								if(insert_index == runend_index)
									METADATA_WORD(qf, runends, runend_index) &= ~(1ULL << ((runend_index%SLOTS_PER_BLOCK)%64));
							}
						}
						if(get_slot(qf, runstart_index+1) > hash_remainder){
							insert_index = runstart_index+1;
							uint64_t empty_slot_index = find_first_empty_slot(qf, insert_index);
							shift_remainders(qf, insert_index, empty_slot_index);

							//set_slot(qf, insert_index, least_significant_bit_mask+1);
							set_slot(qf, insert_index, 0);
							shift_runends(qf, insert_index, empty_slot_index-1, 1);
							
							for (uint64_t i = hash_bucket_index / SLOTS_PER_BLOCK + 1; i <=
									 empty_slot_index/SLOTS_PER_BLOCK; i++) {
								if (get_block(qf, i)->offset < BITMASK(8*sizeof(qf->blocks[0].offset)))
									get_block(qf, i)->offset++;
								//assert(get_block(qf, i)->offset != 0);
							}
							if(insert_index == runend_index)
								METADATA_WORD(qf, runends, runend_index) &= ~(1ULL << ((runend_index%SLOTS_PER_BLOCK)%64));
						}
					}
				}else{//add a new remainder
					insert_index = runstart_index;
					uint64_t empty_slot_index = find_first_empty_slot(qf, runend_index+1);
					shift_remainders(qf, insert_index, empty_slot_index);

					set_slot(qf, insert_index, hash_remainder);
					shift_runends(qf, insert_index, empty_slot_index-1, 1);
					
					for (uint64_t i = hash_bucket_index / SLOTS_PER_BLOCK + 1; i <=
							 empty_slot_index/SLOTS_PER_BLOCK; i++) {
						if (get_block(qf, i)->offset < BITMASK(8*sizeof(qf->blocks[0].offset)))
							get_block(qf, i)->offset++;
						//assert(get_block(qf, i)->offset != 0);
					}
				}
			}
		}else{
			insert_index = runstart_index;
			uint64_t empty_slot_index = find_first_empty_slot(qf, runend_index+1);
			shift_remainders(qf, insert_index, empty_slot_index);

			set_slot(qf, insert_index, hash_remainder);
			shift_runends(qf, insert_index, empty_slot_index-1, 1);
			
			for (uint64_t i = hash_bucket_index / SLOTS_PER_BLOCK + 1; i <=
					 empty_slot_index/SLOTS_PER_BLOCK; i++) {
				if (get_block(qf, i)->offset < BITMASK(8*sizeof(qf->blocks[0].offset)))
					get_block(qf, i)->offset++;
				//assert(get_block(qf, i)->offset != 0);
			}
			
			METADATA_WORD(qf, runends, insert_index) |= 1ULL <<
				(insert_index % 64);
			METADATA_WORD(qf, occupieds, hash_bucket_index) |= 1ULL <<
				(hash_bucket_block_offset % 64);
		}
	}

	if (lock) {
		qf_unlock(qf, hash_bucket_index, true);
	}

	return true;
}

static inline bool insert1_advance(QF *qf, __uint128_t hash, bool lock, bool spin, bool& isNew)
{
	uint64_t hash_remainder           = hash & BITMASK(qf->metadata->bits_per_slot);
	uint64_t hash_bucket_index        = hash >> qf->metadata->bits_per_slot;
	uint64_t hash_bucket_block_offset = hash_bucket_index % SLOTS_PER_BLOCK;
	uint64_t least_significant_bit_mask = 1ULL << (qf->metadata->bits_per_slot-1);
	uint64_t counter_mask = least_significant_bit_mask - 1;
	uint64_t slot_mask = (1ULL<<(qf->metadata->bits_per_slot))-1;

	if (lock) {
		if (!qf_lock(qf, hash_bucket_index, spin, true))
			return false;
	}
	if (is_empty(qf, hash_bucket_index) /* might_be_empty(qf, hash_bucket_index) && runend_index == hash_bucket_index */) {
		isNew = true;
		METADATA_WORD(qf, runends, hash_bucket_index) |= 1ULL <<
			(hash_bucket_block_offset % 64);
		set_slot(qf, hash_bucket_index, hash_remainder);
		METADATA_WORD(qf, occupieds, hash_bucket_index) |= 1ULL <<
			(hash_bucket_block_offset % 64);
		
		/*modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);*/
		/*modify_metadata(qf, &qf->metadata->noccupied_slots, 1);*/
		/*modify_metadata(qf, &qf->metadata->nelts, 1);*/
	} else {
		uint64_t runend_index              = run_end(qf, hash_bucket_index);
		//int operation = 0; /* Insert into empty bucket */
		uint64_t insert_index = runend_index + 1;
		//uint64_t new_value = hash_remainder;

		/* printf("RUNSTART: %02lx RUNEND: %02lx\n", runstart_index, runend_index); */

		uint64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf,
																																	 hash_bucket_index
																																	 - 1) + 1;
		if(is_occupied(qf, hash_bucket_index)){
			uint64_t current_remainder = get_slot(qf, runstart_index);
			uint64_t tmp;
			//find where the remainder is
			while(current_remainder < hash_remainder && runstart_index < runend_index){
				tmp = get_slot(qf, runstart_index+1);
				if(tmp > current_remainder){
					//occuring only once
					current_remainder = tmp;
					runstart_index++;
					continue;
				}else{
					//occuring more than once
					if(tmp==0){
						runstart_index += 2;
						tmp = get_slot(qf, runstart_index);
					}else{
						runstart_index ++;
					}
					while(tmp & least_significant_bit_mask){
						runstart_index++;
						tmp = get_slot(qf, runstart_index);
					}
					runstart_index++;
					current_remainder = get_slot(qf, runstart_index);
				}
			}
			if(runstart_index == runend_index){
				if(current_remainder < hash_remainder){
					runstart_index++;
				}
			}
			if(runstart_index > runend_index){//largest new remainder
				isNew = true;
				insert_index = runend_index+1;
				uint64_t empty_slot_index;
				if(is_empty(qf, insert_index)){
					empty_slot_index = insert_index;
					METADATA_WORD(qf, runends, insert_index-1) &= ~(1ULL << (((insert_index-1) %SLOTS_PER_BLOCK) %64));
					METADATA_WORD(qf, runends, insert_index)   |= 1ULL << ((insert_index % SLOTS_PER_BLOCK) % 64);
				}else{
					empty_slot_index = find_first_empty_slot(qf, insert_index);
					shift_remainders(qf, insert_index, empty_slot_index);
					shift_runends(qf, insert_index-1, empty_slot_index-1, 1);
					METADATA_WORD(qf, runends, insert_index-1) &= ~(1ULL << (((insert_index-1) %SLOTS_PER_BLOCK) %64));
				}
				for (uint64_t i = hash_bucket_index / SLOTS_PER_BLOCK + 1; i <=
						 empty_slot_index/SLOTS_PER_BLOCK; i++) {
					if (get_block(qf, i)->offset < BITMASK(8*sizeof(qf->blocks[0].offset)))
						get_block(qf, i)->offset++;
					//assert(get_block(qf, i)->offset != 0);
				}
				set_slot(qf, insert_index, hash_remainder);
			}else{//remainder inserted within
				//extend counter of existing remainder
				if(current_remainder == hash_remainder){
					isNew = false;
					//locate the end of the counter
					uint64_t counter_end=runstart_index+1;
					if(runstart_index == runend_index){
						counter_end = runstart_index;
					}else{
						tmp=get_slot(qf, counter_end);
						if(tmp > hash_remainder){
							counter_end = runstart_index;
						}else{
							if(tmp == 0){
								counter_end++;
								tmp = get_slot(qf, counter_end);
							}
							while(tmp & least_significant_bit_mask){
								counter_end++;
								tmp = get_slot(qf, counter_end);
							}
						}
					}
					//update the counter
					if(runstart_index == counter_end){//counter being 1
						if(hash_remainder == 0){
							//shift to get two empty slots after insert_index
							insert_index = counter_end+1;
							uint64_t empty_slot_index1 = find_first_empty_slot(qf, insert_index);
							uint64_t empty_slot_index2 = find_first_empty_slot(qf, empty_slot_index1+1);
							if(empty_slot_index2 != empty_slot_index1+1){
								shift_remainders(qf, empty_slot_index1+1, empty_slot_index2);
								shift_runends(qf, empty_slot_index1+1, empty_slot_index2-1, 1);
								//update offset
								for(uint64_t i = empty_slot_index1/SLOTS_PER_BLOCK+1; i <= empty_slot_index2/SLOTS_PER_BLOCK; i++){
									if(get_block(qf, i)->offset < BITMASK(8*sizeof(qf->blocks[0].offset))){
										get_block(qf, i)->offset++;
									}
								}
							}	
							shift_remainders_by2(qf, insert_index, empty_slot_index1+1);
							set_slot(qf, insert_index, 0);
							set_slot(qf, insert_index+1, 1);
							//shift_runends(qf, insert_index, empty_slot_index1-1, 2);
							if(insert_index > runend_index){
								shift_runends(qf, insert_index-1, empty_slot_index1-1, 2);
								//METADATA_WORD(qf, runends, insert_index-1) &= ~(1ULL << (((insert_index-1)%SLOTS_PER_BLOCK)%64));
							}else{
								shift_runends(qf, insert_index, empty_slot_index1-1, 2);
								//METADATA_WORD(qf, runends, (insert_index+1)) &= ~(1ULL << (((insert_index+1)%SLOTS_PER_BLOCK)%64));
							}
							//if(abs(runend_index-insert_index) < 2){
							if((runend_index>insert_index?(runend_index-insert_index):(insert_index - runend_index)) < 2){
								METADATA_WORD(qf, runends, runend_index) &= ~(1ULL << ((runend_index%SLOTS_PER_BLOCK)%64));
							}

							//update offset
							for(uint64_t i = hash_bucket_index/SLOTS_PER_BLOCK+1; i <= empty_slot_index1/SLOTS_PER_BLOCK; i++){
								get_block(qf, i)->offset = increaseOffsetBy2(get_block(qf, i)->offset);
							}
							if(empty_slot_index1%SLOTS_PER_BLOCK == SLOTS_PER_BLOCK-1){
								uint64_t i = empty_slot_index1/SLOTS_PER_BLOCK+1;
								if(get_block(qf, i)->offset == 0){
									get_block(qf, i)->offset =1;
								}
							}
						}else{//non-zero remainder occurring once
							insert_index = counter_end + 1;
							uint64_t empty_slot_index = find_first_empty_slot(qf, insert_index);
							shift_remainders(qf, insert_index, empty_slot_index);

							set_slot(qf, insert_index, 1);
							if(insert_index > runend_index){
								shift_runends(qf, runend_index, empty_slot_index-1, 1);
								METADATA_WORD(qf, runends, runend_index) &= ~(1ULL << ((runend_index%SLOTS_PER_BLOCK)%64));
							}else{
								shift_runends(qf, insert_index, empty_slot_index-1, 1);
							}
							
							for (uint64_t i = hash_bucket_index / SLOTS_PER_BLOCK + 1; i <=
									 empty_slot_index/SLOTS_PER_BLOCK; i++) {
								if (get_block(qf, i)->offset < BITMASK(8*sizeof(qf->blocks[0].offset)))
									get_block(qf, i)->offset++;
								//assert(get_block(qf, i)->offset != 0);
							}
						}
					}else{//counter being >1
						uint64_t zeroStart = get_slot(qf, runstart_index+1)==0?1:0;//whether counter starts with zero
						//increase the counter from lowest slot to highest slot
						uint64_t digit;
						bool advance=true;
						digit = get_slot(qf, counter_end);
						if(digit == counter_mask){
							advance = true;
							set_slot(qf, counter_end, 0);
						}else{
							advance = false;
							set_slot(qf, counter_end, digit+1);
						}
						for(uint64_t i = counter_end-1; i > (runstart_index+zeroStart) && advance; i--){
							digit = get_slot(qf, i);
							if(digit == slot_mask){
								set_slot(qf, i, least_significant_bit_mask);
								advance= true;
							}else{
								set_slot(qf, i, digit+1);
								advance = false;
								break;
							}
						}
						if(advance){
							if(zeroStart){
								if(hash_remainder > least_significant_bit_mask){
									set_slot(qf, runstart_index+1, least_significant_bit_mask+1);
								}else{
									insert_index = runstart_index+2;
									uint64_t empty_slot_index = find_first_empty_slot(qf, insert_index);
									shift_remainders(qf, insert_index, empty_slot_index);

									set_slot(qf, insert_index, least_significant_bit_mask+1);
									shift_runends(qf, insert_index, empty_slot_index-1, 1);
									
									for (uint64_t i = hash_bucket_index / SLOTS_PER_BLOCK + 1; i <=
											 empty_slot_index/SLOTS_PER_BLOCK; i++) {
										if (get_block(qf, i)->offset < BITMASK(8*sizeof(qf->blocks[0].offset)))
											get_block(qf, i)->offset++;
										//assert(get_block(qf, i)->offset != 0);
									}
									if(insert_index == runend_index)
										METADATA_WORD(qf, runends, runend_index) &= ~(1ULL << ((runend_index%SLOTS_PER_BLOCK)%64));
								}
							}else{
								insert_index = runstart_index+1;		
								uint64_t empty_slot_index = find_first_empty_slot(qf, insert_index);
								shift_remainders(qf, insert_index, empty_slot_index);

								set_slot(qf, insert_index, least_significant_bit_mask+1);
								shift_runends(qf, insert_index, empty_slot_index-1, 1);
								
								for (uint64_t i = hash_bucket_index / SLOTS_PER_BLOCK + 1; i <=
										 empty_slot_index/SLOTS_PER_BLOCK; i++) {
									if (get_block(qf, i)->offset < BITMASK(8*sizeof(qf->blocks[0].offset)))
										get_block(qf, i)->offset++;
									//assert(get_block(qf, i)->offset != 0);
								}
								if(insert_index == runend_index)
									METADATA_WORD(qf, runends, runend_index) &= ~(1ULL << ((runend_index%SLOTS_PER_BLOCK)%64));
							}
						}
						if(get_slot(qf, runstart_index+1) > hash_remainder){
							insert_index = runstart_index+1;
							uint64_t empty_slot_index = find_first_empty_slot(qf, insert_index);
							shift_remainders(qf, insert_index, empty_slot_index);

							//set_slot(qf, insert_index, least_significant_bit_mask+1);
							set_slot(qf, insert_index, 0);
							shift_runends(qf, insert_index, empty_slot_index-1, 1);
							
							for (uint64_t i = hash_bucket_index / SLOTS_PER_BLOCK + 1; i <=
									 empty_slot_index/SLOTS_PER_BLOCK; i++) {
								if (get_block(qf, i)->offset < BITMASK(8*sizeof(qf->blocks[0].offset)))
									get_block(qf, i)->offset++;
								//assert(get_block(qf, i)->offset != 0);
							}
							if(insert_index == runend_index)
								METADATA_WORD(qf, runends, runend_index) &= ~(1ULL << ((runend_index%SLOTS_PER_BLOCK)%64));
						}
					}
				}else{//add a new remainder
					isNew=true;
					insert_index = runstart_index;
					uint64_t empty_slot_index = find_first_empty_slot(qf, runend_index+1);
					shift_remainders(qf, insert_index, empty_slot_index);

					set_slot(qf, insert_index, hash_remainder);
					shift_runends(qf, insert_index, empty_slot_index-1, 1);
					
					for (uint64_t i = hash_bucket_index / SLOTS_PER_BLOCK + 1; i <=
							 empty_slot_index/SLOTS_PER_BLOCK; i++) {
						if (get_block(qf, i)->offset < BITMASK(8*sizeof(qf->blocks[0].offset)))
							get_block(qf, i)->offset++;
						//assert(get_block(qf, i)->offset != 0);
					}
				}
			}
		}else{
			isNew=true;
			insert_index = runstart_index;
			uint64_t empty_slot_index = find_first_empty_slot(qf, runend_index+1);
			if(insert_index != empty_slot_index){
				shift_remainders(qf, insert_index, empty_slot_index);
				shift_runends(qf, insert_index, empty_slot_index-1, 1);
			}
			set_slot(qf, insert_index, hash_remainder);
			for (uint64_t i = hash_bucket_index / SLOTS_PER_BLOCK + 1; i <=
					 empty_slot_index/SLOTS_PER_BLOCK; i++) {
				if (get_block(qf, i)->offset < BITMASK(8*sizeof(qf->blocks[0].offset)))
					get_block(qf, i)->offset++;
				//assert(get_block(qf, i)->offset != 0);
			}
			
			METADATA_WORD(qf, runends, insert_index) |= 1ULL <<
				(insert_index % 64);
			METADATA_WORD(qf, occupieds, hash_bucket_index) |= 1ULL <<
				(hash_bucket_block_offset % 64);
		}
	}

	if (lock) {
		qf_unlock(qf, hash_bucket_index, true);
	}

	return true;
}
static inline bool insert(QF *qf, __uint128_t hash, uint64_t count, bool lock,
													bool spin)
{
	uint64_t hash_remainder           = hash & BITMASK(qf->metadata->bits_per_slot);
	uint64_t hash_bucket_index        = hash >> qf->metadata->bits_per_slot;
	uint64_t hash_bucket_block_offset = hash_bucket_index % SLOTS_PER_BLOCK;
	/*uint64_t hash_bucket_lock_offset  = hash_bucket_index % NUM_SLOTS_TO_LOCK;*/

	if (lock) {
		if (!qf_lock(qf, hash_bucket_index, spin, false))
			return false;
	}

	uint64_t runend_index             = run_end(qf, hash_bucket_index);
	
	/* Empty slot */
	if (might_be_empty(qf, hash_bucket_index) && runend_index ==
			hash_bucket_index) {
		METADATA_WORD(qf, runends, hash_bucket_index) |= 1ULL <<
			(hash_bucket_block_offset % 64);
		set_slot(qf, hash_bucket_index, hash_remainder);
		METADATA_WORD(qf, occupieds, hash_bucket_index) |= 1ULL <<
			(hash_bucket_block_offset % 64);
		
		/*modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);*/
		/*modify_metadata(qf, &qf->metadata->noccupied_slots, 1);*/
		/*modify_metadata(qf, &qf->metadata->nelts, 1);*/
		/* This trick will, I hope, keep the fast case fast. */
		if (count > 1) {
			insert(qf, hash, count - 1, false, false);
		}
	} else { /* Non-empty slot */
		uint64_t new_values[67];
		int64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf,
																																	hash_bucket_index
																																	- 1) + 1;

		if (!is_occupied(qf, hash_bucket_index)) { /* Empty bucket, but its slot is occupied. */
			uint64_t *p = encode_counter(qf, hash_remainder, count, &new_values[67]);
			insert_replace_slots_and_shift_remainders_and_runends_and_offsets(qf, 
																																				0, 
																																				hash_bucket_index, 
																																				runstart_index, 
																																				p, 
																																				&new_values[67] - p, 
																																				0);
			/*modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);*/
		} else { /* Non-empty bucket */

			uint64_t current_remainder, current_count, current_end;

			/* Find the counter for this remainder, if one exists. */
			current_end = decode_counter(qf, runstart_index, &current_remainder,
																	 &current_count);
			while (current_remainder < hash_remainder && !is_runend(qf, current_end)) {
				runstart_index = current_end + 1;
				current_end = decode_counter(qf, runstart_index, &current_remainder,
																		 &current_count);	
			}

			/* If we reached the end of the run w/o finding a counter for this remainder,
				 then append a counter for this remainder to the run. */
			if (current_remainder < hash_remainder) {
				uint64_t *p = encode_counter(qf, hash_remainder, count, &new_values[67]);
				insert_replace_slots_and_shift_remainders_and_runends_and_offsets(qf, 
																																					1, /* Append to bucket */
																																					hash_bucket_index, 
																																					current_end + 1, 
																																					p, 
																																					&new_values[67] - p, 
																																					0);
				/*modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);*/
				/* Found a counter for this remainder.  Add in the new count. */
			} else if (current_remainder == hash_remainder) {
				uint64_t *p = encode_counter(qf, hash_remainder, current_count + count, &new_values[67]);
				insert_replace_slots_and_shift_remainders_and_runends_and_offsets(qf, 
																																					is_runend(qf, current_end) ? 1 : 2, 
																																					hash_bucket_index, 
																																					runstart_index, 
																																					p, 
																																					&new_values[67] - p, 
																																					current_end - runstart_index + 1);
				/* No counter for this remainder, but there are larger
					 remainders, so we're not appending to the bucket. */
			} else {
				uint64_t *p = encode_counter(qf, hash_remainder, count, &new_values[67]);
				insert_replace_slots_and_shift_remainders_and_runends_and_offsets(qf, 
																																					2, /* Insert to bucket */
																																					hash_bucket_index, 
																																					runstart_index, 
																																					p, 
																																					&new_values[67] - p, 
																																					0);
				/*modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);*/
			}
		}
		METADATA_WORD(qf, occupieds, hash_bucket_index) |= 1ULL << (hash_bucket_block_offset % 64);
		
		/*modify_metadata(qf, &qf->metadata->nelts, count);*/
	}

	if (lock) {
		qf_unlock(qf, hash_bucket_index, false);
	}

	return true;
}

static inline bool insert_advance(QF *qf, __uint128_t hash, uint64_t count, bool lock,
													bool spin, bool& isNew)
{
	uint64_t hash_remainder           = hash & BITMASK(qf->metadata->bits_per_slot);
	uint64_t hash_bucket_index        = hash >> qf->metadata->bits_per_slot;
	uint64_t hash_bucket_block_offset = hash_bucket_index % SLOTS_PER_BLOCK;
	/*uint64_t hash_bucket_lock_offset  = hash_bucket_index % NUM_SLOTS_TO_LOCK;*/

	if (lock) {
		if (!qf_lock(qf, hash_bucket_index, spin, false))
			return false;
	}

	uint64_t runend_index             = run_end(qf, hash_bucket_index);
	
	/* Empty slot */
	
	//if (might_be_empty(qf, hash_bucket_index) && runend_index ==
	//		hash_bucket_index) {
	//	isNew=true;
	//	METADATA_WORD(qf, runends, hash_bucket_index) |= 1ULL <<
	//		(hash_bucket_block_offset % 64);
	//	set_slot(qf, hash_bucket_index, hash_remainder);
	//	METADATA_WORD(qf, occupieds, hash_bucket_index) |= 1ULL <<
	//		(hash_bucket_block_offset % 64);
		
		/*modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);*/
		/*modify_metadata(qf, &qf->metadata->noccupied_slots, 1);*/
		/*modify_metadata(qf, &qf->metadata->nelts, 1);*/
		/* This trick will, I hope, keep the fast case fast. */
	//	if (count > 1) {
	//		return insert_advance(qf, hash, count - 1, false, false, isNew);
	//	}
	//} else { /* Non-empty slot */
		uint64_t new_values[67];
		int64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf,
																																	hash_bucket_index
																																	- 1) + 1;

		if (!is_occupied(qf, hash_bucket_index)) { /* Empty bucket, but its slot is occupied. */
			isNew=true;
			uint64_t *p = encode_counter(qf, hash_remainder, count, &new_values[67]);
			insert_replace_slots_and_shift_remainders_and_runends_and_offsets(qf, 
																																				0, 
																																				hash_bucket_index, 
																																				runstart_index, 
																																				p, 
																																				&new_values[67] - p, 
																																				0);
			/*modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);*/
		} else { /* Non-empty bucket */

			uint64_t current_remainder, current_count, current_end;

			/* Find the counter for this remainder, if one exists. */
			current_end = decode_counter(qf, runstart_index, &current_remainder,
																	 &current_count);
			while (current_remainder < hash_remainder && !is_runend(qf, current_end)) {
				runstart_index = current_end + 1;
				current_end = decode_counter(qf, runstart_index, &current_remainder,
																		 &current_count);	
			}

			/* If we reached the end of the run w/o finding a counter for this remainder,
				 then append a counter for this remainder to the run. */
			if (current_remainder < hash_remainder) {
				isNew=true;
				uint64_t *p = encode_counter(qf, hash_remainder, count, &new_values[67]);
				insert_replace_slots_and_shift_remainders_and_runends_and_offsets(qf, 
																																					1, /* Append to bucket */
																																					hash_bucket_index, 
																																					current_end + 1, 
																																					p, 
																																					&new_values[67] - p, 
																																					0);
				/*modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);*/
				/* Found a counter for this remainder.  Add in the new count. */
			} else if (current_remainder == hash_remainder) {
				isNew=false;
				uint64_t *p = encode_counter(qf, hash_remainder, current_count + count, &new_values[67]);
				insert_replace_slots_and_shift_remainders_and_runends_and_offsets(qf, 
																																					is_runend(qf, current_end) ? 1 : 2, 
																																					hash_bucket_index, 
																																					runstart_index, 
																																					p, 
																																					&new_values[67] - p, 
																																					current_end - runstart_index + 1);
				/* No counter for this remainder, but there are larger
					 remainders, so we're not appending to the bucket. */
			} else {
				isNew=true;
				uint64_t *p = encode_counter(qf, hash_remainder, count, &new_values[67]);
				insert_replace_slots_and_shift_remainders_and_runends_and_offsets(qf, 
																																					2, /* Insert to bucket */
																																					hash_bucket_index, 
																																					runstart_index, 
																																					p, 
																																					&new_values[67] - p, 
																																					0);
				/*modify_metadata(qf, &qf->metadata->ndistinct_elts, 1);*/
			}
		}
		METADATA_WORD(qf, occupieds, hash_bucket_index) |= 1ULL << (hash_bucket_block_offset % 64);
		
		/*modify_metadata(qf, &qf->metadata->nelts, count);*/
	//}

	if (lock) {
		qf_unlock(qf, hash_bucket_index, false);
	}

	return true;
}

inline static void _remove(QF *qf, __uint128_t hash, uint64_t count)
{
	uint64_t hash_remainder           = hash & BITMASK(qf->metadata->bits_per_slot);
	uint64_t hash_bucket_index        = hash >> qf->metadata->bits_per_slot;
	uint64_t current_remainder, current_count, current_end;
	uint64_t new_values[67];

	/* Empty bucket */
	if (!is_occupied(qf, hash_bucket_index))
		return;

	uint64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf, hash_bucket_index - 1) + 1;
	uint64_t original_runstart_index = runstart_index;
	int only_item_in_the_run = 0;

	/*Find the counter for this remainder, if one exists.*/
	current_end = decode_counter(qf, runstart_index, &current_remainder, &current_count);
	while (current_remainder < hash_remainder && !is_runend(qf, current_end)) {
		runstart_index = current_end + 1;
		current_end = decode_counter(qf, runstart_index, &current_remainder, &current_count);
	}
	/* remainder not found in the given run */
	if (current_remainder != hash_remainder)
		return;
	
	if (original_runstart_index == runstart_index && is_runend(qf, current_end))
		only_item_in_the_run = 1;

	/* endode the new counter */
	uint64_t *p = encode_counter(qf, hash_remainder,
															 count > current_count ? 0 : current_count - count,
															 &new_values[67]);
	remove_replace_slots_and_shift_remainders_and_runends_and_offsets(qf,
																																		only_item_in_the_run,
																																		hash_bucket_index,
																																		runstart_index,
																																		p,
																																		&new_values[67] - p,
																																		current_end - runstart_index + 1);

	// update the nelements.
	/*modify_metadata(qf, &qf->metadata->nelts, -count);*/
	/*qf->metadata->nelts -= count;*/
}

/***********************************************************************
 * Code that uses the above to implement key-value-counter operations. *
 ***********************************************************************/

void qf_init(QF *qf, uint64_t nslots, uint64_t key_bits, uint64_t value_bits,
						 bool mem, const char * path, uint32_t seed)
{
	uint64_t num_slots, xnslots, nblocks;
	uint64_t key_remainder_bits, bits_per_slot;
	uint64_t size;

	assert(popcnt(nslots) == 1); /* nslots must be a power of 2 */
	num_slots = nslots;
	assert(popcnt(nslots) == 1); /* nslots must be a power of 2 */
	xnslots = nslots + 10*sqrt((double)nslots);
	nblocks = (xnslots + SLOTS_PER_BLOCK - 1) / SLOTS_PER_BLOCK;
	key_remainder_bits = key_bits;
	while (nslots > 1) {
		assert(key_remainder_bits > 0);
		key_remainder_bits--;
		nslots >>= 1;
	}

	bits_per_slot = key_remainder_bits + value_bits;
	assert (BITS_PER_SLOT == 0 || BITS_PER_SLOT == qf->metadata->bits_per_slot);
	assert(bits_per_slot > 1);
#if BITS_PER_SLOT == 8 || BITS_PER_SLOT == 16 || BITS_PER_SLOT == 32 || BITS_PER_SLOT == 64
	size = nblocks * sizeof(qfblock);
#else
	size = nblocks * (sizeof(qfblock) + SLOTS_PER_BLOCK * bits_per_slot / 8);
#endif

	qf->mem = (qfmem *)calloc(sizeof(qfmem), 1);
	
	if (mem) {
		qf->metadata = (qfmetadata *)calloc(sizeof(qfmetadata), 1);

		qf->metadata->size = size;
		qf->metadata->seed = seed;
		qf->metadata->nslots = num_slots;
		qf->metadata->xnslots = qf->metadata->nslots +
			10*sqrt((double)qf->metadata->nslots);
		qf->metadata->key_bits = key_bits;
		qf->metadata->value_bits = value_bits;
		qf->metadata->key_remainder_bits = key_remainder_bits;
		qf->metadata->bits_per_slot = bits_per_slot;

		qf->metadata->range = qf->metadata->nslots;
		qf->metadata->range <<= qf->metadata->bits_per_slot;
		qf->metadata->nblocks = (qf->metadata->xnslots + SLOTS_PER_BLOCK - 1) /
			SLOTS_PER_BLOCK;
		qf->metadata->nelts = 0;
		qf->metadata->ndistinct_elts = 0;
		qf->metadata->noccupied_slots = 0;
		qf->metadata->num_locks = (qf->metadata->xnslots/NUM_SLOTS_TO_LOCK)+2;

		qf->blocks = (qfblock *)calloc(size, 1);

	} else {

		qf->mem->fd = open(path, O_RDWR | O_CREAT | O_TRUNC, S_IRWXU);
		if (qf->mem->fd < 0) {
			perror("Couldn't open file:\n");
			exit(EXIT_FAILURE);
		}
		
		/* prashantpandey: Commenting out fallocate call to preallocate space for
		 * the file on disk because fallocate is not supported on MAC OS. Revisit
		 * it later. */
		/*int ret;*/
		/*ret = fallocate(qf->mem->fd, 0, 0, size+sizeof(qfmetadata));*/
		/*if (ret < 0) {*/
			/*perror("Couldn't fallocate file:\n");*/
			/*exit(EXIT_FAILURE);*/
		/*}*/
		qf->metadata = (qfmetadata *)mmap(NULL, size+sizeof(qfmetadata), PROT_READ |
																			PROT_WRITE, MAP_SHARED, qf->mem->fd, 0);

		qf->metadata->seed = seed;
		qf->metadata->nslots = num_slots;
		qf->metadata->xnslots = qf->metadata->nslots +
														10*sqrt((double)qf->metadata->nslots);
		qf->metadata->key_bits = key_bits;
		qf->metadata->value_bits = value_bits;
		qf->metadata->key_remainder_bits = key_remainder_bits;
		qf->metadata->bits_per_slot = bits_per_slot;

		qf->metadata->range = qf->metadata->nslots;
		qf->metadata->range <<= qf->metadata->bits_per_slot;
		qf->metadata->nblocks = (qf->metadata->xnslots + SLOTS_PER_BLOCK - 1) /
			SLOTS_PER_BLOCK;
		qf->metadata->nelts = 0;
		qf->metadata->ndistinct_elts = 0;
		qf->metadata->noccupied_slots = 0;
		qf->metadata->num_locks = (qf->metadata->xnslots/NUM_SLOTS_TO_LOCK)+2;

		qf->blocks = (qfblock *)(qf->metadata + 1);
	}
	
	/* initialize all the locks to 0 */
	qf->mem->metadata_lock = 0;
	qf->mem->locks = (volatile int *)calloc(qf->metadata->num_locks,
																					sizeof(volatile int));
#ifdef LOG_WAIT_TIME
	qf->mem->wait_times = (wait_time_data* )calloc(qf->metadata->num_locks+1,
																						sizeof(wait_time_data));
#endif
}

/* The caller should call qf_init on the dest QF before calling this function. 
 */
void qf_copy(QF *dest, QF *src)
{
	memcpy(dest->mem, src->mem, sizeof(qfmem));
	memcpy(dest->metadata, src->metadata, sizeof(qfmetadata));
	memcpy(dest->blocks, src->blocks, src->metadata->size);
}

/* free up the memory if the QF is in memory.
 * else unmap the mapped memory from pagecache.
 *
 * It does not delete the file on disk for on-disk QF.
 */
void qf_destroy(QF *qf, bool mem)
{
	assert(qf->blocks != NULL);
	if (mem) {
		free(qf->mem);
		free(qf->metadata);
		free(qf->blocks);
	} else {
	munmap(qf->metadata, qf->metadata->size + sizeof(qfmetadata));
	close(qf->mem->fd);
	}
}

void qf_close(QF *qf)
{
	assert(qf->blocks != NULL);
	munmap(qf->metadata, qf->metadata->size + sizeof(qfmetadata));
	close(qf->mem->fd);
}

/* 
 * Will read the on-disk QF using mmap.
 * Data won't be copied in memory.
 *
 */
void qf_read(QF *qf, const char *path)
{
	struct stat sb;
	int ret;

	qf->mem = (qfmem *)calloc(sizeof(qfmem), 1);
	qf->mem->fd = open(path, O_RDWR, S_IRWXU);
	if (qf->mem->fd < 0) {
		perror("Couldn't open file:\n");
		exit(EXIT_FAILURE);
	}

	ret = fstat (qf->mem->fd, &sb);
	if ( ret < 0) {
		perror ("fstat");
		exit(EXIT_FAILURE);
	}

	if (!S_ISREG (sb.st_mode)) {
		fprintf (stderr, "%s is not a file\n", path);
		exit(EXIT_FAILURE);
	}

	qf->metadata = (qfmetadata *)mmap(NULL, sb.st_size, PROT_READ | PROT_WRITE, MAP_SHARED,
																qf->mem->fd, 0);

	qf->blocks = (qfblock *)(qf->metadata + 1);
}

void qf_reset(QF *qf)
{
	assert(popcnt(nslots) == 1); /* nslots must be a power of 2 */

	qf->metadata->nelts = 0;
	qf->metadata->ndistinct_elts = 0;
	qf->metadata->noccupied_slots = 0;

#ifdef LOG_WAIT_TIME
	memset(qf->mem->wait_times, 0, (qf->metadata->num_locks+1)*sizeof(wait_time_data));
#endif
#if BITS_PER_SLOT == 8 || BITS_PER_SLOT == 16 || BITS_PER_SLOT == 32 || BITS_PER_SLOT == 64
	memset(qf->blocks, 0, qf->metadata->nblocks* sizeof(qfblock));
#else
	memset(qf->blocks, 0, qf->metadata->nblocks*(sizeof(qfblock) + SLOTS_PER_BLOCK *
																		 qf->metadata->bits_per_slot / 8));
#endif
}

void qf_serialize(const QF *qf, const char *filename)
{
	FILE *fout;
	fout = fopen(filename, "wb+");
	if (fout == NULL) {
		perror("Error opening file for serializing\n");
		exit(EXIT_FAILURE);
	}

	fwrite(qf->metadata, sizeof(qfmetadata), 1, fout);

	/* we don't serialize the locks */
	fwrite(qf->blocks, qf->metadata->size, 1, fout);

	fclose(fout);
}

void qf_deserialize(QF *qf, const char *filename)
{
	FILE *fin;
	fin = fopen(filename, "rb");
	if (fin == NULL) {
		perror("Error opening file for deserializing\n");
		exit(EXIT_FAILURE);
	}

	qf->mem = (qfmem *)calloc(sizeof(qfmem), 1);
	qf->metadata = (qfmetadata *)calloc(sizeof(qfmetadata), 1);

	fread(qf->metadata, sizeof(qfmetadata), 1, fin);

	/* initlialize the locks in the QF */
	qf->metadata->num_locks = (qf->metadata->xnslots/NUM_SLOTS_TO_LOCK)+2;
	qf->mem->metadata_lock = 0;
	/* initialize all the locks to 0 */
	qf->mem->locks = (volatile int *)calloc(qf->metadata->num_locks, sizeof(volatile int));

	qf->blocks = (qfblock *)calloc(qf->metadata->size, 1);
	fread(qf->blocks, qf->metadata->size, 1, fin);

	fclose(fin);
}

bool qf_insert(QF *qf, uint64_t key, uint64_t value, uint64_t count, bool
							 lock, bool spin)
{
	/*uint64_t hash = (key << qf->metadata->value_bits) | (value & BITMASK(qf->metadata->value_bits));*/
	if (count == 1)
		return insert1(qf, key, lock, spin);
	else
		return insert(qf, key, count, lock, spin);
}

bool qf_insert_advance(QF *qf, uint64_t key, uint64_t value, uint64_t count, bool
							 lock, bool spin, bool& isNew)
{
	/*uint64_t hash = (key << qf->metadata->value_bits) | (value & BITMASK(qf->metadata->value_bits));*/
	if (count == 1)
		return insert1_advance(qf, key, lock, spin, isNew);
	else
		return insert_advance(qf, key, count, lock, spin, isNew);
}

uint64_t qf_count_key_value(const QF *qf, uint64_t key, uint64_t value)
{
	__uint128_t hash = key;
	uint64_t hash_remainder   = hash & BITMASK(qf->metadata->bits_per_slot);
	int64_t hash_bucket_index = hash >> qf->metadata->bits_per_slot;

	if (!is_occupied(qf, hash_bucket_index))
		return 0;

	int64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf,
																																hash_bucket_index-1)
		+ 1;
	if (runstart_index < hash_bucket_index)
		runstart_index = hash_bucket_index;

	/* printf("MC RUNSTART: %02lx RUNEND: %02lx\n", runstart_index, runend_index); */

	uint64_t current_remainder, current_count, current_end;
	do {
		current_end = decode_counter(qf, runstart_index, &current_remainder,
																 &current_count);
		if (current_remainder == hash_remainder)
			return current_count;
		runstart_index = current_end + 1;
	} while (!is_runend(qf, current_end));

	return 0;
}

/* initialize the iterator at the run corresponding
 * to the position index
 */
bool qf_iterator(QF *qf, QFi *qfi, uint64_t position)
{
	assert(position < qf->metadata->nslots);
	if (!is_occupied(qf, position)) {
		uint64_t block_index = position;
		uint64_t idx = bitselect(get_block(qf, block_index)->occupieds[0], 0);
		if (idx == 64) {
			while(idx == 64 && block_index < qf->metadata->nblocks) {
				block_index++;
				idx = bitselect(get_block(qf, block_index)->occupieds[0], 0);
			}
		}
		position = block_index * SLOTS_PER_BLOCK + idx;
	}

	qfi->qf = qf;
	qfi->num_clusters = 0;
	qfi->run = position;
	qfi->current = position == 0 ? 0 : run_end(qfi->qf, position-1) + 1;
	if (qfi->current < position)
		qfi->current = position;

#ifdef LOG_CLUSTER_LENGTH
	qfi->c_info = (cluster_data* )calloc(qf->metadata->nslots/32, sizeof(cluster_data));
	qfi->cur_start_index = position;
	qfi->cur_length = 1;
#endif

	if (qfi->current >= qf->metadata->nslots)
		return false;
	return true;
}

int qfi_get(QFi *qfi, uint64_t *key, uint64_t *value, uint64_t *count)
{
	if(qfi_end(qfi)){
		return 1;
	}
	assert(qfi->current < qfi->qf->metadata->nslots);

	uint64_t current_remainder, current_count;
	decode_counter(qfi->qf, qfi->current, &current_remainder, &current_count);

	*key = (qfi->run << qfi->qf->metadata->bits_per_slot) | current_remainder;
	*value = 0;   // for now we are not using value
	*count = current_count; 
	
	//qfi->qf->metadata->ndistinct_elts++;
	//qfi->qf->metadata->nelts += current_count;

	/*qfi->current = end_index;*/ 		//get should not change the current index
																		//of the iterator
	return 0;
}

int qfi_next(QFi *qfi)
{
	if (qfi_end(qfi))
		return 1;
	else {
		/* move to the end of the current counter*/
		uint64_t current_remainder, current_count;
		qfi->current = decode_counter(qfi->qf, qfi->current, &current_remainder,
																	&current_count);
		
		if (!is_runend(qfi->qf, qfi->current)) {
			qfi->current++;
#ifdef LOG_CLUSTER_LENGTH
			qfi->cur_length++;
#endif
			if (qfi->current > qfi->qf->metadata->nslots)
				return 1;
			return 0;
		}
		else {
#ifdef LOG_CLUSTER_LENGTH
			/* save to check if the new current is the new cluster. */
			uint64_t old_current = qfi->current;
#endif
			uint64_t block_index = qfi->run / SLOTS_PER_BLOCK;
			uint64_t rank = bitrank(get_block(qfi->qf, block_index)->occupieds[0],
															qfi->run % SLOTS_PER_BLOCK);
			uint64_t next_run = bitselect(get_block(qfi->qf,
																							block_index)->occupieds[0],
																		rank);
			if (next_run == 64) {
				rank = 0;
				while (next_run == 64 && block_index < qfi->qf->metadata->nblocks) {
					block_index++;
					next_run = bitselect(get_block(qfi->qf, block_index)->occupieds[0],
															 rank);
				}
			}
			if (block_index == qfi->qf->metadata->nblocks) {
				/* set the index values to max. */
				qfi->run = qfi->current = qfi->qf->metadata->xnslots;
				return 1;
			}
			qfi->run = block_index * SLOTS_PER_BLOCK + next_run;
			qfi->current++;
			if (qfi->current < qfi->run)
				qfi->current = qfi->run;
#ifdef LOG_CLUSTER_LENGTH
			if (qfi->current > old_current + 1) { /* new cluster. */
				if (qfi->cur_length > 10) {
					qfi->c_info[qfi->num_clusters].start_index = qfi->cur_start_index;
					qfi->c_info[qfi->num_clusters].length = qfi->cur_length;
					qfi->num_clusters++;
				}
				qfi->cur_start_index = qfi->run;
				qfi->cur_length = 1;
			} else {
				qfi->cur_length++;
			}
#endif
			return 0;
		}
	}
}

//inline 
int qfi_end(QFi *qfi)
{
	if (qfi->current >= qfi->qf->metadata->xnslots /*&& is_runend(qfi->qf, qfi->current)*/)
		return 1;
	else
		return 0;
}

/*
 * Merge qfa and qfb into qfc 
 */
/*
 * iterate over both qf (qfa and qfb)
 * simultaneously 
 * for each index i 
 * min(get_value(qfa, ia) < get_value(qfb, ib))  
 * insert(min, ic) 
 * increment either ia or ib, whichever is minimum.  
 */
void qf_merge(QF *qfa, QF *qfb, QF *qfc) 
{
	QFi qfia, qfib;
	qf_iterator(qfa, &qfia, 0);
	qf_iterator(qfb, &qfib, 0);

	uint64_t keya, valuea, counta, keyb, valueb, countb;
	qfi_get(&qfia, &keya, &valuea, &counta);
	qfi_get(&qfib, &keyb, &valueb, &countb);
	do {
		if (keya < keyb) {
			qf_insert(qfc, keya, valuea, counta, true, true);
			qfi_next(&qfia);
			qfi_get(&qfia, &keya, &valuea, &counta);
		}
		else {
			qf_insert(qfc, keyb, valueb, countb, true, true);
			qfi_next(&qfib);
			qfi_get(&qfib, &keyb, &valueb, &countb);
		}
	} while(!qfi_end(&qfia) && !qfi_end(&qfib));

	if (!qfi_end(&qfia)) {
		do {
			qfi_get(&qfia, &keya, &valuea, &counta);
			qf_insert(qfc, keya, valuea, counta, true, true);
		} while(!qfi_next(&qfia));
	}
	if (!qfi_end(&qfib)) {
		do {
			qfi_get(&qfib, &keyb, &valueb, &countb);
			qf_insert(qfc, keyb, valueb, countb, true, true);
		} while(!qfi_next(&qfib));
	}

	return;
}

/*
 * Merge an array of qfs into the resultant QF
 */
void qf_multi_merge(QF *qf_arr[], int nqf, QF *qfr)
{
	int i;
	QFi qfi_arr[nqf];
	int flag = 0;
	int smallest_i = 0;
	uint64_t smallest_key = UINT64_MAX;
	for (i=0; i<nqf; i++) {
		qf_iterator(qf_arr[i], &qfi_arr[i], 0);
	}

	while (!flag) {
		uint64_t keys[nqf];
		uint64_t values[nqf];
		uint64_t counts[nqf];
		for (i=0; i<nqf; i++)
			qfi_get(&qfi_arr[i], &keys[i], &values[i], &counts[i]);
		
		do {
			smallest_key = UINT64_MAX;
			for (i=0; i<nqf; i++) {
				if (keys[i] < smallest_key) {
					smallest_key = keys[i]; smallest_i = i;
				}
			}
			qf_insert(qfr, keys[smallest_i], values[smallest_i], counts[smallest_i],
								true, true);
			qfi_next(&qfi_arr[smallest_i]);
			qfi_get(&qfi_arr[smallest_i], &keys[smallest_i], &values[smallest_i],
							&counts[smallest_i]);
		} while(!qfi_end(&qfi_arr[smallest_i]));

		/* remove the qf that is exhausted from the array */
		if (smallest_i < nqf-1)
			memmove(&qfi_arr[smallest_i], &qfi_arr[smallest_i+1],
							(nqf-smallest_i-1)*sizeof(qfi_arr[0]));
		nqf--;
		if (nqf == 1)
			flag = 1;
	}
	if (!qfi_end(&qfi_arr[0])) {
		do {
			uint64_t key, value, count;
			qfi_get(&qfi_arr[0], &key, &value, &count);
			qf_insert(qfr, key, value, count, true, true);
		} while(!qfi_next(&qfi_arr[0]));
	}

	return;
}

/* find cosine similarity between two QFs. */
uint64_t qf_inner_product(QF *qfa, QF *qfb)
{
	uint64_t acc = 0;
	QFi qfi;
	QF *qf_mem, *qf_disk;

	// create the iterator on the larger QF.
	if (qfa->metadata->size > qfb->metadata->size) {
		qf_mem = qfb;
		qf_disk = qfa;
	} else {
		qf_mem = qfa;
		qf_disk = qfb;
	}

	qf_iterator(qf_disk, &qfi, 0);
	do {
		uint64_t key = 0, value = 0, count = 0;
		uint64_t count_mem;
		qfi_get(&qfi, &key, &value, &count);
		if ((count_mem = qf_count_key_value(qf_mem, key, 0)) > 0) {
			acc += count*count_mem;
		}
	} while (!qfi_next(&qfi));

	return acc;
}

/* find cosine similarity between two QFs. */
void qf_intersect(QF *qfa, QF *qfb, QF *qfr)
{
	QFi qfi;
	QF *qf_mem, *qf_disk;

	// create the iterator on the larger QF.
	if (qfa->metadata->size > qfb->metadata->size) {
		qf_mem = qfb;
		qf_disk = qfa;
	} else {
		qf_mem = qfa;
		qf_disk = qfb;
	}

	qf_iterator(qf_disk, &qfi, 0);
	do {
		uint64_t key = 0, value = 0, count = 0;
		qfi_get(&qfi, &key, &value, &count);
		if (qf_count_key_value(qf_mem, key, 0) > 0)
			qf_insert(qfr, key, value, count, false, false);
	} while (!qfi_next(&qfi));
}

/* magnitude of a QF. */
uint64_t qf_magnitude(QF *qf)
{
	return sqrt(qf_inner_product(qf, qf));
}

void qf_print_metadata(QF *qf){}
/*
	cout<<"#metadata"<<endl
		<<"size: "<<qf->metadata->size<<endl
		<<"seed: "<<qf->metadata->seed<<endl
		<<"nslots: "<<qf->metadata->nslots<<endl
		<<"xnslots: "<<qf->metadata->xnslots<<endl
		<<"key_bits: "<<qf->metadata->key_bits<<endl
		<<"value_bits: "<<qf->metadata->value_bits<<endl
		<<"key_remainder_bits: "<<qf->metadata->key_remainder_bits<<endl
		<<"bits_per_slots: "<<qf->metadata->bits_per_slot<<endl
		//<<"range: "<<qf->metadata->range<<endl
		<<"nelts: "<<qf->metadata->nelts<<endl
		<<"ndistinct_elts: "<<qf->metadata->ndistinct_elts<<endl
		<<"noccupied_slots: "<<qf->metadata->noccupied_slots<<endl
		<<"num_locks: "<<qf->metadata->num_locks<<endl;
}
*/

//CHRISTINA
uint64_t firstBlock_offsetNotMax(const QF *qf, uint64_t start_block){
	while(get_block(qf, start_block)->offset == 255){
		start_block++;
	}
	return start_block;
}

void qf_clean_singleton(const QF* qf, uint64_t start_bucket_id, uint64_t end_bucket_id, uint64_t* removed_elts){
	assert(start_bucket_id < end_bucket_id);
	uint64_t least_significant_bit_mask = 1ULL << (qf->metadata->bits_per_slot-1);
	//lock the beginging and end block
	uint64_t last_empty_slot, bucket_idx, run_start, run_end, insert_idx, block_idx, offset;
	last_empty_slot = bucket_idx = run_start = insert_idx = start_bucket_id;
	uint64_t remainder, tmp;
	while(bucket_idx <= end_bucket_id){
		//while(insert_idx < bucket_idx){
			//not efficient enough
		//	set_slot(qf, insert_idx, 0);
		//	insert_idx++;
		//}
		last_empty_slot = insert_idx;
		run_end = run_start;
		while(true){
			if(is_runend(qf, run_end)){
				(*removed_elts)++;
				break;
			}
			remainder = get_slot(qf, run_end);
			tmp = get_slot(qf, ++run_end);
			if(remainder >= tmp){
				if(insert_idx+1== run_end){
					//skip the counter part
					if(tmp==0){
						run_end ++;
					}
					tmp = get_slot(qf, run_end);
					while(tmp & least_significant_bit_mask){
						tmp = get_slot(qf, ++run_end);
					}
					insert_idx = run_end+1;
				}else{
					set_slot(qf, insert_idx++, remainder);
					if(tmp==0){
						set_slot(qf, insert_idx++, 0);
						run_end++;
					}
					tmp = get_slot(qf, run_end);
					while(tmp & least_significant_bit_mask){
						set_slot(qf, insert_idx++, tmp);
						tmp=get_slot(qf, ++run_end);
					}
					set_slot(qf, insert_idx++, tmp);
				}
				if(is_runend(qf, run_end)){
					//(*removed_elts)++;
					break;
				}
				run_end ++;
			}else{
				(*removed_elts)++;
			}
		}
		if(last_empty_slot == insert_idx){
			METADATA_WORD(qf, occupieds, bucket_idx) &= ~(1ULL<<((bucket_idx%SLOTS_PER_BLOCK)%64));	
			METADATA_WORD(qf, runends, run_end) &= ~(1ULL<<((run_end%SLOTS_PER_BLOCK)%64));		
		}else if(run_end+1 != insert_idx){
			METADATA_WORD(qf, runends, run_end) &= ~(1ULL<<((run_end%SLOTS_PER_BLOCK)%64));		
			METADATA_WORD(qf, runends, insert_idx-1) |= (1ULL<<(((insert_idx-1)%SLOTS_PER_BLOCK)%64));	
		}
		run_start = run_end +1;
		bucket_idx++;
		while(!is_occupied(qf, bucket_idx) && bucket_idx <= end_bucket_id){
			bucket_idx++;
		}
		while(insert_idx < bucket_idx){
			//not efficient enough
			set_slot(qf, insert_idx, 0);
			insert_idx++;
		}
	}


	//update offset
	for(block_idx = start_bucket_id/SLOTS_PER_BLOCK+1; block_idx <= end_bucket_id/SLOTS_PER_BLOCK; block_idx++){
		offset = block_offset_strict(qf, block_idx);
		if(offset > BLOCK_OFFSET_MAX){
			get_block(qf, block_idx)->offset = BLOCK_OFFSET_MAX;
		}else{
			get_block(qf, block_idx)->offset = offset;
		}
	}
}

void qf_clean_singleton_discrete(const QF* qf, uint64_t start_bucket_id, uint64_t end_bucket_id, uint64_t* removed_elts){
	uint64_t start, end;
	start = start_bucket_id;//find_first_nonempty_slot(qf, start_bucket_id);
	while(start < end_bucket_id){
		end = find_first_empty_slot(qf, start+1)-1;
		qf_clean_singleton(qf, start, end, removed_elts);
		start = find_first_nonempty_slot(qf, end+1);
	}
}

void qf_clean_singleton_with_lock_atStart(const QF* qf, uint64_t start_bucket_id, uint64_t end_bucket_id, uint64_t* removed_elts){
	uint64_t least_significant_bit_mask = 1ULL << (qf->metadata->bits_per_slot-1);
	assert(start_bucket_id < end_bucket_id);
	//lock the beginging and end block
	uint64_t last_empty_slot, bucket_idx, run_start, run_end, insert_idx, block_idx, offset, start_block_idx, end_block_idx, count;
	last_empty_slot = bucket_idx = run_start = insert_idx = start_bucket_id;
	block_idx = bucket_idx/SLOTS_PER_BLOCK;

	start_block_idx = start_bucket_id/SLOTS_PER_BLOCK;
	end_block_idx = end_bucket_id/SLOTS_PER_BLOCK;
	
	uint64_t remainder, tmp;	
	//lock for the first block;
	qf_spin_lock(&qf->mem->locks[start_bucket_id/NUM_SLOTS_TO_LOCK], true);
	//process the first block;
	while(bucket_idx <= end_bucket_id && block_idx == start_block_idx){
		last_empty_slot = insert_idx;
		run_end = run_start;
		while(true){
			if(is_runend(qf, run_end)){
				(*removed_elts)++;
				break;
			}
			remainder = get_slot(qf, run_end);
			tmp = get_slot(qf, ++run_end);
			if(remainder >= tmp){
				if(insert_idx+1== run_end){
					//skip the counter part
					if(tmp==0){
						run_end ++;
					}
					tmp = get_slot(qf, run_end);
					while(tmp & least_significant_bit_mask){
						tmp = get_slot(qf, ++run_end);
					}
					insert_idx = run_end+1;
				}else{
					set_slot(qf, insert_idx++, remainder);
					if(tmp==0){
						set_slot(qf, insert_idx++, 0);
						run_end++;
					}
					tmp = get_slot(qf, run_end);
					while(tmp & least_significant_bit_mask){
						set_slot(qf, insert_idx++, tmp);
						tmp=get_slot(qf, ++run_end);
					}
					set_slot(qf, insert_idx++, tmp);
				}
				if(is_runend(qf, run_end)){
					//(*removed_elts)++;
					break;
				}
				run_end ++;
			}else{
				(*removed_elts)++;
			}
		}
		if(last_empty_slot == insert_idx){
			METADATA_WORD(qf, occupieds, bucket_idx) &= ~(1ULL<<((bucket_idx%SLOTS_PER_BLOCK)%64));	
			METADATA_WORD(qf, runends, run_end) &= ~(1ULL<<((run_end%SLOTS_PER_BLOCK)%64));		
		}else if(run_end+1 != insert_idx){
			METADATA_WORD(qf, runends, run_end) &= ~(1ULL<<((run_end%SLOTS_PER_BLOCK)%64));		
			METADATA_WORD(qf, runends, insert_idx-1) |= (1ULL<<(((insert_idx-1)%SLOTS_PER_BLOCK)%64));	
		}
		run_start = run_end +1;
		bucket_idx++;
		while(!is_occupied(qf, bucket_idx) && bucket_idx <= end_bucket_id){
			bucket_idx++;
		}		
		while(insert_idx < bucket_idx){
			//not efficient enough
			set_slot(qf, insert_idx, 0);
			insert_idx++;
		}
		block_idx = bucket_idx/SLOTS_PER_BLOCK;
	}
	qf_spin_unlock(&qf->mem->locks[start_bucket_id/NUM_SLOTS_TO_LOCK]);
	
	while(bucket_idx <= end_bucket_id){
		//block_idx = bucket_idx/SLOTS_PER_BLOCK;
		last_empty_slot = insert_idx;
		run_end = run_start;
		while(true){
			if(is_runend(qf, run_end)){
				(*removed_elts)++;
				break;
			}
			remainder = get_slot(qf, run_end);
			tmp = get_slot(qf, ++run_end);
			if(remainder >= tmp){
				if(insert_idx+1 == run_end){
					//skip the counter part
					if(tmp==0){
						run_end ++;
					}
					tmp = get_slot(qf, run_end);
					while(tmp & least_significant_bit_mask){
						tmp = get_slot(qf, ++run_end);
					}
					insert_idx = run_end+1;
				}else{
					set_slot(qf, insert_idx++, remainder);
					if(tmp==0){
						set_slot(qf, insert_idx++, 0);
						run_end++;
					}
					tmp = get_slot(qf, run_end);
					while(tmp & least_significant_bit_mask){
						set_slot(qf, insert_idx++, tmp);
						tmp=get_slot(qf, ++run_end);
					}
					set_slot(qf, insert_idx++, tmp);
				}
				if(is_runend(qf, run_end)){
					//(*removed_elts)++;
					break;
				}
				run_end ++;
			}else{
				(*removed_elts)++;
			}
		}
		if(last_empty_slot == insert_idx){
			METADATA_WORD(qf, occupieds, bucket_idx) &= ~(1ULL<<((bucket_idx%SLOTS_PER_BLOCK)%64));	
			METADATA_WORD(qf, runends, run_end) &= ~(1ULL<<((run_end%SLOTS_PER_BLOCK)%64));		
		}else if(run_end+1 != insert_idx){
			METADATA_WORD(qf, runends, run_end) &= ~(1ULL<<((run_end%SLOTS_PER_BLOCK)%64));		
			METADATA_WORD(qf, runends, insert_idx-1) |= (1ULL<<(((insert_idx-1)%SLOTS_PER_BLOCK)%64));	
		}
		run_start = run_end +1;
		bucket_idx++;
		while(!is_occupied(qf, bucket_idx) && bucket_idx <= end_bucket_id){
			bucket_idx++;
		}
		while(insert_idx < bucket_idx){
			//not efficient enough
			set_slot(qf, insert_idx, 0);
			insert_idx++;
		}
	}
	//update offset
	qf_spin_lock(&qf->mem->locks[start_bucket_id/NUM_SLOTS_TO_LOCK], true);
	for(block_idx = start_bucket_id/SLOTS_PER_BLOCK+1; block_idx <= end_bucket_id/SLOTS_PER_BLOCK; block_idx++){
		offset = block_offset_strict(qf, block_idx);
		if(offset > BLOCK_OFFSET_MAX){
			get_block(qf, block_idx)->offset = BLOCK_OFFSET_MAX;
		}else{
			get_block(qf, block_idx)->offset = offset;
		}
	}
	qf_spin_unlock(&qf->mem->locks[start_bucket_id/NUM_SLOTS_TO_LOCK]);
}

uint64_t popcnt_runends(const QF* qf){
	uint64_t counts = 0;
	for(int x = 0; x < qf->metadata->nblocks; x++){
		counts += popcnt(get_block(qf, x)->runends[0]);
	}
	return counts;
}
uint64_t popcnt_occupieds(const QF* qf){
	uint64_t counts = 0;
	for(int x = 0; x <= qf->metadata->nblocks; x++){
		counts += popcnt(get_block(qf, x)->occupieds[0]);
	}
	return counts;
}
bool check_offset(const QF* qf){
	uint64_t offset, offset_real;
	for(int x = 1; x < qf->metadata->nblocks; x++){
		 offset_real = block_offset_strict(qf, x);
		 offset = get_block(qf, x)->offset;
		 if(offset_real <= BLOCK_OFFSET_MAX){
				if(offset != offset_real){
					printf("Block %d: offset.%lu offset_real.%lu\n", x, offset, offset_real);
					return false;
				}
		 }
	}
	return true;
}
#ifdef GRAPH_TRAVERSE
bool qf_is_traveled(const QF *qf, uint64_t index){
	return (METADATA_WORD(qf, traveled, index) >> ((index % SLOTS_PER_BLOCK) % 64)) & 1ULL;
}

void qf_set_traveled(const QF* qf, uint64_t index){
	uint64_t block_idx = index/SLOTS_PER_BLOCK;
	uint64_t block_offset = (index % SLOTS_PER_BLOCK) % 64;
	get_block(qf, block_idx)->traveled[block_offset/64] |= (1ULL<<(block_offset%64));
	//qf->blocks[block_idx].traveled[0] |= (1ULL<<block_offset);
}

int qfi_next_untraveled(QFi * qfi){
	int isEnd = qfi_next(qfi);
	while(!isEnd && qf_is_traveled(qfi->qf, qfi->current)){
		isEnd = qfi_next(qfi);
	}
	return isEnd;
}

//return whether is traveled (true) or not traveled (false)
//and then set it to be traveled
bool qf_count_key_value_set_traveled(const QF *qf, uint64_t key, uint64_t value, uint64_t* count){
	__uint128_t hash = key;
	uint64_t hash_remainder   = hash & BITMASK(qf->metadata->bits_per_slot);
	int64_t hash_bucket_index = hash >> qf->metadata->bits_per_slot;

	if (!is_occupied(qf, hash_bucket_index)){
		*count = 0;
		return false;
	}

	int64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf,
																																hash_bucket_index-1)
		+ 1;
	if (runstart_index < hash_bucket_index)
		runstart_index = hash_bucket_index;

	/* printf("MC RUNSTART: %02lx RUNEND: %02lx\n", runstart_index, runend_index); */

	uint64_t current_remainder, current_count, current_end;
	do {
		current_end = decode_counter(qf, runstart_index, &current_remainder,
																 &current_count);
		if (current_remainder == hash_remainder){
			*count = current_count;
			if(qf_is_traveled(qf, runstart_index)){
				return true;
			}else{
				qf_set_traveled(qf, runstart_index);
				return false;
			}
		}
		runstart_index = current_end + 1;
	} while (!is_runend(qf, current_end));

	*count = 0;
	return false;
}

//return whether is traveled
//and set it to be traveled
bool qf_count_key_value_is_traveled(const QF *qf, uint64_t key, uint64_t value, uint64_t* count){
	__uint128_t hash = key;
	uint64_t hash_remainder   = hash & BITMASK(qf->metadata->bits_per_slot);
	int64_t hash_bucket_index = hash >> qf->metadata->bits_per_slot;

	if (!is_occupied(qf, hash_bucket_index)){
		*count = 0;
		return false;
	}

	int64_t runstart_index = hash_bucket_index == 0 ? 0 : run_end(qf,
																																hash_bucket_index-1)
		+ 1;
	if (runstart_index < hash_bucket_index)
		runstart_index = hash_bucket_index;

	/* printf("MC RUNSTART: %02lx RUNEND: %02lx\n", runstart_index, runend_index); */

	uint64_t current_remainder, current_count, current_end;
	do {
		current_end = decode_counter(qf, runstart_index, &current_remainder,
																 &current_count);
		if (current_remainder == hash_remainder){
			*count = current_count;
			return qf_is_traveled(qf, runstart_index);
		}
		runstart_index = current_end + 1;
	} while (!is_runend(qf, current_end));

	*count = 0;
	return false;
}

#endif

//CHRSTINA END
