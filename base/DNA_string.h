#pragma once

#ifndef DNA_STRING_H
#define DNA_STRING_H

#include "global.h"
#include "nthash.hpp"

namespace DNA{
	const char bases[4]={'A', 'C', 'G', 'T'};
	//Map bases to 2 bit representation: A|a:00; C|c:01; G|g:10; T|t:11; N|n:arbitrary.
	//Every 8 bits store 4 bases.
	const uint8_t base_mask[4]={0xC0, 0x30, 0x0C, 0x03};
	const uint8_t base2bits[]={0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,//0x0F
													0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,//0x1F
													0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,//0x2F
													0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,//0x3F
													0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,//0x4F
													0x00, 0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,//0x5F
													0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,//0x6F
													0x00, 0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,//0x7F
													0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,//0x8F
													0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,//0x9F
													0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,//0xAF
													0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,//0xBF
													0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,//0xCF
													0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,//0xDF
													0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,//0xEF
													0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};//0xFF
	const char bits2base[]={'A', 'C', 'G', 'T', 'C', ' ', ' ', ' ', 'G', ' ', ' ', ' ', 'T', ' ', ' ', ' ',//0x0F
													'C', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x1F
													'G', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x2F
													'T', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x3F
													'C', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x4F
													' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x5F
													' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x6F
													' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x7F
													'G', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x8F
													' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x9F
													' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0xAF
													' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0xBF
													'T', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0xCF
													' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0xDF
													' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0xEF
													' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '};//0xFF
	const char bits2RCbase[]={'T', 'G', 'C', 'A', 'G', ' ', ' ', ' ', 'C', ' ', ' ', ' ', 'A', ' ', ' ', ' ',//0x0F
													'G', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x1F
													'C', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x2F
													'A', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x3F
													'G', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x4F
													' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x5F
													' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x6F
													' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x7F
													'C', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x8F
													' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x9F
													' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0xAF
													' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0xBF
													'A', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0xCF
													' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0xDF
													' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0xEF
													' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '};//0xFF											
	//const char bits2base[]={'A','C','G','T'};
}

class DNAString{
public:
	DNAString();

	DNAString(vector<uint8_t> d, size_t dna_base_num);

	DNAString(const DNAString& dnaStr);

	DNAString& operator= (const DNAString& dnaStr);

	DNAString& operator= (const string& str);

	const uint8_t* data() const;

	char operator[](size_t n);

	DNAString& clear();

	DNAString(const string& str);

	DNAString make_DNAString(const string& str);

	DNAString& append(const char& base);

	//remove the left most base
	DNAString& pop();

	DNAString& operator +=(const char& base);

	size_t dna_base_num() const;

	size_t length() const;

	size_t data_size() const;

	DNAString& RC();

	DNAString get_RC();

	//pos is 0-based coordinate
	void replace(const size_t pos, const char c);

	//start_pos is 0-based coordinate
	DNAString substr(size_t start_pos, size_t len=string::npos) const;

	string to_str() const;

	bool is_palindrome() const;

	bool is_hairpin(size_t len=0) const;

	//return true if kmer satisfy (A|T|G|C)*
	bool is_simple() const;

	friend bool operator== (const DNAString& lhs, const DNAString& rhs);
	friend bool operator!= (const DNAString& lhs, const DNAString& rhs);
	friend bool operator< (const DNAString& lhs, const DNAString& rhs);
	friend ostream& operator<<( ostream &output, const DNAString& dnaStr);

private:
	vector<uint8_t> _data;
	size_t _dna_base_num;
	int8_t _offset;//Number of bits shift to the left in order to store the next base. 
};


DNAString operator + (const DNAString& str, const char& base);
DNAString operator + (const char& base, const DNAString& str);

//canonical nthash
uint64_t NTPC64(const DNAString& dnaStr, const unsigned k, uint64_t& fhVal, uint64_t& rhVal);

#endif
