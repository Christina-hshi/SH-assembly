/*
 * ============================================================================
 *  
 *        Authors: Christina SHI <hshi@cse.cuhk.edu.hk>
 *                 Kevin Yip <kevinyip@cse.cuhk.edu.hk>
 *
 * ============================================================================
 */

#include "DNA_string.h"

DNAString::DNAString():_dna_base_num(0), _offset(0){}

DNAString::DNAString(vector<uint8_t> d, size_t dna_base_num):_data(d.begin(), d.end()), _dna_base_num(dna_base_num){
	this->_offset = (d.size()*4 - dna_base_num - 1)*2;
}

DNAString::DNAString(const DNAString& dnaStr){
	this->_data.assign(dnaStr._data.begin(), dnaStr._data.end());
	this->_dna_base_num = dnaStr._dna_base_num;
	this->_offset = dnaStr._offset;
}

DNAString& DNAString::operator= (const DNAString& dnaStr){
	this->_data.assign(dnaStr._data.begin(), dnaStr._data.end());
	this->_dna_base_num = dnaStr._dna_base_num;
	this->_offset = dnaStr._offset;
	return *this;
}

DNAString& DNAString::operator= (const string& str){
	this->clear();
	this->_dna_base_num = str.size();
	size_t anchor, buf_anchor, buf_size;
	if(_dna_base_num % 4){
		buf_size = str.size()/4+1;
	}else{
		buf_size = str.size()/4;
	} 
	this->_offset = (buf_size*4-this->_dna_base_num-1)*2;
	_data.resize(buf_size, 0x00);
	anchor=buf_anchor=0;
	while(anchor < str.size()){
		uint8_t c=0x00;
		for(int idx = 0; idx < 4 && anchor < str.size(); idx++, anchor++){
			//cout<<str[anchor]<<":"<<DNA::base2bits[anchor]<<(6-idx*2)<<endl;
			c |= (DNA::base2bits[str[anchor]]<<(6-idx*2));
		}
		_data[buf_anchor++] =c;
	}
	return *this;
}

const uint8_t* DNAString::data() const{
	return this->_data.data();
}

char DNAString::operator[](size_t n){
	assert(n < this->_dna_base_num);
	return DNA::bits2base[this->_data[n/4] & DNA::base_mask[n%4]];
}

DNAString& DNAString::clear(){
	this->_data.clear();
	this->_dna_base_num = 0;
	this->_offset = 0;
	return *this;
}

DNAString::DNAString(const string& str):_dna_base_num(str.size()){
	size_t anchor, buf_anchor, buf_size;
	if(_dna_base_num % 4){
		buf_size = str.size()/4+1;
	}else{
		buf_size = str.size()/4;
	} 
	this->_offset = (buf_size*4-this->_dna_base_num-1)*2;
	_data.resize(buf_size, 0x00);
	anchor=buf_anchor=0;
	while(anchor < str.size()){
		uint8_t c=0x00;
		for(int idx = 0; idx < 4 && anchor < str.size(); idx++, anchor++){
			//cout<<str[anchor]<<":"<<DNA::base2bits[anchor]<<(6-idx*2)<<endl;
			c |= (DNA::base2bits[str[anchor]]<<(6-idx*2));
		}
		_data[buf_anchor++] =c;
	}
}

DNAString DNAString::make_DNAString(const string& str){
	DNAString tmp;
	tmp._dna_base_num = str.size();
	
	size_t anchor, buf_anchor, buf_size;
	if(tmp._dna_base_num % 4){
		buf_size = str.size()/4+1;
	}else{
		buf_size = str.size()/4;
	} 
	tmp._offset = (buf_size*4 - tmp._dna_base_num-1)*2;
	tmp._data.resize(buf_size, 0x00);
	anchor=buf_anchor=0;
	while(anchor < str.size()){
		uint8_t c=0x00;
		for(int idx = 0; idx < 4 && anchor < str.size(); idx++, anchor++){
			//cout<<str[anchor]<<":"<<DNA::base2bits[anchor]<<(6-idx*2)<<endl;
			c |= (DNA::base2bits[str[anchor]]<<(6-idx*2));
		}
		tmp._data[buf_anchor++] =c;
	}
	return tmp;
}

DNAString& DNAString::append(const char& base){
	if(this->_offset == -2){
		_data.push_back(DNA::base2bits[base]<<6);
		this->_offset = 4;
	}else{
		_data.back() = (_data.back() | (DNA::base2bits[base]<<this->_offset));
		this->_offset -= 2;
	}
	_dna_base_num++;
	return *this;
}

//remove the left most base
DNAString& DNAString::pop(){
	int last_byte = this->_dna_base_num/4;
	if(this->_dna_base_num % 4 == 0){
		last_byte--;
	}
	for(int x = 0; x< last_byte; x++){
		this->_data[x] <<= 2;
		this->_data[x] |= this->_data[x+1]>>6;
	}
	this->_data[last_byte] <<= 2;
	if(this->_offset == 4){
		this->_offset = -2;
		this->_data.resize(this->_data.size()-1);
	}else{
		this->_offset += 2;
	}
	this->_dna_base_num --;
	return *this;
}

DNAString& DNAString::operator +=(const char& base){
	this->append(base);
	return *this;
}

size_t DNAString::dna_base_num() const{
	return _dna_base_num;
}

size_t DNAString::length() const{
	return _dna_base_num;
}

size_t DNAString::data_size() const{
	return this->_data.size();
}

DNAString& DNAString::RC(){
	vector<uint8_t> _rc_data(_data.size(), 0x00);
	size_t _rc_slot_idx; _rc_slot_idx = 0;
	size_t base_2_process = this->_dna_base_num;
	size_t backward_slot_idx = (base_2_process-1)/4;
	int backward_offset = (base_2_process % 4 + 3)%4;

	while(base_2_process > 0){
		uint8_t tmp=0x00;
		for(int x = 0; x< 4 && base_2_process > 0; x++, base_2_process--){
			if(backward_offset == -1){
				backward_slot_idx--;
				backward_offset = 3;
			}
			uint8_t extract = ((~(this->_data[backward_slot_idx])) & DNA::base_mask[backward_offset]);
			
			if(backward_offset == x){
				tmp |= extract;
			}else if(backward_offset < x){
				tmp |= (extract >> ((x - backward_offset)*2));
			}else{
				tmp |= (extract << ((backward_offset - x)*2));
			}
			backward_offset--;
		}
		_rc_data[_rc_slot_idx++] = tmp;
	}
	this->_data.assign(_rc_data.begin(), _rc_data.end());

	return *this;
}

DNAString DNAString::get_RC(){
	DNAString tmp(*this);
	tmp.RC();
	return tmp;
}

//pos is 0-based coordinate
void DNAString::replace(const size_t pos, const char c){
	assert(pos < this->_dna_base_num);
	uint8_t d = this->_data[pos/4];
	d &= ~(DNA::base_mask[pos%4]);
	d |= (DNA::base2bits[c]<<((3 - pos%4)*2));
	this->_data[pos/4] = d;
}

//start_pos is 0-based coordinate
DNAString DNAString::substr(size_t start_pos, size_t len) const{
	assert(start_pos < this->_dna_base_num);

	if(len == string::npos){
		len = this->_dna_base_num - start_pos;
	}else if(start_pos + len > this->_dna_base_num){
		len = this->_dna_base_num - start_pos;
	}

	size_t byte_num = len/4;
	if(len%4){
		byte_num++;
	} 

	vector<uint8_t> data(byte_num, 0x00);
	size_t data_slot_idx = 0;
	size_t _data_slot_idx = start_pos/4;
	size_t _data_slot_offset = start_pos%4;

	if(start_pos==0){
		while(data_slot_idx <= (len-1)/4 ){
			data[data_slot_idx] = this->_data[data_slot_idx];
			data_slot_idx++;
		}
	}else{
		while(_data_slot_idx <= (start_pos+len-1)/4){
			uint8_t slot = (this->_data[_data_slot_idx]<<(_data_slot_offset*2));
			if((_data_slot_idx+1) <= (start_pos+len-1)/4 && _data_slot_idx < this->_data.size()){
				slot |= (this->_data[_data_slot_idx+1] >> (8 - _data_slot_offset*2));
			}
			data[data_slot_idx++] = slot;
			_data_slot_idx++;
		}
	}
	uint8_t mask = 0x00;
	for(int x = 0; x<len%4; x++){
		mask |= DNA::base_mask[x];
	}
	data[data_slot_idx-1] &= mask;//reset the rest of the bits

	return DNAString(data, len);
}

string DNAString::to_str() const{
	string str(this->_dna_base_num, ' ');
	size_t base_idx=0;
	for(size_t slot=0; slot < this->_data.size(); slot++){
		for(size_t x=0; x < 4 && base_idx < _dna_base_num ; x++, base_idx++){
			str[base_idx]=DNA::bits2base[(this->_data[slot]) & DNA::base_mask[x]];
		}
	}
	return str;
}

bool DNAString::is_palindrome() const{
	if(this->_dna_base_num % 2){
		return false;
	}
	return is_hairpin(this->_dna_base_num/2);
}

bool DNAString::is_hairpin(size_t len) const{
	if(len==0){
		len = this->_dna_base_num / 2;
	}else{
		assert(len <= this->_dna_base_num / 2);
	}
	size_t forward_block_idx, backward_block_idx;
	size_t forward_offset, backward_offset;//index(0-3) of slot to compare inside each data block.
	forward_block_idx = 0; backward_block_idx = this->_data.size()-1;
	forward_offset = 0; backward_offset = (this->_dna_base_num-1) % 4;
	for(int x = 0; x < len; x++){
		// cout<<x<<":"<<endl;
		// cout<<forward_block_idx<<":"<<forward_offset<<":"<<DNA::bits2base[DNA::base_mask[forward_offset] & this->_data[forward_block_idx]]<<endl;
		// cout<<backward_block_idx<<":"<<backward_offset<<":"<<DNA::bits2RCbase[DNA::base_mask[backward_offset] & this->_data[backward_block_idx]]<<endl;
		if(DNA::bits2base[DNA::base_mask[forward_offset] & this->_data[forward_block_idx]] != DNA::bits2RCbase[DNA::base_mask[backward_offset] & this->_data[backward_block_idx]]){
			return false;
		}
		if(forward_offset == 3){
			forward_block_idx++;
			forward_offset=0;	
		}else{
			forward_offset++;
		}
		if(backward_offset == 0){
			backward_block_idx--;
			backward_offset=3;
		}else{
			backward_offset--;
		}
	}
	return true;
}

//return true if kmer satisfy (A|T|G|C)*
bool DNAString::is_simple() const{
	assert(this->_data.size()>0);

	uint8_t repeats[4] = {0x00, 0x55, 0xAA, 0xFF};
	size_t block_num = this->_data.size();
	if(block_num==1){
		for(int x=0; x < 4; x++){
			if(!((repeats[x] ^ this->_data[0]) >> (this->_offset+2))){
				return true;
			}
		}
	}else{
		int pattern_idx = 0;
		bool flag=false;
		for(; pattern_idx<4; pattern_idx++){
			if(!(repeats[pattern_idx] ^ this->_data[0])){
				flag=true;
				break;
			}
		}
		if(flag){
			size_t block_idx=1;
			while(block_idx < block_num-1){
				if(repeats[pattern_idx] != this->_data[block_idx]){
					return false;
				}
				block_idx++;
			}
			if(!((repeats[pattern_idx] ^ this->_data[block_idx]) >> (this->_offset+2))){
				return true;
			}
		}
	}
	return false;
}

bool operator== (const DNAString& lhs, const DNAString& rhs){
	//cout<<"called"<<endl;
	if(lhs._data == rhs._data && lhs._dna_base_num == rhs._dna_base_num){
		return true;
	}
	return false;
}

bool operator!= (const DNAString& lhs, const DNAString& rhs){
	return !(lhs == rhs);
}

bool operator< (const DNAString& lhs, const DNAString& rhs){
	if(lhs._data < rhs._data){
		return true;
	}
	return false;
}

std::ostream& operator<<(std::ostream &output, const DNAString& dnaStr){
	// output<<dnaStr._data.size()<<endl;
	// for(const auto ele:dnaStr._data){
	// 	output<<ele<<endl;
	// }
	output<<dnaStr.to_str();
	return output;
}

DNAString operator + (const DNAString& str, const char& base){
	DNAString tmp(str);
	tmp.append(base);
	return tmp;
}
DNAString operator + (const char& base, const DNAString& str){
	DNAString tmp(string(1, base) + str.to_str());
	return tmp;
}

//canonical nthash
uint64_t NTPC64(const DNAString& dnaStr, const unsigned k, uint64_t& fhVal, uint64_t& rhVal) {
	fhVal=0, rhVal=0;
	const uint8_t* data = dnaStr.data();
	size_t slot = 0;
  for(unsigned i=0; i<dnaStr.length(); ) {
  	for(int x = 0; x<4 && i < dnaStr.length(); x++, i++){
  		char c = DNA::bits2base[data[slot]&DNA::base_mask[x]];
    	fhVal ^= msTab[(unsigned char)c][(k-1-i)%64];
    	rhVal ^= msTab[(unsigned char)c & cpOff][i%64];
  	}
  	slot++;
  }
  return (rhVal<fhVal)? rhVal : fhVal;
}