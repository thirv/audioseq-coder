#include "filebs.hpp"
using namespace std;


void write_frbs(const vector<bool> &in, ofstream &ofs) {
	//auto num_bytes = (unsigned int)ceil(in.size()/8.0);
	//ofs.write(reinterpret_cast<char*>(&num_bytes),sizeof(num_bytes));
	char byte = 0;
	auto i = 0;
	for(auto bit:in) {
		if(bit)
			byte |= 1<<i;
		i++;
		if(i==8) {
			ofs.put(byte);
			byte = 0;
			i = 0;
		}
	}
	if(i>0)
		ofs.put(byte);
}


vector<bool> read_frbs(ifstream &ifs, unsigned int num_bytes) {
	//unsigned int num_bytes;
	//ifs.read(reinterpret_cast<char*>(&num_bytes),sizeof(num_bytes)); 
	char byte;
	vector<bool> bs(8*num_bytes);
	for(unsigned int i = 0;i<num_bytes;i++) {
		ifs.get(byte);
		for(unsigned int i2 = 0;i2<8;i2++)
			bs[8*i+i2] = (byte & (1<<i2))!=0;
	}
	return bs;
}