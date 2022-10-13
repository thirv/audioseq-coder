#pragma once
#include <vector>
#include <complex>

struct settings {
	std::vector<std::vector<float>> in_buffer;
	std::vector<std::vector<float>> out_buffer;
	std::vector<float> mdct_win;
	std::vector<std::complex<float>> tw;
	unsigned int nro_ch;
	unsigned int frame_size;
	unsigned int fs;
	unsigned int mdct_size;
    unsigned int nro_bands;
    std::vector<unsigned int> blims;
    std::vector<unsigned int> bwidths;
    std::vector<unsigned int> prek;
	std::vector<unsigned int> band_spread;
	unsigned int copy_lim;
	unsigned int nocopy_lim;
	unsigned int target_rate;
    unsigned int bytes_fr;
	unsigned int rate_fr;
	std::vector<std::vector<int>> pre_norms_enc;
	std::vector<std::vector<int>> pre_norms_dec;
    std::vector<std::vector<unsigned long long>> pvqN;
	unsigned int fr;
	std::vector<float> frew;
};