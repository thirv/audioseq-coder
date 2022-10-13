#pragma once
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <complex>

/*
 Modified Discrete Cosine Transform
 */

std::vector<float> generate_win(unsigned int size);

void fft(std::vector<std::complex<float>>& zs,int sign);

std::vector<std::complex<float>> mdct_init(int mdct_size);

std::vector<float> mdct(const std::vector<float> &pcm,const std::vector<std::complex<float>> &tw,unsigned int mdct_size);

std::vector<float> imdct(const std::vector<float> &mdct_in,const std::vector<std::complex<float>> &tw,unsigned int frame_size);