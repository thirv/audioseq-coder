#pragma once
#include <vector>
#include <numeric>
#include <cmath>    
#include <cstdlib>

/*
 Pyramid Vector Quantization
 */

unsigned long long pvq_N2(unsigned int l,unsigned int K);

unsigned long long pvq_N3(unsigned int l, unsigned int K, std::vector<std::vector<unsigned long long>> out);

std::vector<std::vector<unsigned long long>> pvq_Npre(unsigned int l, std::vector<unsigned int>);

std::vector<int> pvq_quantize(std::vector<float> x,unsigned int k);

unsigned long long pvq_encode(const std::vector<int> &x,unsigned int l,unsigned int k,const std::vector<std::vector<unsigned long long>> &pvqN);

std::vector<int> pvq_decode(unsigned long long b,unsigned int l,unsigned int k,const std::vector<std::vector<unsigned long long>> &pvqN);