#pragma once

#include <fstream>
#include "settings.hpp"

/*
 Bitsream Generation
 */

std::vector<bool> generate_frbs(const std::vector<std::vector<float>> &fre, settings &s);

std::vector<std::vector<float>> parse_frbs(const std::vector<bool> &fr_bs, settings &s);

std::vector<unsigned int> adjust_k(unsigned int target_rate, const settings &s);

std::vector<std::vector<float>> get_bandrates(const std::vector<std::vector<float>> &enes,const settings &s);

int get_k(float rate, unsigned int l, const std::vector<std::vector<unsigned long long>> &pvqN);