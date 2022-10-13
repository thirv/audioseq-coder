#pragma once
#include <vector>
#include <fstream>
#include <cmath>

/*
 Bitsream read/write
 */

void write_frbs(const std::vector<bool> &in,std::ofstream &ofs);

std::vector<bool> read_frbs(std::ifstream &ifs, unsigned int num_bytes);