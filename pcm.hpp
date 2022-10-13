#pragma once
#include <vector>
#include <complex>
#include <cmath>
#include <iostream>
#include <fstream>

/*
 Wav read/write
 */

struct pcm_header {
    std::vector<std::uint8_t> bytes;
    unsigned int nro_ch, fs, nro_bytes, nro_fr, last_size;
};

pcm_header get_header(unsigned int frame_size, std::ifstream &f);

std::vector<std::uint8_t> header_bytes(unsigned int nro_ch, unsigned int fs, unsigned int nro_frames, unsigned int frame_size);

std::vector<std::vector<float>> get_pcm(unsigned int nro_ch, unsigned int frame_size, std::ifstream &f);

void put_pcm(const std::vector<std::vector<float>> &pcm, unsigned int nro_ch, unsigned int frame_size, std::ofstream &f);