#pragma once
#include "settings.hpp"

/*
 Perceptual Coding
 */

settings init_coder(unsigned int frame_size, unsigned int fs);

void reset_inbuff(unsigned int nro_ch, settings &s);

void reset_outbuff(unsigned int nro_ch, settings &s);

std::vector<float> get_perw(unsigned int frame_size, unsigned int fs);

std::vector<bool> encoder(const std::vector<std::vector<float>> &pcm, settings &s);

std::vector<std::vector<float>> decoder(std::vector<bool> fr_bs, settings &s);

