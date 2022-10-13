#include "bs.hpp"
#include "pvq.hpp"
#include <cmath>
#include <numeric>
#include <algorithm>
#include <map>
using namespace std;

const map<int,vector<bool>> enchuff = {
	{-20,{1,1,0,1,1,0}},
	{-19,{1,1,0,0,0,0}},
	{-18,{1,0,0,1,1,0}},
	{-17,{1,0,0,0,0,0}},
	{-16,{0,1,0,1,1,0}},
	{-15,{0,1,0,0,0,0}},
	{-14,{0,0,0,1,1,0}},
	{-13,{0,0,0,0,1,0}},
	{-12,{0,0,0,0,0,0}},
	{-11,{1,1,1,1,0}},
	{-10,{1,1,1,0,0}},
	{-9,{1,1,0,0,1}},
	{-8,{1,0,1,1,0}},
	{-7,{1,0,1,0,0}},
	{-6,{1,0,0,0,1}},
	{-5,{0,1,1,1,0}},
	{-4,{0,1,1,0,0}},
	{-3,{0,1,0,0,1}},
	{-2,{0,0,1,1,0}},
	{-1,{0,0,1,0,0}},
	{0,{0,0,0,1,0}},
	{1,{0,0,1,0,1}},
	{2,{0,0,1,1,1}},
	{3,{0,1,0,1,0}},
	{4,{0,1,1,0,1}},
	{5,{0,1,1,1,1}},
	{6,{1,0,0,1,0}},
	{7,{1,0,1,0,1}},
	{8,{1,0,1,1,1}},
	{9,{1,1,0,1,0}},
	{10,{1,1,1,0,1}},
	{11,{1,1,1,1,1}},
	{12,{0,0,0,0,0,1}},
	{13,{0,0,0,0,1,1}},
	{14,{0,0,0,1,1,1}},
	{15,{0,1,0,0,0,1}},
	{16,{0,1,0,1,1,1}},
	{17,{1,0,0,0,0,1}},
	{18,{1,0,0,1,1,1}},
	{19,{1,1,0,0,0,1}},
	{20,{1,1,0,1,1,1}}};

const map<vector<bool>,int> dechuff5 = {
	{{1,1,1,1,0},-11},
	{{1,1,1,0,0},-10},
	{{1,1,0,0,1},-9},
	{{1,0,1,1,0},-8},
	{{1,0,1,0,0},-7},
	{{1,0,0,0,1},-6},
	{{0,1,1,1,0},-5},
	{{0,1,1,0,0},-4},
	{{0,1,0,0,1},-3},
	{{0,0,1,1,0},-2},
	{{0,0,1,0,0},-1},
	{{0,0,0,1,0},0},
	{{0,0,1,0,1},1},
	{{0,0,1,1,1},2},
	{{0,1,0,1,0},3},
	{{0,1,1,0,1},4},
	{{0,1,1,1,1},5},
	{{1,0,0,1,0},6},
	{{1,0,1,0,1},7},
	{{1,0,1,1,1},8},
	{{1,1,0,1,0},9},
	{{1,1,1,0,1},10},
	{{1,1,1,1,1},11}};

const map<vector<bool>,int> dechuff6 = {
	{{1,1,0,1,1,0},-20},
	{{1,1,0,0,0,0},-19},
	{{1,0,0,1,1,0},-18},
	{{1,0,0,0,0,0},-17},
	{{0,1,0,1,1,0},-16},
	{{0,1,0,0,0,0},-15},
	{{0,0,0,1,1,0},-14},
	{{0,0,0,0,1,0},-13},
	{{0,0,0,0,0,0},-12},
	{{0,0,0,0,0,1},12},
	{{0,0,0,0,1,1},13},
	{{0,0,0,1,1,1},14},
	{{0,1,0,0,0,1},15},
	{{0,1,0,1,1,1},16},
	{{1,0,0,0,0,1},17},
	{{1,0,0,1,1,1},18},
	{{1,1,0,0,0,1},19},
	{{1,1,0,1,1,1},20}};


vector<bool> generate_frbs(const vector<vector<float>> &fre, settings &s)
{
	vector<vector<int>> env(s.nro_ch,vector<int>(s.nro_bands));
	vector<vector<float>> enes(s.nro_ch,vector<float>(s.nro_bands));
	vector<bool> out;
    out.reserve(s.rate_fr);
	//auto ene_tot = 1e-30f;
	for(unsigned int ch = 0; ch<s.nro_ch; ch++) {
		for(unsigned int b = 0; b<s.nro_bands; b++) {
			env[ch][b] = (int)round(10*log10(inner_product(fre[ch].begin()+s.blims[b],
				fre[ch].begin()+s.blims[b+1],fre[ch].begin()+s.blims[b],0.0f)));
			auto diff_db = min(20,max(-20,env[ch][b]-s.pre_norms_enc[ch][b]));
			auto cw = enchuff.at(diff_db);
			out.insert(out.end(),cw.begin(),cw.end());
			s.pre_norms_enc[ch][b] += diff_db;
			enes[ch][b] = pow(10.0f,(s.pre_norms_enc[ch][b])/10.0f); // use quantized energies for calculations
			//ene_tot += enes[ch][b];
		}
	}
    auto newk = adjust_k(floor((s.rate_fr-out.size())/s.nro_ch), s);
	for(unsigned int ch = 0; ch<s.nro_ch; ch++) {
		for(unsigned int b = 0; b<s.nro_bands; b++) {
			vector<float> band(fre[ch].begin()+s.blims[b],fre[ch].begin()+s.blims[b+1]);
			auto l = (unsigned int)band.size();
			if(b<s.copy_lim) { // todo some other (copylocation?) analysis for higher bands
				if(newk[b]>l/2 | b<=s.nocopy_lim) { // dont use too few pulses
					auto val = pvq_encode(pvq_quantize(band,newk[b]),l,newk[b],s.pvqN);
					for(auto i = 0; i<ceil(log2(s.pvqN[l][newk[b]])); i++) {
						out.push_back(val&1);
						val >>= 1;
					}
				}
			}
		}
	}
    while (out.size()<s.rate_fr)
        out.push_back(0); // make sure out size has desired num of bytes
	return out;
}


vector<vector<float>> parse_frbs(const vector<bool> &fr_bs, settings &s)
{
	vector<vector<int>> env(s.nro_ch,vector<int>(s.nro_bands));
	vector<vector<float>> enes(s.nro_ch,vector<float>(s.nro_bands));
	//auto ene_tot = 0.0f;
	unsigned int bi = 0;
	for(unsigned int ch = 0; ch<s.nro_ch; ch++) {
		for(unsigned int b = 0; b<s.nro_bands; b++) {
			vector<bool> cw5(fr_bs.begin()+bi,fr_bs.begin()+bi+5);
			auto it = dechuff5.find(cw5);
			if(it!=dechuff5.end())
				bi += 5;
			else {
				vector<bool> cw6(fr_bs.begin()+bi,fr_bs.begin()+bi+6);
				it = dechuff6.find(cw6);
				bi += 6;
			}
			env[ch][b] = it->second+s.pre_norms_dec[ch][b];
			enes[ch][b] = pow(10.0f,env[ch][b]/10.0f);
			//ene_tot += enes[ch][b];
			s.pre_norms_dec[ch][b] = env[ch][b];
		}
	}
    auto newk = adjust_k(floor((s.rate_fr-bi)/s.nro_ch), s);
	vector<vector<float>> fre(s.nro_ch,vector<float>(s.frame_size,0));
	for(unsigned int ch = 0; ch<s.nro_ch; ch++) {
		for(unsigned int b = 0; b<s.nro_bands; b++) {
			auto l = s.blims[b+1]-s.blims[b];
			if(((b<s.copy_lim) & (newk[b]>l/2)) | (b<=s.nocopy_lim)) {
				unsigned long long val = 0;
				for(unsigned long long i = 0; i<ceil(log2(s.pvqN[l][newk[b]])); i++) {
					if(fr_bs[bi])
						val |= (unsigned long long) 1 << i;
					bi += 1;
				}
				auto band = pvq_decode(val,l,newk[b],s.pvqN);
                copy(band.begin(),band.end(),fre[ch].begin()+s.blims[b]);
				//copy(s.dbg_fre[ch].begin()+s.blims[b],s.dbg_fre[ch].begin()+s.blims[b+1],fre[ch].begin()+s.blims[b]);
			}
			else {
				copy(fre[ch].begin()+s.blims[b]-l,fre[ch].begin()+s.blims[b],fre[ch].begin()+s.blims[b]);
				//copy(s.dbg_fre[ch].begin()+s.blims[b],s.dbg_fre[ch].begin()+s.blims[b+1],fre[ch].begin()+s.blims[b]);
			}
			auto scale = sqrt(enes[ch][b]/inner_product(fre[ch].begin()+s.blims[b],
						 fre[ch].begin()+s.blims[b+1],fre[ch].begin()+s.blims[b],1e-30f));
			for(auto i = 0; i<l; i++)
				fre[ch][s.blims[b]+i] *= scale;
		}
	}
	return fre;
}


vector<unsigned int> adjust_k(unsigned int target_rate, const settings &s)
{
    // assume each change is less than 2 bits?
    auto newk = s.prek;
    unsigned int rate = 0;
    for (unsigned int b=0; b<s.copy_lim; b++)
        rate += ceil(log2(s.pvqN[s.bwidths[b]][s.prek[b]]));
    if (target_rate>rate) {
        auto i = 0;
        while (true) {
            rate -= ceil(log2(s.pvqN[s.bwidths[i]][newk[i]])) - ceil(log2(s.pvqN[s.bwidths[i]][newk[i]+1]));
            if (target_rate<rate)
                break;
            newk[i] += 1;
            i = (i+1) % s.copy_lim;
        }
    }
    else if (target_rate<rate) {
        auto i = 0;
        while (true) {
            rate -= ceil(log2(s.pvqN[s.bwidths[s.copy_lim-i-1]][newk[s.copy_lim-i-1]])) -
                ceil(log2(s.pvqN[s.bwidths[s.copy_lim-i-1]][newk[s.copy_lim-i-1]-1]));
            if (target_rate>rate)
                break;
            newk[s.copy_lim-i-1] -= 1;
            i = (i+1) % s.copy_lim;
        }
    }
    return newk;
}


vector<vector<float>> get_bandrates(const vector<vector<float>> &enes, const settings &s) //not used yet/anymore, error from decode
{
	vector<vector<float>> bandrates(s.nro_ch,vector<float>(s.nro_bands));
	vector<vector<float>> w(s.nro_ch,vector<float>(s.nro_bands));
	for(unsigned int ch = 0; ch<s.nro_ch; ch++) {
		auto rate_ch = s.rate_fr / s.nro_ch; //TODO between ch ene dist
		for(unsigned int b = 0; b<s.nro_bands; b++) {
			auto ene_spread = accumulate(enes[ch].begin()+b-s.band_spread[b],enes[ch].begin()+b+1+s.band_spread[b],1e-30f);
			w[ch][b] = enes[ch][b] / ene_spread; 
		}
		auto sum_wr = 1e-30f;
		for(auto i:w[ch])
			sum_wr += i * rate_ch / s.nro_bands;
		for(unsigned int b = 0; b<s.nro_bands; b++) {
			w[ch][b] *= s.rate_fr / sum_wr; // some error here
			bandrates[ch][b] = w[ch][b] * rate_ch / s.nro_bands;
		}
	}
	return bandrates;
}


int get_k(float rate, unsigned int l, const vector<vector<unsigned long long>> &pvqN) //not used yet/anymore
{
	int k = 0;
	while(k<pvqN[l].size()-1) {
		if(rate>log2(pvqN[l][k]))
			k += 1;
		else
			break;
	}
	return k;
}