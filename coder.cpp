#include "coder.hpp"
#include "mdct.hpp"
#include "pvq.hpp"
#include "bs.hpp"
#include <functional>
#include <algorithm>
using namespace std;


settings init_coder(unsigned int frame_size, unsigned int fs)
{
	settings s;
	s.mdct_win = generate_win(2*frame_size);
	s.frame_size = frame_size;
	s.fs = fs;
	s.mdct_size = 2 * frame_size;
	s.tw = mdct_init(2*frame_size);
	s.blims = {0,4,8,12,16,24,32,40,56,72,88,104,120,136,152,168,184,200,216,232,248,264,280,296,
		   312,328,344,360,376,392,408,424,440,544,672,1024};
    	for (unsigned int i=1; i<s.blims.size(); i++){
        	s.bwidths.push_back(s.blims[i]-s.blims[i-1]);
	}
    	s.prek = {16,14,12,10,16,14,12,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,
		  16,16,16,16,0,0,0};
	s.nro_bands = (unsigned int)s.blims.size() - 1;
	s.band_spread = {0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,0,0,0}; // not used yet, how many neighbor bands are grouped together in sperad analysis
	s.copy_lim = 32;//; // 18 35 copy bins above and including this band (0-based index)
	s.nocopy_lim = 4;//; // never copy below and including this band (0-based index)
	s.target_rate = 64000; //nominal target in bits/s TODO as input
	s.bytes_fr = ceil(s.target_rate*s.frame_size/fs/8.0);
	s.rate_fr = s.bytes_fr * 8;
	s.pvqN = pvq_Npre(26,vector<unsigned int>(27,26));//pvq_Npre(16,{16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16});
	s.frew = get_perw(frame_size, fs);//.resize(frame_size,1.0f);//
	return s;
}


void reset_inbuff(unsigned int nro_ch, settings &s)
{
    s.nro_ch = nro_ch;
    s.in_buffer.resize(nro_ch,vector<float>(2*s.frame_size,0));
    s.pre_norms_enc.resize(nro_ch,vector<int>(s.nro_bands,-60));
    s.fr = 0;
}


void reset_outbuff(unsigned int nro_ch, settings &s)
{
    s.nro_ch = nro_ch;
    s.out_buffer.resize(nro_ch,vector<float>(s.frame_size,0));
    s.pre_norms_dec.resize(nro_ch,vector<int>(s.nro_bands,-60));
    s.fr = 0;
}


vector<float> get_perw(unsigned int frame_size, unsigned int fs)
{
	vector<float> w(frame_size);
	for(auto i = 0; i<frame_size; i++) {
		auto bin_fre = (i+0.5f) * fs / 2 / frame_size * 0.001f;
		/*auto hear_thres = (3.64f*pow(bin_fre,-0.8f))-(6.5f*exp(-0.6f*pow(bin_fre-3.3f,2.0f)))+(0.001f*pow(bin_fre,4.0f));*/ // official low cutting
		auto hear_thres = (3.64f*pow(bin_fre,-0.4f)) - (4.0f*exp(-0.6f*pow(bin_fre,8.0f))) - 
			              (6.5f*exp(-0.6f*pow(bin_fre-3.3f,2.0f))) + (0.001f*pow(bin_fre,4.0f));
		w[i] = 1 / pow(10,(min(40.0f,hear_thres)/20));
	}
	auto max_val = *max_element(w.begin(),w.end());
	for(auto i:w)
		i /= max_val;
	return w;
}


vector<bool> encoder(const vector<vector<float>> &pcm, settings &s)
{
	vector<vector<float>> fre(s.nro_ch,vector<float>(s.frame_size));
	for(unsigned int i = 0; i<s.nro_ch; i++) {
		copy(pcm[i].begin(),pcm[i].end(),s.in_buffer[i].begin()+s.frame_size);
		transform(s.in_buffer[i].begin(),s.in_buffer[i].end(),
			s.mdct_win.begin(),s.in_buffer[i].begin(),multiplies<float>());
		fre[i] = mdct(s.in_buffer[i],s.tw,s.mdct_size);
		transform(fre[i].begin(),fre[i].end(),s.frew.begin(),fre[i].begin(),multiplies<float>());
		copy(pcm[i].begin(),pcm[i].end(),s.in_buffer[i].begin());
	}
	return generate_frbs(fre,s);
}


vector<vector<float>> decoder(vector<bool> fr_bs, settings &s)
{
	auto fre = parse_frbs(fr_bs,s);
	vector<vector<float>> pcm(s.nro_ch,vector<float>(s.frame_size));
	for(unsigned int i = 0; i < s.nro_ch; i++) {
		transform(fre[i].begin(),fre[i].end(),s.frew.begin(),fre[i].begin(),divides<float>());
		auto mdct_out = imdct(fre[i],s.tw,s.mdct_size);
		transform(mdct_out.begin(),mdct_out.end(),
			s.mdct_win.begin(),mdct_out.begin(),multiplies<float>());
		transform(mdct_out.begin(),mdct_out.begin()+s.frame_size,
			s.out_buffer[i].begin(),s.out_buffer[i].begin(),plus<float>());
		pcm[i] = s.out_buffer[i];
		copy(mdct_out.begin()+s.frame_size,mdct_out.end(),s.out_buffer[i].begin());
	}
	return pcm;
}
