#include "mdct.hpp"
using namespace std;


vector<float> generate_win(unsigned int size)
{
	vector<float> win(size);
	for (unsigned int i=0; i<size; i++)
		win[i] = (float) sin(M_PI/size*(i+0.5));
	return win;
}


void fft(vector<complex<float>>& zs,int sign)
{
	auto s = zs.size();
	auto s2 = s/2;
	unsigned int j = 0;
	auto pis = M_PI * sign;

	for(unsigned int i = 0; i<s-1; i++) {
		if(i<j)
			swap(zs[i],zs[j]);
		auto m = s2;
		j ^= m;
		while((j&m)==0) {
			m /= 2;
			j ^= m;
		}
	}

	for(unsigned int j = 1; j<s; j *= 2)
		for(unsigned int m = 0; m<j; m++) {
			auto t = pis * m/j;
			complex<float> w((float)cos(t),(float)sin(t));
			for(unsigned int i = m;i<s;i += 2*j) {
				auto zi = zs[i];
				auto t = w * zs[i+j];
				zs[i] = zi+t;
				zs[i+j] = zi-t;
			}
		}
}


vector<complex<float>> mdct_init(int mdct_size)
{
	vector<complex<float>> twiddle(mdct_size/4);
	auto omega = 2 * M_PI / mdct_size;
	auto alpha = omega/8;
	auto scale = sqrt(sqrt(2.0f/mdct_size));
	for(auto i = 0; i<mdct_size/4; i++) {
		twiddle[i] = complex<float>(scale*(float)cos(omega*i+alpha),-scale*(float)sin(omega*i+alpha)); //- sign eases complext product in use
	}
	return twiddle;
}


vector<float> mdct(const vector<float> &pcm,const vector<complex<float>> &tw,unsigned int mdct_size)
{
	int n4 = mdct_size / 4;
	int n2 = 2*n4;
	int n34 = 3*n4;
	int n54 = 5*n4;
	vector<complex<float>> x(n4);
	vector<float> out(n2);

	for(int n = 0; n<n4; n += 2)
		x[n/2] = complex<float>(pcm[n34-1-n]+pcm[n34+n],pcm[n4+n]-pcm[n4-1-n]) * tw[n/2]; //x*conj(y)
	for(int n = n4; n<n2; n += 2)
		x[n/2] = complex<float>(pcm[n34-1-n]-pcm[-n4+n],pcm[n4+n]+pcm[n54-1-n]) * tw[n/2];

	fft(x,-1);

	for(int n = 0; n<n2; n += 2) {
		auto w = x[n/2]*tw[n/2];
		out[n] = -real(w);
		out[n2-1-n] = imag(w);
	}

	return out; // scale by sqrt(n2)?
}


vector<float> imdct(const vector<float> &mdct_in,const vector<complex<float>>& tw,unsigned int mdct_size)
{
	int n4 = mdct_size/4;
	int n2 = 2*n4;
	int n34 = 3*n4;
	int n54 = 5*n4;
	vector<complex<float>> x(n4);
	vector<float> out(mdct_size);

	for(auto n = 0;n<n2;n += 2)
		x[n/2] = complex<float>(-2,0) * complex<float>(mdct_in[n],mdct_in[n2-1-n]) * tw[n/2]; //x*conj(y)// todo inv freq weighting	

	fft(x,-1);

	for(auto n = 0;n<n4;n += 2) {
		auto w = x[n/2]*tw[n/2];
		float r1 = real(w);
		float i1 = -imag(w);
		out[n34-1-n] = r1;
		out[n34+n] = r1;
		out[n4+n] = i1;
		out[n4-1-n] = -i1;
	}
	for(auto n = n4;n<n2;n += 2) {
		auto w = x[n/2]*tw[n/2];
		float r1 = real(w);
		float i1 = -imag(w);
		out[n34-1-n] = r1;
		out[-n4+n] = -r1;
		out[n4+n] = i1;
		out[n54-1-n] = i1;
	}

	return out;

}