#include "pvq.hpp"
using namespace std;


unsigned long long pvq_N2(unsigned int l,unsigned int K)
{
	if(K==0)
		return 1;
	else if(l==0)
		return 0;
	else if(l==1)
		return 2;
	else if(K==1)
		return 2*l;
	else
		return pvq_N2(l-1,K)+pvq_N2(l-1,K-1)+pvq_N2(l,K-1);
}


unsigned long long pvq_N3(unsigned int l, unsigned int K, vector<vector<unsigned long long>> out)
{
	if(K==0)
		return 1;
	else if(l==0)
		return 0;
	else if(l==1)
		return 2;
	else if(K==1)
		return 2*l;
	else {
		unsigned long long v1, v2, v3;
		v1 = K<out[l-1].size() ? out[l-1][K] : pvq_N3(l-1,K,out);
		v2 = (K-1)<out[l-1].size() ? out[l-1][K-1] : pvq_N3(l-1,K-1,out);
		v3 = (K-1)<out[l].size() ? out[l][K-1] : pvq_N3(l,K-1,out);
		return v1 + v2 + v3;
	}
}


vector<vector<unsigned long long>> pvq_Npre(unsigned int l, vector<unsigned int> k)
{
	vector<vector<unsigned long long>> out(l+1);
	for(unsigned int li=0;li<=l;li++){
		for(unsigned int ki=0;ki<=k[li];ki++){
			out[li].push_back(pvq_N3(li,ki,out));
        	}
    	}
    	return out;
}


vector<int> pvq_quantize(vector<float> x,unsigned int K)
{
	auto n1 = accumulate(x.begin(),x.end(),1e-30f,[](float tot,float i) {return tot+fabs(i);});
	for(auto &i:x)
		i *= K / n1;
	vector<int> xr(x.size());
	for(int i = 0;i<x.size();i++)
		xr[i] = (int)round(x[i]);
	auto tk = (unsigned int)accumulate(xr.begin(),xr.end(),0,[](int tot,int i) {return tot+abs(i);});
	while(tk<K) {
		auto mi = 0;
		auto me = fabs(x[0]) - fabs(xr[0]);
		for(int i = 1;i<x.size();i++) {
			auto met = fabs(x[i]) - fabs(xr[i]);
			if(me<met) {
				mi = i;
				me = met;
			}
		}
		xr[mi] += ((0<x[mi])-(x[mi]<0));
		tk += 1;
	}
	while(tk>K) {
		auto mi = 0;
		auto me = fabs(x[0]) - fabs(xr[0]);
		for(int i = 1;i<x.size();i++) {
			auto met = fabs(x[i]) - fabs(xr[i]);
			if(me>met) {
				mi = i;
				me = met;
			}
		}
		xr[mi] -= ((0<x[mi])-(x[mi]<0));
		tk -= 1;
	}
	return xr;
}


unsigned long long pvq_encode(const vector<int> &x,unsigned int l,unsigned int k,const vector<vector<unsigned long long>> &pvqN)
{
	unsigned long long b = 0;
	auto i = 0;
	while(k>0) {
		auto xa = abs(x[i]);
		auto xs = (0 < x[i])-(x[i] < 0); // sign function
		if(xa==1)
			b += pvqN[l-1][k] + (1-xs) / 2 * pvqN[l-1][k-1];
		else if(xa>1) {
			unsigned long long tmp = 0;
			b += 2 * accumulate(next(pvqN[l-1].begin(),k-xa+1),next(pvqN[l-1].begin(),k),tmp);
			b += pvqN[l-1][k] + (1-xs) / 2 * pvqN[l-1][k-xa];
		}
		k -= xa;
		l -= 1;
		i += 1;
	}
	return b;
}


vector<int> pvq_decode(unsigned long long b,unsigned int l,unsigned int k,const vector<vector<unsigned long long>> &pvqN)
{
	vector<int> x(l);
	unsigned long long xb = 0;
	auto i = 0;
	while(k>0) {
		if(b==xb) {
			x[i] = 0;
			break;
		}
		else if(b<xb+pvqN[l-1][k]) {
			x[i] = 0;
		}
		else {
			xb += pvqN[l-1][k];
			unsigned int j = 1;
			while(b>=xb+2*pvqN[l-1][k-j]) {
				xb += 2*pvqN[l-1][k-j];
				j += 1;
			}
			if(xb<=b && b<xb+pvqN[l-1][k-j])
				x[i] = (int)j;
			else if(b>=xb+pvqN[l-1][k-j]) {
				x[i] = -(int)j;
				xb += pvqN[l-1][k-j];
			}
			k -= j;
		}
		l -= 1;
		i += 1;
	}
	if(k>0)
		x.back() = k-abs(x[i]);
	return x;
}
