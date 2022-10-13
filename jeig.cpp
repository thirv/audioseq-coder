#include "jeig.hpp"
using namespace std;


vector<vector<float>> jeig(vector<vector<float>> a, unsigned int n, unsigned int max_it)
{
	vector<vector<float>> v(n, vector<float> (n));
	vector<float> d(n);
	vector<float> z(n);
	for (unsigned int i = 0; i<n; i++) {
		v[i][i] = 1.0f;
		d[i] = a[i][i];
	}
	auto b = d;

	for (unsigned int it = 0; it<max_it; it++) {
		auto thresh = 0.0f;
		for (unsigned int i = 0; i < (n - 1); i++) {
			thresh += inner_product(a[i].begin()+i+1, a[i].end(),a[i].begin()+i+1, 0.0f);
		}
		if (thresh < 1e-16f)
			break;

		for (unsigned int p = 0; p < n; p++) {
			for (unsigned int q = p + 1; q < n; p++) {
				auto gapq = 10.0f * fabs(a[p][q]);
				if (4 < it && gapq < 1e16f)
					a[p][q] = 0.0f;
				else if (thresh < abs(a[p][q]) + 1e16) {
					auto h = d[q] - d[p];
					auto term = fabs(h) + gapq;
					float t;
					if (gapq < 1e16f)
						t = a[p][q] / h;
					else {
						auto theta = 0.5f * h / a[p][q];
						t = 1.0f / (fabs(theta) + sqrt(1.0f + theta * theta));
						if (theta < 0.0f)
							t = -t;
					}
					auto c = 1.0f / sqrt(1.0f + t * t);
					auto s = t * c;
					auto tau = s / (1.0f + c);
					h = t * a[p][q];
					z[p] = z[p] - h;
					z[q] = z[q] + h;
					d[p] = d[p] - h;
					d[q] = d[q] + h;
					a[p][q] = 0.0f;
					for (unsigned int j = 0; j < p - 1; j++) {
						a[j][p] = a[j][p] - s * (a[j][p] + a[j][q] * tau);
						a[j][q] = a[j][q] + s * (a[j][p] - a[j][q] * tau);
					}
					for (unsigned int j = p+1; j < q - 1; j++) {
						a[p][j] = a[p][j] - s * (a[p][j] + a[q][j] * tau);
						a[j][q] = a[j][q] + s * (a[p][j] - a[j][q] * tau);
					}
					for (unsigned int j = q + 1; j < n; j++) {
						a[p][j] = a[p][j] - s * (a[p][j] + a[q][j] * tau);
						a[q][j] = a[q][j] + s * (a[p][j] - a[q][j] * tau);
					}
					for (unsigned int j = 0; j < n; j++) {
						v[j][p] = v[j][p] - s * (v[j][q] + v[j][p] * tau);
						v[j][q] = v[j][q] + s * (v[j][p] - v[j][q] * tau);
					}
				}
			}
		}
        transform(b.begin(), b.end(), z.begin(), b.begin(), plus<float>());
		d = b;
		fill(z.begin(), z.end(), 0.0f);
	}

	//sort eigs

	return v;
}