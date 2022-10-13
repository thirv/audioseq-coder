#pragma once
#include <vector>
#include <numeric>
#include <algorithm>
#include <functional>
#include <cmath>

/*
 Jacobi Eigenvalue Algoritm
 */

std::vector<std::vector<float>> jeig(std::vector<std::vector<float>> a, unsigned int n, unsigned int max_it);
