// Copyright (c) 2019 Wen Luo <wenluo2016@u.northwestern.edu>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "functions_Random_vec.h"
#include "erf_inv.h"

#include <random>
#include <chrono>
#include <cmath>

void SetUniformDistribution(const double &lower_bd, const double &upper_bd,
                              vector<double> &vec){
    default_random_engine generator;
    std::uniform_real_distribution<double> distribution(lower_bd, upper_bd);
    for (auto it = vec.begin(); it != vec.end(); ++it){
        (*it) = distribution(generator);
    }
}

void SetBoundariesWeibull(const double &w_modulus, const double &x_0,
                              vector<double> &vec_t, vector<double> &vec_c,
                              const double &ratio_c_t){
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    std::weibull_distribution<double> distribution(w_modulus, x_0);
    for (size_t i = 0; i != vec_t.size(); ++i){
        vec_t[i] = distribution(generator);
        vec_c[i] = - ratio_c_t * vec_t[i];
    }
}


/** Weibull(modulus=10, x0=100), Gaussian(mean=100, StdDev=16) grafted at x=50.12122, F=9.99996e-4 **/
void SetBoundariesGaussWeibull(vector<double> &vec_t, vector<double> &vec_c,
                               const double &ratio_c_t){
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution(0, 1);
    for (size_t i = 0; i != vec_t.size(); ++i){
        vec_t[i] = distribution(generator);
        if (vec_t[i] > 9.99996e-4)
            vec_t[i] = 100 - 16*sqrt(2)*erfinv<double> (1-(vec_t[i]-8.78845006974339e-5)/4.999560577496512e-1 );
        else
            vec_t[i] = 100 * exp(0.1*log(-log(1-vec_t[i])));
        vec_c[i] = - ratio_c_t * vec_t[i];
    }
}
