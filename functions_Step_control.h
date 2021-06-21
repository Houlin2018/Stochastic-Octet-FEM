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

#ifndef FUNCTIONS_STEP_CONTROL_H_INCLUDED
#define FUNCTIONS_STEP_CONTROL_H_INCLUDED

#include <vector>
#include <Eigen/Sparse>

using namespace std;
using std::vector;

void GetLoadRatio_Elastic(const vector<double> &bound_plus,
                          const vector<double> &bound_minus,
                          const vector<double> &res_Ke,
                          const vector<double> &stress,
                          vector<double> &load_ratio);

inline void MulVec(vector<double> &vec, const double &mult){
    for(auto it = vec.begin(); it != vec.end(); ++it)
        *it *= mult;
}

inline void AddVec(vector<double> &vec, const vector<double> &vec_add){
    for(size_t i = 0; i != vec.size(); ++i)
        vec[i] += vec_add[i];
}

void UpdateSPK(const size_t &c_ele,
               const double &len,
               const double &k_old, const double &k_new,
               const vector<vector<double> > &coords,
               const vector<vector<size_t> > &conn_flat,
               Eigen::SparseMatrix<double> &spK_glob);



#endif // FUNCTIONS_STEP_CONTROL_H_INCLUDED
