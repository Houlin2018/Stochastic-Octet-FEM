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

#include "functions_Step_control.h"

#include <stdexcept>
#include <algorithm>
#include <limits>
//#include <iostream>

void GetLoadRatio_Elastic(const vector<double> &bound_plus,
                          const vector<double> &bound_minus,
                          const vector<double> &res_Ke,
                          const vector<double> &stress,
                          vector<double> &load_ratio){
    const double infinity = std::numeric_limits<double>::infinity();
    if( bound_plus.size() != load_ratio.size() ){
            throw std::invalid_argument( "Output vector does NOT have the same size as inputs!" );
    }

    for(size_t i = 0; i != bound_plus.size(); ++i){
        if (res_Ke[i] == 0)
            load_ratio[i] = infinity;
        else{
            if (stress[i] >= 0)
                load_ratio[i] = bound_plus[i] / stress[i];
            else
                load_ratio[i] = bound_minus[i] / stress[i];
        }
    }
    for(size_t i = 0; i != load_ratio.size(); ++i){
        if (load_ratio[i] <= 0)
            throw std::invalid_argument( "load_ratio contains NEGATIVE values. Results are INVALID!" );
    }
}

void UpdateSPK(const size_t &c_ele,
               const double &len,
               const double &k_old, const double &k_new,
               const vector<vector<double> > &coords,
               const vector<vector<size_t> > &conn_flat,
               Eigen::SparseMatrix<double> &spK_glob){
    size_t r_ind, c_ind;
    double cx, cy, cz;
    cx = (coords[ conn_flat[c_ele][1] ][0] - coords[ conn_flat[c_ele][0] ][0]) / len;
    cy = (coords[ conn_flat[c_ele][1] ][1] - coords[ conn_flat[c_ele][0] ][1]) / len;
    cz = (coords[ conn_flat[c_ele][1] ][2] - coords[ conn_flat[c_ele][0] ][2]) / len;
    vector<double> dir_cos{cx, cy, cz, -cx, -cy, -cz};
    /* update entries */
    double Ke;
    for (size_t node_i = 0; node_i != 2; ++node_i){ // node_i, node_j: node loops
        for (size_t node_j = 0; node_j != 2; ++node_j){
            /* row/column index of the spK_glob */
            r_ind = 3 * conn_flat[c_ele][node_i];
            c_ind = 3 * conn_flat[c_ele][node_j];
            /* construct spK_glob */
            Ke = (k_new - k_old) * dir_cos[3*node_i]*dir_cos[3*node_j];
            if( Ke != 0)
                    spK_glob.coeffRef(r_ind, c_ind) += Ke;

            Ke = (k_new - k_old) * dir_cos[3*node_i+1]*dir_cos[3*node_j];
            if( Ke != 0){
                    spK_glob.coeffRef(r_ind+1, c_ind) += Ke;
                    spK_glob.coeffRef(r_ind, c_ind+1) += Ke;
            }

            Ke = (k_new - k_old) * dir_cos[3*node_i+2]*dir_cos[3*node_j];
            if( Ke != 0){
                    spK_glob.coeffRef(r_ind+2, c_ind) += Ke;
                    spK_glob.coeffRef(r_ind, c_ind+2) += Ke;
            }

            Ke = (k_new - k_old) * dir_cos[3*node_i+1]*dir_cos[3*node_j+1];
            if( Ke != 0)
                    spK_glob.coeffRef(r_ind+1, c_ind+1) += Ke;

            Ke = (k_new - k_old) * dir_cos[3*node_i+1]*dir_cos[3*node_j+2];
            if( Ke != 0){
                    spK_glob.coeffRef(r_ind+1, c_ind+2) += Ke;
                    spK_glob.coeffRef(r_ind+2, c_ind+1) += Ke;
            }

            Ke = (k_new - k_old) * dir_cos[3*node_i+2]*dir_cos[3*node_j+2];
            if( Ke != 0)
                    spK_glob.coeffRef(r_ind+2, c_ind+2) += Ke;
        }
    }
}
