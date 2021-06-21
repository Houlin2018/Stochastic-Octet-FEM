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

#ifndef FUNCTIONS_SOLVER_H_INCLUDED
#define FUNCTIONS_SOLVER_H_INCLUDED

#include <vector>
#include <Eigen/Sparse>

using namespace std;
using std::vector;


void SetStiff_glob(const size_t &n_ele,
                   const double &elas, const double &area, const double &len,
                   const vector<vector<double>> &coords,
                   const vector<vector<size_t> > &conn_flat,
                   Eigen::SparseMatrix<double> &spK_glob);   // output

void SetDispLoadDomain(const double &x_coord,
                       const vector<vector<double>> &coords,
                       vector <size_t> &ind_RF, vector<size_t> &ind_u, vector<size_t> &ind_f);  // outputs

void Set_u0(const size_t &size_u, const double &coord_x,
            const vector<vector<double>> &coords,
            const vector<size_t> &ind_u,
            Eigen::VectorXd &u_0);   // output

void RearrangeStiff(const size_t &size_vec_1, const size_t &size_vec_2,
                    const Eigen::VectorXi &perm_Vec,
                    const Eigen::SparseMatrix<double> &spK_glob,
                    Eigen::SparseMatrix<double> &spK_00,  // output #1
                    Eigen::SparseMatrix<double> &spK_10,  // output #2
                    Eigen::SparseMatrix<double> &spK_11); // output #3

inline bool srt_by_c0 (const vector<double> &vec_i, const vector<double> &vec_j) {
    return (vec_i[0] < vec_j[0]);
}

void Unsort(const size_t &size_u, size_t size_f,
            const vector<size_t> &ind_u, const vector<size_t> &ind_f,
            const Eigen::VectorXd &v_0, const Eigen::VectorXd &v_1,
            Eigen::VectorXd &result);   // output

void RunElasticSolver(const Eigen::VectorXd &u_0, const Eigen::VectorXd &f_1,
                        const vector<size_t> &ind_u, const vector<size_t> &ind_f,
                        const Eigen::VectorXi &perm_Vec,
                        const Eigen::SparseMatrix<double> &spK_glob,
                        Eigen::VectorXd &u_glob, Eigen::VectorXd &f_glob); // outputs

void GetStressStrain_Elastic(const double &elas, const double &area, const double &len,
                             const vector<double> &res_Ke,
                             const vector<vector<size_t> > &conn_flat,
                             const vector<vector<double> > &coords,
                             const Eigen::VectorXd &u_glob,
                             vector<double> &strain, vector<double> &stress); // outputs

#endif // FUNCTIONS_SOLVER_H_INCLUDED
