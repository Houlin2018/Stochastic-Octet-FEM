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

#include "functions_Kernel.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

void SetStiff_glob(const size_t &n_ele,
                          const double &elas, const double &area, const double &len,
                          const vector<vector<double> > &coords,
                          const vector<vector<size_t> > &conn_flat,
                          Eigen::SparseMatrix<double> &spK_glob){
    vector<Eigen::Triplet<double> > tri_list;
    tri_list.reserve(6*n_ele);
    size_t r_ind, c_ind;
    double cx, cy, cz;
    double Ke; // ENTRY value of element stiffness matrix
    for(auto it = conn_flat.cbegin(); it != conn_flat.cend(); ++it){ // it: element loop
        /* calculate direction cosine of an element */
        cx = (coords[ (*it)[1] ][0] - coords[ (*it)[0] ][0]) / len;
        cy = (coords[ (*it)[1] ][1] - coords[ (*it)[0] ][1]) / len;
        cz = (coords[ (*it)[1] ][2] - coords[ (*it)[0] ][2]) / len;
        vector<double> dir_cos{cx, cy, cz, -cx, -cy, -cz};
        /* assembly */
        for (size_t node_i = 0; node_i != 2; ++node_i){ // node_i, node_j: node loops
            for (size_t node_j = 0; node_j != 2; ++node_j){
                /* row/column index of the spK_glob */
                r_ind = 3 * (*it)[node_i];
                c_ind = 3 * (*it)[node_j];
                /* construct spK_glob */
                Ke = elas * area * dir_cos[3*node_i]*dir_cos[3*node_j] / len;
                if( Ke != 0)
                        tri_list.push_back(Eigen::Triplet<double> (r_ind, c_ind, Ke) );

                Ke = elas * area * dir_cos[3*node_i+1]*dir_cos[3*node_j] / len;
                if( Ke != 0){
                        tri_list.push_back(Eigen::Triplet<double> (r_ind+1,c_ind,Ke) );
                        tri_list.push_back(Eigen::Triplet<double> (r_ind,c_ind+1,Ke) );
                }

                Ke = elas * area * dir_cos[3*node_i+2]*dir_cos[3*node_j] / len;
                if( Ke != 0){
                        tri_list.push_back(Eigen::Triplet<double> (r_ind+2,c_ind,Ke) );
                        tri_list.push_back(Eigen::Triplet<double> (r_ind,c_ind+2,Ke) );
                }

                Ke = elas * area * dir_cos[3*node_i+1]*dir_cos[3*node_j+1] / len;
                if( Ke != 0)
                        tri_list.push_back(Eigen::Triplet<double> (r_ind+1,c_ind+1,Ke) );

                Ke = elas * area * dir_cos[3*node_i+1]*dir_cos[3*node_j+2] / len;
                if( Ke != 0){
                        tri_list.push_back(Eigen::Triplet<double> (r_ind+1,c_ind+2,Ke) );
                        tri_list.push_back(Eigen::Triplet<double> (r_ind+2,c_ind+1,Ke) );
                }

                Ke = elas * area * dir_cos[3*node_i+2]*dir_cos[3*node_j+2] / len;
                if( Ke != 0)
                        tri_list.push_back(Eigen::Triplet<double> (r_ind+2,c_ind+2,Ke) );

            }
        }
    }
    spK_glob.setFromTriplets(tri_list.begin(), tri_list.end());

    //spK_glob.makeCompressed();
}

void SetDispLoadDomain(const double &x_coord, // x_coord of face with prescribed nonzero disp.
                       const vector<vector<double>> &coords,
                       vector <size_t> &ind_RF,
                       vector<size_t> &ind_u,
                       vector<size_t> &ind_f){
    size_t node_i = 0;
    for (auto it = coords.cbegin(); it != coords.cend(); ++it){
        // DOF with 0 and Nonzero displacements:
        if ( (*it)[0] == 0){ // ux = 0
            ind_u.push_back(3*node_i);
        }
        else if ( (*it)[0] == x_coord){ // ux = x_coord
            ind_u.push_back(3*node_i);
            ind_RF.push_back(3*node_i);
        }
        else {
            ind_f.push_back(3*node_i);
        }

        // DOF with 0 displacements:
        if ( (*it)[1] == 0 ){ // uy = 0
            ind_u.push_back(3*node_i+1);
        }
        else{
            ind_f.push_back(3*node_i+1);
        }

        if ( (*it)[2] == 0 ) {// uz = 0
            ind_u.push_back(3*node_i+2);
        }
        else{
            ind_f.push_back(3*node_i+2);
        }
        ++node_i;
    }
}

void Set_u0(const size_t &size_u, const double &coord_x,
            const vector<vector<double>> &coords,
            const vector<size_t> &ind_u,
            Eigen::VectorXd &u_0){
    for (size_t i = 0; i != size_u; ++i){
        if ( ind_u[i] % 3 == 0 ){
            if (coords[ind_u[i]/3 ][0] == coord_x){
                u_0[i] = 1;
//                cout << "dbg_DOF = " << ind_u[i] << endl;
//                cout << "dbg_Node: " << ind_u[i]/3 << ",\t" << "coord_X: " << coords[ind_u[i]/3][0] << endl;
            }
        }
    }
}

/* rearrange spK_glob */
void RearrangeStiff(const size_t &size_vec_1, const size_t &size_vec_2,
                     const Eigen::VectorXi &perm_Vec,
                     const Eigen::SparseMatrix<double> &spK_glob,
                     Eigen::SparseMatrix<double> &spK_00,  // output #1
                     Eigen::SparseMatrix<double> &spK_10,  // output #2
                     Eigen::SparseMatrix<double> &spK_11){ // output #3
    // construct permutation matrix sp_Perm: switching columns of spK; sp_Perm^T ... rows of spK
    vector<Eigen::Triplet<double> > tri_list;
    tri_list.reserve(size_vec_1+size_vec_2); // size_vec_1+size_vec_2 = n_DOF
    Eigen::SparseMatrix<double> sp_Perm(size_vec_1+size_vec_2,size_vec_1+size_vec_2);
    for(size_t i = 0; i != size_vec_1+size_vec_2; ++i)
        tri_list.push_back(Eigen::Triplet<double> (perm_Vec(i),i,1) );
    sp_Perm.setFromTriplets(tri_list.begin(), tri_list.end());

    // rearrange K_glob
    Eigen::SparseMatrix<double> spK_sorted(size_vec_1+size_vec_2,size_vec_1+size_vec_2);
    spK_sorted = sp_Perm.transpose() * spK_glob * sp_Perm;

    // construct output block matrices
    spK_00 = spK_sorted.topLeftCorner(size_vec_1,size_vec_1);
    spK_10 = spK_sorted.bottomLeftCorner(size_vec_2,size_vec_1);
    spK_11 = spK_sorted.bottomRightCorner(size_vec_2,size_vec_2);
}

void Unsort(const size_t &size_u, size_t size_f,
            const vector<size_t> &ind_u, const vector<size_t> &ind_f,
            const Eigen::VectorXd &v_0, const Eigen::VectorXd &v_1,
            Eigen::VectorXd &result){
    vector<vector<double> > srt(size_f+size_u);
    for (size_t i = 0; i != size_u; ++i)
        srt[i] = {(double) ind_u[i], v_0(i)};
    for (size_t i = 0; i != size_f; ++i)
        srt[size_u+i] = {(double) ind_f[i], v_1(i)};

    sort(srt.begin(), srt.end(), srt_by_c0);

    for(size_t i = 0; i != size_u+size_f; ++i)
        result(i) = srt[i][1];
}


void RunElasticSolver(const Eigen::VectorXd &u_0,
                      const Eigen::VectorXd &f_1,
                      const vector<size_t> &ind_u,
                      const vector<size_t> &ind_f,
                      const Eigen::VectorXi &perm_Vec,
                      const Eigen::SparseMatrix<double> &spK_glob,
                      Eigen::VectorXd &u_glob,
                      Eigen::VectorXd &f_glob){
    size_t size_u = ind_u.size();
    size_t size_f = ind_f.size();
    /* rearrange spK_glob */
    Eigen::SparseMatrix<double> spK_00(size_u,size_u), spK_10(size_f,size_u), spK_11(size_f,size_f);
    RearrangeStiff(size_u, size_f, perm_Vec, spK_glob, spK_00, spK_10, spK_11);

    /* solve for u_1 */
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower> cg;
    cg.compute(spK_11);
    Eigen::VectorXd u_1 = cg.solve( - spK_10 * u_0);
    /* f_0 */
    Eigen::VectorXd f_0 = spK_00 * u_0 + spK_10.transpose() * u_1;

    /* construct u_glob and f_glob by ind_u and ind_f (inverse sorting) */
    Unsort(size_u, size_f, ind_u, ind_f, u_0, u_1, u_glob);
    Unsort(size_u, size_f, ind_u, ind_f, f_0, f_1, f_glob);
}

void GetStressStrain_Elastic(const double &elas, const double &area, const double &len,
                             const vector<double> &res_Ke,
                             const vector<vector<size_t> > &conn_flat,
                             const vector<vector<double> > &coords,
                             const Eigen::VectorXd &u_glob,
                             vector<double> &strain, vector<double> &stress){
    size_t node_i, node_j;
    double len_new;
    for(size_t i = 0; i != conn_flat.size(); ++i){
        node_i = conn_flat[i][0];
        node_j = conn_flat[i][1];
        // len_original is the same for all
        len_new = pow(coords[node_j][0]+u_glob[3*node_j  ] - coords[node_i][0]-u_glob[3*node_i  ], 2 )
                + pow(coords[node_j][1]+u_glob[3*node_j+1] - coords[node_i][1]-u_glob[3*node_i+1], 2 )
                + pow(coords[node_j][2]+u_glob[3*node_j+2] - coords[node_i][2]-u_glob[3*node_i+2], 2 );

        len_new = pow(len_new,0.5);
        strain[i] = ( len_new - len ) / len;

        if (res_Ke[i] == 0)
            stress[i] = 0;
        else
            stress[i] = elas * strain[i];
    }
}
