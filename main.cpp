// This is a finite element solver for the stochastic failure of 3D octect-truss
// system under uniaxial tension. Matrix calculations rely on Eigen linear algebra
// library.
//
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

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <Eigen/Sparse>
#include <algorithm>
#include <fstream>

#include "functions_Geometry.h"
#include "functions_Kernel.h"
#include "functions_Post-processing.h"
#include "functions_Random_vec.h"
#include "functions_Step_control.h"
#include "erf_inv.h"

using namespace std;
using std::vector;
//using namespace Eigen;

//// Debug only --------------//
//#include <chrono>
//using namespace std::chrono;
////--------------------------//


/** ----------------------------------------- Unit System ------------------------------------------- *
  * length:       10^-6 m                                                                             *
  * force:        10^-6 N                                                                             *
  * pressure:     MPa (10^6 Pa)                                                                       *
  * ------------------------------------------------------------------------------------------------- */

/** --------------------------------------- Model Parameters ---------------------------------------- *
  * Young's Modulus:              10^4 MPa                    Link Length:            2 micron        *
  * Mean Tensile Strength:      ~ 10^2 MPa                    Link Section Diameter:  0.5 micron      *
  * ------------------------------------------------------------------------------------------------- */

int main()
{
    /**------------------------------Geometry Definition & Preprocessing-----------------------------**/
    const size_t nx = 8, ny = 4, nz = 4;
    const double l_block = sqrt(2); // 10^-6 m
    const double len = sqrt(2) * l_block;
    /* nodal coordinates */
    vector<vector<double> > coords;
    coords.reserve((nx+3)*(ny+3)*(nz+3)/2);
    SetCoords(nx, ny, nz, l_block, coords);
    /* total number of nodes */
    const size_t n_node = coords.size();

    /* connectivity matrix */
    vector<vector<vector<size_t> > > conn(6);
    /* flattened conn */
    vector<vector<size_t>> conn_flat;
    conn_flat.reserve(6*nx*ny*nz);
    SetConnectivity(len, coords, conn, conn_flat);

//    for(auto it = conn_flat.cbegin(); it != conn_flat.cend(); ++it){
//        cout << (*it)[0] << ", " << (*it)[1] << endl;
//    }

    /* total number of elements */
    const size_t n_ele = conn_flat.size();

    /* write mesh information to VTK_XML file */
//    WriteMeshXML("octet_lattice_CELL.vtu", n_node, n_ele, coords, conn_flat);


    /**-------------------------------------End of Pre-processing-------------------------------------**/
    /**-----------------------------------------------------------------------------------------------**/
    /**--------------------------------------------Kernel---------------------------------------------**/
    /* element stiffness matrix - entry form */
    const double elas = 10000, area = 3.1415926 * 0.5 * 0.5 / 4;

    /* ind_u: list of DOF index of delta_u_0 (prescribed disp.)
     * ind_f: list of DOF index of delta_f_1 (prescribed load) */
    vector<size_t> ind_u, ind_f, ind_RF;
    ind_u.reserve( (nx+1)*(ny+1)+2*(ny+1)*(nz+1)+(nx+1)*(nz+1) );
    ind_f.reserve( 3*n_node - nx*ny+2*ny*nz+nx*nz );
    ind_RF.reserve( (nx+1)*(ny+1) );
    SetDispLoadDomain(nx*l_block, coords, ind_RF, ind_u, ind_f);

    /* global u_0 vector (trial_u0) */
    const size_t size_u = ind_u.size();
    Eigen::VectorXd trial_u0 = Eigen::ArrayXd::Zero(size_u);
    Set_u0(size_u, nx*l_block, coords, ind_u, trial_u0);

    /* global f_1 vector - ALL ZEROS */
    const size_t size_f = ind_f.size();
    Eigen::VectorXd delta_f_1 = Eigen::ArrayXd::Zero(size_f);

    /* permutation vector */
    Eigen::VectorXi perm_Vec(size_u+size_f);
    for (size_t i = 0; i != size_u+size_f; ++i){
        if (i < size_u)
            perm_Vec(i) = ind_u[i];
        else
            perm_Vec(i) = ind_f[i-size_u];
    }

    const size_t n_run = 5;
    const size_t n_batch = 100;
//    const size_t batch_offset = 0*n_run;
    for(size_t run = 0; run != n_run; ++run){
        /* string index of runs for output */
//        string str_run = std::to_string(run+batch_offset);
//        str_run = string(6 - str_run.length(), '0') + str_run;

        /* global stiffness matrix spK*/
        Eigen::SparseMatrix<double> spK_glob(3*n_node,3*n_node); // matrix dimension: dimension of geometry * n_node
        SetStiff_glob(n_ele, elas, area, len, coords, conn_flat, spK_glob);

        /* Variables initialization */
        Eigen::VectorXd u_glob = Eigen::ArrayXd::Zero(size_u+size_f);
        Eigen::VectorXd f_glob = Eigen::ArrayXd::Zero(size_u+size_f);

        vector<double> strain(n_ele), stress(n_ele);
        vector<double> sig_ratio(n_ele);

        // residual element stiffness (damage indicator)
        vector<double> res_Ke(n_ele, elas*area/len);
        // critical elements
        size_t cri_ele_ind_new = n_ele + 10;

        // reaction forces
        double RF = 0;

        // history variables
        size_t n_step = 12*ny*nz;
        Eigen::MatrixXd load_disp = Eigen::MatrixXd::Constant(n_step, 2, 0);

        // random number generation
        vector<double> bound_tensile(n_ele), bound_compressive(n_ele);
//        SetBoundariesWeibull(5, 100, bound_tensile, bound_compressive, 10);
        SetBoundariesGaussWeibull(bound_tensile, bound_compressive, 10);

//        cout << "Bounds = " << endl;
//        for(auto it = bound_tensile.cbegin(); it != bound_tensile.cend(); ++it)
//            cout << *it << endl;

        /** ~~~~~~ Entering step loops ~~~~~~ **/
        for(size_t step = 0; step != n_step; ++step){
//            cout << "------------------------------------" << endl;
//            cout << "STEP = " << step << endl;

            RunElasticSolver(0.01*trial_u0, delta_f_1, ind_u, ind_f, perm_Vec, spK_glob, u_glob, f_glob);

            GetStressStrain_Elastic(elas, area, len, res_Ke, conn_flat, coords, u_glob, strain, stress);

            /* Get loading multiplier */
            GetLoadRatio_Elastic(bound_tensile, bound_compressive, res_Ke, stress, sig_ratio);
//            cout << "RATIO = " << endl;
//            for(auto it = sig_ratio.cbegin(); it != sig_ratio.cend(); ++it)
//                cout << '\t' << *it << endl;

            /* Critical element */
            auto cri_ele_it = std::min_element(sig_ratio.cbegin(), sig_ratio.cend());
            cri_ele_ind_new = distance(sig_ratio.cbegin(), cri_ele_it);
//            cout << "index = " << cri_ele_ind_new << endl;

            /* Correct Field Variables */
            u_glob *= (*cri_ele_it);
            f_glob *= (*cri_ele_it);
            MulVec(strain, *cri_ele_it);
            MulVec(stress, *cri_ele_it);

            /* Update SPK */
            UpdateSPK(cri_ele_ind_new, len, elas*area/len, 0, coords, conn_flat, spK_glob);

            /* Update res_Ke */
            res_Ke[cri_ele_ind_new] = 0;

            /**--------------------------**/
            /**      Record Output       **/
            /**--------------------------**/
            // load-displacement curve
            RF = 0;
            for (auto it = ind_RF.cbegin(); it != ind_RF.cend(); ++it){
                RF += f_glob(*it);
            }
            load_disp(step,0) = u_glob(ind_RF[0]);
            load_disp(step,1) = RF;

            /**--------------------------**/
            /**   Write to Output File   **/
            /**--------------------------**/
//            // Field Variables
//            WriteMeshXML("Brittle_T_"+str_run+"/odb_octet_"+std::to_string(step)+".vtu",
//                         n_node, n_ele, coords, conn_flat);
//            WriteNodeDataXML("Brittle_T_"+str_run+"/odb_octet_"+std::to_string(step)+".vtu",
//                             n_node, u_glob, f_glob);
//            WriteElementDataXML("Brittle_T_"+str_run+"/odb_octet_"+std::to_string(step)+".vtu",
//                                n_ele, strain, stress, res_Ke, bound_tensile);
//            WriteEndXML("Brittle_T_"+str_run+"/odb_octet_"+std::to_string(step)+".vtu");
            /**--------------------------**/
            /* Break Condition */
            if(load_disp(step,1)/load_disp(step,0) < 0.01 * load_disp(0,1)/load_disp(0,0))
                break;

        }/** ~~~~~~ End of step loops ~~~~~~ **/
//        /* Load Displacement Curve */
//        ofstream dbg_f;
//        //dbg_f.open("Brittle_T_"+str_run+"/LoadDisp_Lattice_"+str_run+".dat",ios::out|ios::trunc);
//        dbg_f.open("Brittle_T_results/LoadDisp_octet_"+str_run+".dat",ios::out|ios::trunc);
//        dbg_f << load_disp << endl;
//        dbg_f.close();
        /* maximum loads vector */
        ofstream dbg_f;
        dbg_f.open("Brittle_T_results/PeakLoads_octet_"+std::to_string(n_batch)+".dat",ios::out|ios::app);
        dbg_f << load_disp.colwise().maxCoeff()(1) << endl;
        dbg_f.close();

//        /* first strength vector */
//        dbg_f.open("Brittle_T_results/FirstLoads_octet_"+std::to_string(n_batch)+".dat",ios::out|ios::app);
//        dbg_f << load_disp(0,1) << endl;
//        dbg_f.close();
//
//        /* second strength vector */
//        dbg_f.open("Brittle_T_results/SecondLoads_octet_"+std::to_string(n_batch)+".dat",ios::out|ios::app);
//        dbg_f << load_disp(1,1) << endl;
//        dbg_f.close();
//
//        /* third strength vector */
//        dbg_f.open("Brittle_T_results/ThirdLoads_octet_"+std::to_string(n_batch)+".dat",ios::out|ios::app);
//        dbg_f << load_disp(2,1) << endl;
//        dbg_f.close();
    }
    /* Update new coordinates *///    cout << "minimum_ratio = " << *critical_ele_it << endl;
//    vector<vector<double> > coords_new;
//    get_stress_strain(elas, len, conn_flat, coords, u_glob);

    /**-----------------------------------------End of Kernel------------------------------------------**/
    /**------------------------------------------------------------------------------------------------**/
    /**-------------------------------------------Debugging--------------------------------------------**/

//    cout << "u_1 = " << endl;
//    cout << u_1 << endl;
//    cout << "u_glob = " << endl;
//    cout << u_glob << endl;


//    cout << "Sum of forces = " << temp << endl;
//    ofstream dbg_f;
//    dbg_f.open("dbg_B.dat",ios::out|ios::trunc);
//    dbg_f << '\n' << "Sparse Global Stiffness Matrix:" << '\n'
//          << Eigen::MatrixXd(spK_glob) << endl;
//    dbg_f << load_disp << endl;
//    dbg_f << '\n' << "permutation vector:" << '\n' << perm_Vec << endl;
//    dbg_f.close();

//    for (auto it = ind_u.cbegin(); it != ind_u.cend(); ++it){
//        cout << *it << '\n';
//    }
//    cout << endl;


//    auto start = high_resolution_clock::now();
//    auto stop = high_resolution_clock::now();
//    auto duration = duration_cast<microseconds>(stop - start);
//
//    cout << "Time taken by function: "
//         << duration.count() << " 10^-6s " << endl;


    return 0;
}
