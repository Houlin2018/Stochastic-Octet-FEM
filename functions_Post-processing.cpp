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

#include "functions_Post-processing.h"

#include <iomanip>
#include <fstream>

void write_disp_RF(const string &file_name,
                   const size_t &n_node,
                   const Eigen::VectorXd &u_glob,
                   const Eigen::VectorXd &f_glob){
    ofstream result_f;
    result_f.open(file_name,ios::out|ios::app);
    result_f << '\n' << "POINT_DATA " << n_node << " \n"
           << "VECTORS displacement float" << endl;
    for(size_t i = 0; i != n_node; ++i){
        result_f << left << setw(16) << setfill(' ') << u_glob(3*i)
               << left << setw(16) << setfill(' ') << u_glob(3*i+1)
               << left << setw(16) << setfill(' ') << u_glob(3*i+2)
               << endl;
    }

    result_f << '\n' << "VECTORS reaction_force float" << endl;
    for(size_t i = 0; i != n_node; ++i){
        result_f << left << setw(16) << setfill(' ') << f_glob(3*i)
               << left << setw(16) << setfill(' ') << f_glob(3*i+1)
               << left << setw(16) << setfill(' ') << f_glob(3*i+2)
               << endl;
    }
    result_f.close();
}

void WriteNodeDataXML(const string &file_name,
                      const size_t &n_node,
                      const Eigen::VectorXd &u_glob,
                      const Eigen::VectorXd &f_glob){
    ofstream result_f;
    result_f.open(file_name,ios::out|ios::app);
    result_f << "\t\t\t<PointData>\n"
             << "\t\t\t\t<DataArray type=\"Float32\" Name=\"u_glob\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
    for(size_t i = 0; i != n_node; ++i){
        result_f << "\t\t\t\t\t" << left << setw(16) << setfill(' ') << u_glob(3*i)
                                 << left << setw(16) << setfill(' ') << u_glob(3*i+1)
                                 << left << setw(16) << setfill(' ') << u_glob(3*i+2)
                                 << endl;
    }
    result_f << "\t\t\t\t</DataArray>" << endl;

    result_f << "\t\t\t\t<DataArray type=\"Float32\" Name=\"f_glob\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
    for(size_t i = 0; i != n_node; ++i){
        result_f << "\t\t\t\t\t" << left << setw(16) << setfill(' ') << f_glob(3*i)
                                 << left << setw(16) << setfill(' ') << f_glob(3*i+1)
                                 << left << setw(16) << setfill(' ') << f_glob(3*i+2)
                                 << endl;
    }
    result_f << "\t\t\t\t</DataArray>\n"
             << "\t\t\t</PointData>"
             << endl;
    result_f.close();
}

void write_field_variables(const string &file_name,
                   const size_t &n_ele,
                   const vector<double> &strain,
                   const vector<double> &stress){
    ofstream result_f;
    result_f.open(file_name,ios::out|ios::app);
    result_f << '\n' << "CELL_DATA " << n_ele << " \n"
             << "FIELD stress_and_strain 2" << " \n"
             << "strain_field 1 " << n_ele << " float"
             << endl;
    for(size_t i = 0; i != n_ele; ++i){
        result_f << left << setw(16) << setfill(' ') << strain[i]
                 << endl;
    }

    result_f << " \n" << "stress_field 1 " << n_ele << " float"
             << endl;
    for(size_t i = 0; i != n_ele; ++i){
        result_f << left << setw(16) << setfill(' ') << stress[i]
                 << endl;
    }
    result_f.close();
}

void WriteElementDataXML(const string &file_name,
                         const size_t &n_ele,
                         const vector<double> &strain,
                         const vector<double> &stress,
                         const vector<double> &res_Ke,
                         const vector<double> &bound_tensile){
    ofstream result_f;
    result_f.open(file_name,ios::out|ios::app);
    result_f << "\t\t\t<CellData>\n"
             << "\t\t\t\t<DataArray type=\"Float32\" Name=\"strain\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
    for(size_t i = 0; i != n_ele; ++i){
        result_f << "\t\t\t\t\t" << left << setw(16) << setfill(' ') << strain[i]
                                 << endl;
    }
    result_f << "\t\t\t\t</DataArray>" << endl;

    result_f << "\t\t\t\t<DataArray type=\"Float32\" Name=\"stress\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
    for(size_t i = 0; i != n_ele; ++i){
        result_f << "\t\t\t\t\t" << left << setw(16) << setfill(' ') << stress[i]
                                 << endl;
    }
    result_f << "\t\t\t\t</DataArray>" << endl;

    result_f << "\t\t\t\t<DataArray type=\"Float32\" Name=\"Kr\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
    for(size_t i = 0; i != n_ele; ++i){
        result_f << "\t\t\t\t\t" << left << setw(16) << setfill(' ') << res_Ke[i]
                                 << endl;
    }
    result_f << "\t\t\t\t</DataArray>" << endl;

    result_f << "\t\t\t\t<DataArray type=\"Float32\" Name=\"Strength\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
    for(size_t i = 0; i != n_ele; ++i){
        result_f << "\t\t\t\t\t" << left << setw(16) << setfill(' ') << bound_tensile[i]
                                 << endl;
    }

    result_f << "\t\t\t\t</DataArray>\n"
             << "\t\t\t</CellData>"
             << endl;
    result_f.close();
}


void WriteEndXML(const string &file_name){
    ofstream result_f;
    result_f.open(file_name,ios::out|ios::app);
    result_f << "\t\t</Piece>\n"
           << "\t</UnstructuredGrid>\n"
           << "</VTKFile>" << endl;
    result_f.close();
}
