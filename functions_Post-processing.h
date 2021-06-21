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

#ifndef FUNCTIONS_POST_PROCESSING_H_INCLUDED
#define FUNCTIONS_POST_PROCESSING_H_INCLUDED

#include <vector>
#include <Eigen/Sparse>

using namespace std;
using std::vector;

void write_disp_RF(const string &file_name,
                   const size_t &n_node,
                   const Eigen::VectorXd &u_glob,
                   const Eigen::VectorXd &f_glob);

void WriteNodeDataXML(const string &file_name,
                      const size_t &n_node,
                      const Eigen::VectorXd &u_glob,
                      const Eigen::VectorXd &f_glob);

void write_field_variables(const string &file_name,
                   const size_t &n_ele,
                   const vector<double> &strain,
                   const vector<double> &stress);

void WriteElementDataXML(const string &file_name,
                         const size_t &n_ele,
                         const vector<double> &strain,
                         const vector<double> &stress,
                         const vector<double> &res_Ke,
                         const vector<double> &bound_tensile);

void WriteEndXML(const string &file_name);

#endif // FUNCTIONS_POST_PROCESSING_H_INCLUDED
