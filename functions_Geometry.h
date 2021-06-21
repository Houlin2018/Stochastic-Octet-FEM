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

#ifndef FUNCTIONS_GEOMETRY_H_INCLUDED
#define FUNCTIONS_GEOMETRY_H_INCLUDED

//#include <iostream>
#include <vector>
#include <string>


using namespace std;
using std::vector;

void SetCoords(const size_t &nx, const size_t &ny, const size_t &nz, const double &l_block, vector<vector<double>> &coords);

void SetConnectivity(const double &len, const vector<vector<double>> &coords,
                            vector<vector<vector<size_t> > > &conn,
                            vector<vector<size_t> > &conn_flat);
/* -------------------------------------------------------
   conn[0]: XY-plane, + direction (vx * vy > 0)
   conn[1]: XY-plane, - direction (vx * vy < 0)
   conn[2]: YZ-plane, + direction
   conn[3]: YZ-plane, - direction
   conn[4]: XZ-plane, + direction
   conn[5]: XZ-plane, - direction
   ------------------------------------------------------- */

void WriteMeshVTK(const string &file_name, const size_t &n_node, const size_t &n_ele,
             const vector<vector<double>> &coords,
             const vector<vector<size_t> > &conn_flat);

void WriteMeshXML(const string &file_name, const size_t &n_node, const size_t &n_ele,
             const vector<vector<double>> &coords,
             const vector<vector<size_t> > &conn_flat);

#endif // FUNCTIONS_GEOMETRY_H_INCLUDED
