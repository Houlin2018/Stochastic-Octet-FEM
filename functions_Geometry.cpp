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

#include "functions_Geometry.h"

#include <iomanip>
#include <fstream>
#include <cmath>

void SetCoords(const size_t &nx, const size_t &ny, const size_t &nz, const double &l_block, vector<vector<double>> &coords){
    for(size_t k = 0; k != nz+1; ++k){
        for(size_t j = 0; j != ny+1; ++j){
            for(size_t i = 0; i != nx+1; ++i){
                if ( (i+j+k)%2==0 )
                    coords.push_back({i * l_block, j * l_block, k * l_block});
            }
        }
    }
}

void SetConnectivity(const double &len, const vector<vector<double>> &coords,
                            vector<vector<vector<size_t> > > &conn,
                            vector<vector<size_t> > &conn_flat){
    /* -------------------------------------------------------
       conn[0]: XY-plane, + direction (vx * vy > 0)
       conn[1]: XY-plane, - direction (vx * vy < 0)
       conn[2]: YZ-plane, + direction
       conn[3]: YZ-plane, - direction
       conn[4]: XZ-plane, + direction
       conn[5]: XZ-plane, - direction
       ------------------------------------------------------- */
    double dis;
    size_t i = 0;
    // conn
    for (auto it = coords.cbegin(); it != coords.cend(); ++it){
        size_t j = i + 1;
        for (auto jt = it + 1; jt != coords.cend(); ++jt){
            dis = pow((*jt)[0]-(*it)[0],2) + pow((*jt)[1]-(*it)[1],2) + pow((*jt)[2]-(*it)[2],2);
            dis = sqrt(dis);
            if ( dis<1.01*len && dis>0.99*len ){
                if ((*jt)[2]-(*it)[2] == 0){
                    if ( ((*jt)[0]-(*it)[0])*((*jt)[1]-(*it)[1]) > 0 )
                        conn[0].push_back({i,j});
                    else conn[1].push_back({i,j});
                }
                else if ((*jt)[0]-(*it)[0] == 0){
                    if ( ((*jt)[1]-(*it)[1])*((*jt)[2]-(*it)[2]) > 0 )
                        conn[2].push_back({i,j});
                    else conn[3].push_back({i,j});
                }
                else{
                    if ( ((*jt)[0]-(*it)[0])*((*jt)[2]-(*it)[2]) > 0 )
                        conn[4].push_back({i,j});
                    else conn[5].push_back({i,j});
                }
            }
            ++j;
        }
        ++i;
    }

    //conn_flat
    for(size_t i = 0; i != conn.size(); ++i){
        conn_flat.insert(conn_flat.end(), conn[i].begin(), conn[i].end() );
    }
}

void WriteMeshVTK(const string &file_name, const size_t &n_node, const size_t &n_ele,
             const vector<vector<double>> &coords,
             const vector<vector<size_t> > &conn_flat){
    ofstream mesh_f;
    mesh_f.open(file_name,ios::out|ios::trunc);
    mesh_f << "# vtk DataFile Version 2.0 \n"
           << "3D_Lattice odb\n"
           << "ASCII\n"
           << "DATASET UNSTRUCTURED_GRID\n"
           << "POINTS " << n_node << " float"
           << endl;
    // print nodal coordinates
    for(auto it = coords.cbegin(); it != coords.cend(); ++it){
        mesh_f << left << setw(4) << setfill(' ') << (*it)[0]
               << left << setw(4) << setfill(' ') << (*it)[1]
               << left << setw(4) << setfill(' ') << (*it)[2]
               << " \t"
               << endl;
    }
    // print connectivity
    mesh_f << "CELLS " << n_ele << ' ' << n_ele * (2+1) << endl;
    for(auto it = conn_flat.cbegin(); it != conn_flat.cend(); ++it){
        mesh_f << "2 \t";
        mesh_f << left << setw(8) << setfill(' ') << (*it)[0]
               << left << setw(8) << setfill(' ') << (*it)[1]
               << endl;
    }
    // print element type
    mesh_f << "CELL_TYPES " << n_ele << endl;
    for(unsigned i = 0; i != n_ele; ++i)
        mesh_f << "3" << endl;

    mesh_f.close();
}


void WriteMeshXML(const string &file_name, const size_t &n_node, const size_t &n_ele,
             const vector<vector<double>> &coords,
             const vector<vector<size_t> > &conn_flat){
    ofstream mesh_f;
    mesh_f.open(file_name,ios::out|ios::trunc);
    mesh_f << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
           << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n"
           << "\t<UnstructuredGrid>\n"
           << "\t\t<Piece NumberOfPoints=\"" << n_node << "\" NumberOfCells=\"" << n_ele << "\">"
           << endl;

    /* print points */
    mesh_f << "\t\t\t<Points>\n"
           << "\t\t\t\t<DataArray type=\"Float32\" Name=\"Nodes\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
    for(size_t i = 0; i != coords.size(); ++i){
        mesh_f << "\t\t\t\t\t" << left << setw(8) << setfill(' ') << coords[i][0]
               << left << setw(8) << setfill(' ') << coords[i][1]
               << left << setw(8) << setfill(' ') << coords[i][2]
               << endl;
    }
    mesh_f << "\t\t\t\t</DataArray>\n"
           << "\t\t\t</Points>" << endl;

    /* print cells */
    /* 1 */
    mesh_f << "\t\t\t<Cells>\n"
           << "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << endl;
    for(auto it = conn_flat.cbegin(); it != conn_flat.cend(); ++it){
        mesh_f << "\t\t\t\t\t" << left << setw(8) << setfill(' ') << (*it)[0]
               << left << setw(8) << setfill(' ') << (*it)[1]
               << endl;
    }
    mesh_f << "\t\t\t\t</DataArray>" << endl;
    /* 2 */
    mesh_f << "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << endl;
    for(size_t i = 0; i != conn_flat.size(); ++i){
        mesh_f << "\t\t\t\t\t" << left << setw(8) << setfill(' ') << 2*(i+1)
               << endl;
    }
    mesh_f << "\t\t\t\t</DataArray>" << endl;
    /* 3 */
    mesh_f << "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
    for(auto it = conn_flat.cbegin(); it != conn_flat.cend(); ++it){
        mesh_f << "\t\t\t\t\t" << left << setw(8) << setfill(' ') << 3
               << endl;
    }
    mesh_f << "\t\t\t\t</DataArray>" << endl;

    mesh_f << "\t\t\t</Cells>" << endl;

    mesh_f.close();
}
