/////////////////////////////////////////////////////////////
//                                                         //
// Copyright (c) 2003-2017 by The University of Queensland //
// Centre for Geoscience Computing                         //
// http://earth.uq.edu.au/centre-geoscience-computing      //
//                                                         //
// Primary Business: Brisbane, Queensland, Australia       //
// Licensed under the Open Software License version 3.0    //
// http://www.apache.org/licenses/LICENSE-2.0              //
//                                                         //
/////////////////////////////////////////////////////////////

#include "mesh.h" 

// --- Project includes ---
#include "Parallel/MeshReader.h"
#include "Model/MeshData.h"
using esys::lsm::MeshReader;

// --- System includes ---
#include <string>
#include <fstream>
#include <vector>

using std::vector;
using std::ifstream;
using std::ofstream;
using std::endl;

void do_mesh(const string& infilename,const string& outfilename)
{
  // --- setup reader
  MeshReader meshReader(infilename);
  
  // storage ...
  vector<MeshNodeData> node_vector;
  vector<MeshTriData> tri_vector;

  // --- read data in
  // nodes 
  MeshReader::NodeIterator &niter=meshReader.getNodeIterator();
  while(niter.hasNext()){
    node_vector.push_back(niter.next());
  }
  // triangles
  MeshReader::TriIterator &titer=meshReader.getTriIterator();
  while(titer.hasNext()){
    tri_vector.push_back(titer.next());
  }
  
  // write povray file
  ofstream outfile(outfilename.c_str());
    
  outfile << "mesh2 {" << endl;
  outfile << "   vertex_vectors {" << endl;
  outfile << "      " << node_vector.size() ;
  for(vector<MeshNodeData>::iterator iter=node_vector.begin();
      iter!=node_vector.end();
      iter++){
    outfile << ",\n      <" << iter->x << "," << iter->y << "," << iter->z << ">";
  }
  outfile << "\n   }" << endl;
  outfile << "    face_indices {" << endl;
  outfile << "      " << tri_vector.size();
  for(vector<MeshTriData>::iterator iter=tri_vector.begin();
      iter!=tri_vector.end();
      iter++){
    outfile << ",\n      <" << iter->p1 << "," << iter->p2 << "," << iter->p3 << ">";
  }
  outfile << "\n   }" << endl;
  // texture hardwired - replace with CL option
  outfile << "   pigment {rgbf <0.3,0.3,0.3,0.3>}" << endl;
  outfile << "}" << endl;
  
  outfile.close();
}
