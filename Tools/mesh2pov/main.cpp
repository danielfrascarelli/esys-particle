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

/*
  mesh2pov - tool to convert a surface mesh (Finley format) into
  a Povray3.5 mesh2 object
*/

// --- Project includes ---
#include "Parallel/MeshReader.h"
#include "Model/MeshData.h"
using esys::lsm::MeshReader;

// --- System includes ---
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <utility>

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::atoi;
using std::vector;
using std::ofstream;

int main(int argc,char** argv)
{
  string infilename;
  string outfilename;
  bool options_valid=true;
  
  int args_read=1;
  
  while(args_read<argc){
    string option=string(argv[args_read]);
    if(option=="-i"){
      if(argc>args_read){
        infilename=string(argv[args_read+1]);
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-o"){
      if(argc>args_read){
        outfilename=argv[args_read+1];
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else  {
      cerr << "Unknown option " << option << endl;
      options_valid = false;
      break;
    }
  }

  if(options_valid){
    cout << "converting mesh " << infilename << " to POV object " << outfilename << endl;
    // --- setup reader
    MeshReader meshReader(infilename);
 
    // storage ...
    vector<MeshNodeData> node_vector;
    vector<MeshTriData> tri_vector;

    // --- read data in
    // nodes 
    MeshReader::NodeIterator &niter=meshReader.getNodeIterator();
    cout << "will read " << niter.getNumRemaining() << " nodes" << endl;
    while(niter.hasNext()){
      node_vector.push_back(niter.next());
    }
    cout << "have read " << node_vector.size() << " nodes" << endl;  
    // triangles
    MeshReader::TriIterator &titer=meshReader.getTriIterator();
    cout << "will read " << titer.getNumRemaining() << " triangles" << endl;
    while(titer.hasNext()){
      tri_vector.push_back(titer.next());
    }
    cout << "have read " << tri_vector.size() << " triangles" << endl;  

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
    outfile << "   pigment {rgb <0.3,0.6,0.9>}" << endl;
    outfile << "}" << endl;

    outfile.close();
  }
  
  return 0;
}
