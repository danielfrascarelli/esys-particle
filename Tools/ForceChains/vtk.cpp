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

#include "vtk.h"

//--- IO includes ---
#include <fstream>
#include <iostream>

using std::ifstream;
using std::ofstream;
using std::cout;
using std::endl;
using std::flush;

//--- STL includes ---
#include <map>
#include <vector>
#include <utility>

using std::map;
using std::vector;
using std::pair;
using std::make_pair;

//--- project includes ---
#include "vec3.h"

struct idata
{
  int p1;
  int p2;
  Vec3 force;
};

/*!
  Convert a file containing interaction forces in RAW2 format into a vtk-xml file

  \param ifname name of the input file
  \param ofname name of the output file
*/
void convert_to_vtk(const string& ifname,const string& ofname)
{
  map<Vec3,int> rev_pos_map;
  map<int,Vec3> pos_map;
  vector<idata> interactions;
 
  Vec3 ppos1; // particle 1 position
  Vec3 ppos2; // particle 2 position
  Vec3 ipos; // interaction position
  Vec3 force; // interaction force
  double r1,r2; // particle radii;
  int pcnt=0; // particle counter

  // open input file
  ifstream infile(ifname.c_str());

  while(!infile.eof()){// until end of input file
    int pnr1,pnr2;
    // read line
    infile >> ppos1 >> r1 >> ppos2 >> r2 >> ipos >> force;
    if(force.norm()>0) { // if force !=(0,0,0) , otherwise ignore interaction
      // try to find particle 1 in map
      map<Vec3,int>::iterator iter=rev_pos_map.find(ppos1);
      if(iter!=rev_pos_map.end()){ // particle already in
	// cout << " found ppos1 : " << ppos1 << " id: " << iter->second << endl;
	pnr1=iter->second;
      } else { // particle not yet in -> insert & update count
	rev_pos_map.insert(make_pair(ppos1,pcnt));
	pos_map.insert(make_pair(pcnt,ppos1));
	pnr1=pcnt;
	pcnt++;
      }
      // same thing for ppos2
      iter=rev_pos_map.find(ppos2);
      if(iter!=rev_pos_map.end()){ // particle already in
	// cout << " found ppos2 : " << ppos2 << " id: " << iter->second << endl;
	pnr2=iter->second;
      } else { // particle not yet in -> insert & update count
	rev_pos_map.insert(make_pair(ppos2,pcnt));
	pos_map.insert(make_pair(pcnt,ppos2));
	pnr2=pcnt;
	pcnt++;
      }
      // add interaction
      idata new_int;
      new_int.p1=pnr1;
      new_int.p2=pnr2;
      new_int.force=force;
      interactions.push_back(new_int);
    }
  }
  // close input file
  infile.close();    

  // open output file
  ofstream outfile(ofname.c_str());
  // write header
  outfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
  outfile << "<UnstructuredGrid>\n";
  outfile << "<Piece NumberOfPoints=\"" << pos_map.size()<< "\" NumberOfCells=\"" << interactions.size() << "\">\n";

  // write point data 
  outfile << "<Points>\n";
  outfile << "<DataArray NumberOfComponents=\"3\" type=\"Float64\" format=\"ascii\">\n";
  for(map<int,Vec3>::iterator iter=pos_map.begin();
      iter!=pos_map.end();
      iter++){
    outfile << iter->second << endl;
  }  
  outfile << "</DataArray>\n";
  outfile << "</Points>\n";

  // write interaction (cell) data
  outfile << "<Cells>\n";
  outfile << "<DataArray type=\"Int32\" NumberOfComponents=\"1\" Name=\"connectivity\" format=\"ascii\">\n";
  for(vector<idata>::iterator iter=interactions.begin();
      iter!=interactions.end();
      iter++){
    outfile << iter->p1 << " " << iter->p2 << endl;
  }
  outfile << "</DataArray>";
  // offsets
  outfile << "<DataArray type=\"Int32\" NumberOfComponents=\"1\" Name=\"offsets\" format=\"ascii\">\n";
  for (size_t i = 1; i < interactions.size()*2; i+=2) outfile << i+1 << "\n";
  outfile << "</DataArray>\n";
  // element type
  outfile << "<DataArray type=\"UInt8\" NumberOfComponents=\"1\" Name=\"types\" format=\"ascii\">\n";
  const int CELL_LINE_TYPE = 3;
  for (size_t i = 0; i < interactions.size(); i++) outfile << CELL_LINE_TYPE << "\n";
  outfile << "</DataArray>\n";  
  outfile << "</Cells>\n";
  // bond data
  outfile << "<CellData>\n";
  // bond strain 
  outfile << "<DataArray type=\"Float64\" Name=\"force\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(vector<idata>::iterator iter=interactions.begin();
      iter!=interactions.end();
      iter++){
    outfile << iter->force.norm() << endl;
  }
  outfile << "</DataArray>\n";  
  outfile << "</CellData>\n";
  // write footer
  outfile << "</Piece>\n";
  outfile << "</UnstructuredGrid>\n";
  outfile << "</VTKFile>\n";
  // close output file
  outfile.close();

}
