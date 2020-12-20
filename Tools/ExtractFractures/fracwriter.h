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

#ifndef __FRACWRITER_H
#define __FRACWRITER_H

// --- project includes ---
#include "fracframe.h"
#include "Foundation/vec3.h"
#include "Geometry/Plane3D.h"

// --- STL includes ---
#include <vector>
#include <string>

using std::vector;
using std::string;

struct fwdata
{
  Vec3 pos;
  Vec3 normal;
  double size;
  double dist;
  int time;
  int id1,id2;
  int ptag1, ptag2;
  int tag;

  fwdata(const FracFrame::fdata&,int);
};

class FracWriter
{
 private:
  vector<fwdata> m_data;
  map<int,int> m_nbrk_map;

  Vec3 m_c1,m_c2,m_c3,m_c4;
  bool with_plane;
  
  void writePlane(const string&);

 public:
  FracWriter();
  void addData(const vector<FracFrame::fdata>&,int);
  void addPlane(const Plane3D&);
  void write(const string&);
  void writeText(const string&);
  void writeProfile(double,double,int,const string&);
  void writeParticleList(const string&);
};

#endif //__FRACWRITER_H
