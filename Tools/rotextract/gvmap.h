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

#ifndef __GVMAP_H
#define __GVMAP_H

//--- project include ---
#include "../../Foundation/vec3.h"

//--- STL includes ---
#include <map>
#include <utility>
#include <set>
#include <vector>
#include <string>

using std::multimap;
using std::pair;
using std::set;
using std::string;
using std::vector;

/*!
  \struct gdata
  \brief structure used for the "interesting" data of each particle
*/
struct gdata
{
  Vec3 vel;
  Vec3 pos;
  double mass;
};

/*!
  \class GVMap
  \brief data structure of grain velocity/rotation claculations
 */
class GVMap
{
 private:
  multimap<int,gdata> m_map;
  set<int> m_tag_set;
  map<int,Vec3> m_vel;
  map<int,Vec3> m_angvel;
  int m_min_tag;

  typedef multimap<int,gdata>::iterator mm_iter;

 public:
  GVMap(int);
  void read(const string&,bool);
  void calc();
};

#endif //__GVMAP_H
