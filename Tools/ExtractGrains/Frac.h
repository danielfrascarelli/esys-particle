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

#ifndef __FRAC_H
#define __FRAC_H

#include "graph.h"
#include <map>
#include <set>
#include <iostream>
#include <string>

using std::map;
using std::multimap;
using std::set;
using std::ostream;
using std::string;

class Frac
{
 private:
  multimap<int,int> m_old_grain_map;
  multimap<int,int> m_new_grain_to_particle_map;
  map<int,int> m_new_grain_map;
  map<int,set<int> > m_grain_to_grain_map;
  map<int,int> m_tag_map;
  map<int,double> m_old_grain_mass_map;
  map<int,double> m_new_grain_mass_map;

  map<int,Vec3> get_move_vectors(const Graph&);
  map<int,double> get_grain_mass(const Graph&);

 public:
  Frac(const Graph&, const Graph&);

  void writeMassRatio(ostream&,double);
  int writeAllMass(ostream&,double,int,bool with_tag=false);
  void writeAsVtk(const string&,float,const Graph&,bool);
  pair<double,double> getSplitAbrasion(Graph&, double,double);
};
#endif //__FRAC_H
