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

#ifndef __PFCC_H
#define __PFCC_H

// --- project includes ---
#include "../../Foundation/vec3.h"

//--- STL includes ---
#include <utility>
using std::pair;

void pfc_convert(const string&,const string&,const Vec3&,const Vec3&,int,int,int,double);
pair<Vec3,Vec3> read_bbx(const string&);

#endif //__PFCC_H
