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

#ifndef __READSNAP_H
#define __READSNAP_H

// --- STL incldes ---
#include <string>
using std::string;

// --- Project includes ---
#include "graph.h"

Graph readSnap(const string&,int,int,int btag=-1);
void readSnap(const string&,int,int,Graph&,int btag=-1);

Graph readGeo(const string&,int,int,int btag=-1);
void readGeo(const string&,int,int,Graph&,int btag=-1);


#endif // __READSNAP_H
