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

#ifndef __BMESH2DINTERACTION_H
#define __BMESH2DINTERACTION_H

// --- STL includes ---
#include <string>

using std::string;  

// --- project includes ---
#include "Model/BEdge2DInteraction.h"
#include "Model/BCorner2DInteraction.h"
#include "Model/BMesh2DIP.h"

struct BMesh2DInteraction
{
  typedef BMesh2DIP ParameterType;
  typedef BEdge2DInteraction TriIntType;
  typedef BCorner2DInteraction CornerIntType;

  static string getType(){return "BondedMesh2D";};
};

#endif // __BMESH2DINTERACTION_H
