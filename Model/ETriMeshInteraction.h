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

#ifndef __ETRIMESHINTERACTION_H
#define __ETRIMESHINTERACTION_H

// --- project includes ---
#include "Model/ETriangleInteraction.h"
#include "Model/EEdgeInteraction.h"
#include "Model/ECornerInteraction.h"
#include "Model/ETriMeshIP.h"

struct ETriMeshInteraction
{
  typedef ETriMeshIP ParameterType;
  typedef ETriangleInteraction TriIntType;
  typedef EEdgeInteraction EdgeIntType;
  typedef ECornerInteraction CornerIntType;
};

#endif //__ETRIMESHINTERACTION_H
