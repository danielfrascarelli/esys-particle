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
#ifndef __EMESH2DINTERACTION_H
#define __EMESH2DINTERACTION_H

// --- project includes ---
#include "Model/EEdge2DInteraction.h"
#include "Model/ECorner2DInteraction.h"
#include "Model/ETriMeshIP.h"

struct EMesh2DInteraction
{
  typedef ETriMeshIP ParameterType;
  typedef EEdge2DInteraction EdgeIntType;
  typedef ECorner2DInteraction CornerIntType;
};

#endif //__EMESH2DINTERACTION_H
