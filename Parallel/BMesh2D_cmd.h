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

#ifndef __BMESH2D_CMD_H
#define __BMESH2D_CMD_H

// -- project includes --
#include "Parallel/BroadCast_cmd.h"
#include "Model/BMesh2DIP.h"
#include "Model/BTriMeshIP.h"

/*!
  \class BondedMesh2DIGCommand
  \brief command for adding bonded interactions with 2d mesh

  $Revision$
 */
class BondedMesh2DIGCommand : public BroadcastCommand
{
 public:
  BondedMesh2DIGCommand(const MpiRankAndComm&);
  
  void appendMesh2DParam(const BMesh2DIP&);
  void appendGapBuildPrms(const MeshGapBuildPrms&);
  void appendTagBuildPrms(const MeshTagBuildPrms&);
};
#endif //__BMESH2D_CMD_H
