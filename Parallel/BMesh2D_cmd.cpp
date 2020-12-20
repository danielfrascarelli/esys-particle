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

#include "Parallel/BMesh2D_cmd.h"
#include "Parallel/sublattice_cmd.h"
 
BondedMesh2DIGCommand::BondedMesh2DIGCommand(const MpiRankAndComm& globalRankAndComm) 
  : BroadcastCommand(globalRankAndComm, CMD_ADDBONDEDMESH2DIG)
{}
  
void BondedMesh2DIGCommand::appendMesh2DParam(const BMesh2DIP &MeshPrms)
{
  append(MeshPrms.getName().c_str());
  append(MeshPrms.getMeshName().c_str());
  append(MeshPrms.k);
  append(MeshPrms.brk);
}
  
void BondedMesh2DIGCommand::appendGapBuildPrms(const MeshGapBuildPrms &buildPrms)
{
  append(buildPrms.getTypeString().c_str());
  append(buildPrms.m_maxGap);
}

void BondedMesh2DIGCommand::appendTagBuildPrms(const MeshTagBuildPrms &buildPrms)
{
  append(buildPrms.getTypeString().c_str());
  append(buildPrms.m_tag);
  append(buildPrms.m_mask);
}
