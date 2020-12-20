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

#include "Parallel/GetRef_cmd.h"
#include "Parallel/sublattice_cmd.h"

/*!
  constructor for GetNodeRefCommand

  \param globalRankAndComm
  \param meshname
*/
GetNodeRefCommand::GetNodeRefCommand(const MpiRankAndComm& globalRankAndComm,const string& meshname)
  : BroadcastCommand(globalRankAndComm, CMD_GETMESHNODEREF)
{
  append(meshname.c_str());
}

/*!
  constructor for GetFaceRefCommand

  \param globalRankAndComm
  \param meshname
*/
GetFaceRefCommand::GetFaceRefCommand(const MpiRankAndComm& globalRankAndComm,const string& meshname)
  : BroadcastCommand(globalRankAndComm, CMD_GETMESHFACEREF)
{
  append(meshname.c_str());
}
