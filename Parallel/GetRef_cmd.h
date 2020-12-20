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

#ifndef __GETREF_CMD_H
#define __GETREF_CMD_H


// -- project includes --
#include "Parallel/BroadCast_cmd.h"

/*!
  \class GetNodeRefCommand
  \brief command for getting mesh node reference list

  $Revision$
  $Date$
*/
class GetNodeRefCommand : public BroadcastCommand
{
 public:
  GetNodeRefCommand(const MpiRankAndComm&,const string&);
};

/*!
  \class GetFaceRefCommand
  \brief command for getting mesh node reference list

  $Revision$
  $Date$
*/
class GetFaceRefCommand : public BroadcastCommand
{
 public:
  GetFaceRefCommand(const MpiRankAndComm&,const string&);
};

#endif // __GETREF_CMD_H
