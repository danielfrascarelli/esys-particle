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

#include "WallFieldSlave.h"
#include "Foundation/console.h"

/*!
  Constructor

  \param comm the communicator
  \param w the walls
*/
AWallFieldSlave::AWallFieldSlave(TML_Comm* comm)
  : AFieldSlave(comm)
{}

/*!
  add a wall pointer 

  \param wallp 
*/
void AWallFieldSlave::addWall(CWall* wallp)
{
  console.XDebug() << "AWallFieldSlave::addWall()\n";
  m_wall.push_back(wallp);
}
