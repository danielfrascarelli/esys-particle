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

// --- project includes ---
#include "Foundation/console.h"

/*!
  Constructor

  \param comm
  \param wall
  \param rdf the function to read the field from the wall
*/
template <typename WallType>
VectorWallFieldSlave<WallType>::VectorWallFieldSlave(TML_Comm* comm,typename WallType::VectorFieldFunction rdf) 
  : AWallFieldSlave(comm)
{
  console.XDebug() << "VectorWallFieldSlave::VectorWallFieldSlave()\n";
  m_rdf=rdf;
}

/*!
  send data back to master
*/
template <typename WallType>
void VectorWallFieldSlave<WallType>::sendData()
{
  console.XDebug() << "VectorWallFieldSlave::sendData()\n";
  vector<pair<int,Vec3> > data; 
  // get data from wall
  int cnt=0;
  for(typename vector<WallType*>::const_iterator iter=m_wall.begin();
      iter!=m_wall.end();
      iter++){
    data.push_back(make_pair(cnt,((*iter)->*m_rdf)()));
    cnt++;
  }
  // send it to master
  m_comm->send_gather(data,0);
  console.XDebug() << " end VectorWallFieldSlave::sendData()\n";
} 
