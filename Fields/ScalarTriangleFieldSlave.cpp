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

#include "Fields/ScalarTriangleFieldSlave.h"
#include "Fields/field_const.h"
#include "console.h"

/*!
  constructor

  \param comm the TML communicator used for sending the data back to the master
  \param mesh the mesh on which the field is defined
  \param rdf the particle member function to access the data
*/

ScalarTriangleFieldSlave::ScalarTriangleFieldSlave(TML_Comm* comm,TriMesh* mesh,Triangle::ScalarFieldFunction rdf)
  :AFieldSlave(comm)
{
  m_mesh=mesh;
  m_rdf=rdf;
}

/*!
  Send data back to master. Recieve the collection type and 
  call the appropriate sending function
*/
void ScalarTriangleFieldSlave::sendData()
{
  int coll_type;
  m_comm->recv_broadcast(coll_type,0);
   
  switch(coll_type){
  case COLL_TYPE_FULL: SendDataFull();break;
  case COLL_TYPE_FULL_DX: SendDataFullDX();break;	

  default: std::cerr << "unknown collection type" << std::endl;
  }
}

/*!
  Send data as id,value pairs
*/
void ScalarTriangleFieldSlave::SendDataFull()
{
  console.XDebug() << "ScalarTriangleFieldSlave::SendDataFull\n";
  vector<pair<int,double> > data_vec;

  // get data
  data_vec=m_mesh->forAllTrianglesGetIndexed(m_rdf);
  
//  for(auto v:data_vec){
//	console.XDebug() << "data: (" << v.first << " , " << v.second << ")\n"; 
//  }
 
  // send data to master
  m_comm->send_gather(data_vec,0);

  console.XDebug() << "end ScalarTriangleFieldSlave::SendDataFull\n";
}

/*!
  send data in a for saving as DX format 
*/
void ScalarTriangleFieldSlave::SendDataFullDX()
{
  console.XDebug() << "ScalarTriangleFieldSlave::SendDataFullDX() - NOT IMPLEMENTED\n";
}
