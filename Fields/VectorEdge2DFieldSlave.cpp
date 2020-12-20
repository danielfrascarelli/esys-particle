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

#include <mpi.h>
#include "Fields/VectorEdge2DFieldSlave.h"
#include "Fields/field_const.h"
#include "console.h"

/*!
  constructor

  \param comm the TML communicator used for sending the data back to the master
  \param mesh the mesh on which the field is defined
  \param rdf the particle member function to access the data
*/
VectorEdge2DFieldSlave::VectorEdge2DFieldSlave(TML_Comm* comm,Mesh2D* mesh,Edge2D::VectorFieldFunction rdf)
  :AFieldSlave(comm)
{
  m_mesh=mesh;
  m_rdf=rdf;
}


/*!
  Send data back to master. Recieve the collection type and 
  call the appropriate sending function
*/
void VectorEdge2DFieldSlave::sendData()
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
void VectorEdge2DFieldSlave::SendDataFull()
{
  console.XDebug() << "VectorEdgeFieldSlave::SendDataFull\n";
  vector<pair<int,Vec3> > data_vec;

  // get data
  data_vec=m_mesh->forAllEdgesGetIndexed(m_rdf);
 
  // send data to master
  m_comm->send_gather(data_vec,0);

  console.XDebug() << "end VectorEdgeFieldSlave::SendDataFull\n";
}

/*!
  send data in a for saving as DX format 

  \warning NOT IMPLEMENTED
*/
void VectorEdge2DFieldSlave::SendDataFullDX()
{ 
  console.XDebug() << "VectorEdgeFieldSlave::SendDataFullDX() - NOT IMPLEMENTED\n";
}

