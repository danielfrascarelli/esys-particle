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

#ifndef __VECTOR_EDGE2D_FIELD_SLAVE_H
#define __VECTOR_EDGE2D_FIELD_SLAVE_H

// -- project includes --
#include "Fields/FieldSlave.h"
#include "Model/Mesh2D.h"
#include "tml/comm/comm.h"

// == STL includes --
#include <map>
using std::map;

/*!
  \class VectorEdge2DFieldSlave
  \brief Slave part for saving a vector field defined on the edges
  in a given Mesh2D

  \author Steffen Abe
  $Date$
  $Revision$
*/
class VectorEdge2DFieldSlave : public AFieldSlave
{
 private:
  map<int,Vec3>  m_data;

 protected: 
  Mesh2D *m_mesh;
  Edge2D::VectorFieldFunction m_rdf;
  virtual void SendDataFull();
  virtual void SendDataFullDX();

 public:
  VectorEdge2DFieldSlave(TML_Comm*,Mesh2D*,Edge2D::VectorFieldFunction);
  virtual ~VectorEdge2DFieldSlave()
  {
  }
  
  virtual void sendData();
};

#endif //__VECTOR_EDGE2D_FIELD_SLAVE_H
