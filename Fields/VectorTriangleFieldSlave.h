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

#ifndef __VECTOR_TRIANGLE_FIELD_SLAVE_H
#define __VECTOR_TRIANGLE_FIELD_SLAVE_H

// -- project includes --
#include "Fields/FieldSlave.h"
#include "Model/TriMesh.h"
#include "tml/comm/comm.h"

// == STL includes --
#include <map>
using std::map;
using std::multimap;

/*!
  \class VectorTriangleFieldSlave
  \brief Slave part for saving a vector field defined on the triangles
  in a given TriMesh

  \author Steffen Abe
  $Date$
  $Revision$
*/
class VectorTriangleFieldSlave : public AFieldSlave
{
 private:
  map<int,Vec3>  m_data;
 protected: 
  TriMesh *m_mesh;
  Triangle::VectorFieldFunction m_rdf;
  virtual void SendDataFull();
  virtual void SendDataFullDX();

 public:
  VectorTriangleFieldSlave(TML_Comm*,TriMesh*,Triangle::VectorFieldFunction);

  virtual ~VectorTriangleFieldSlave()
  {
  }
  
  virtual void sendData();
};

#endif // __VECTOR_TRIANGLE_FIELD_SLAVE_H
