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

#ifndef __SCALAR_TRIANGLE_FIELD_SLAVE_H
#define __SCALAR_TRIANGLE_FIELD_SLAVE_H

// -- project includes --
#include "Fields/FieldSlave.h"
#include "Model/TriMesh.h"
#include "tml/comm/comm.h"

// == STL includes --
#include <map>
using std::map;
using std::multimap;

/*!
  \class ScalarTriangleFieldSlave
  \brief Slave part for saving a scalar field defined on the triangles
  in a given TriMesh

  \author Steffen Abe
  $Date$
  $Revision$
*/
class ScalarTriangleFieldSlave : public AFieldSlave
{
 private:
  map<int,double>  m_data;

 protected: 
  TriMesh *m_mesh;
  Triangle::ScalarFieldFunction m_rdf;
  virtual void SendDataFull();
  virtual void SendDataFullDX();

 public:
  ScalarTriangleFieldSlave(TML_Comm*,TriMesh*,Triangle::ScalarFieldFunction);

  virtual ~ScalarTriangleFieldSlave()
  {
  }
  
  virtual void sendData();
};

#endif // __SCALAR_TRIANGLE_FIELD_SLAVE_H
