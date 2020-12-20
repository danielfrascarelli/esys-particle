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

#ifndef __SCALARTRIANGLEFIELDMASTER_H
#define __SCALARTRIANGLEFIELDMASTER_H

//--- project includes ---
#include "FieldMaster.h"

//--- STL includes ---
#include <map>
using std::map;

/*!
  \class ScalarTriangleFieldMaster
  \brief Master part of a scalar field which is defined on the triangles in a given triangle mesh

  \author Steffen Abe
  $Revision$
  $Date$
*/
class ScalarTriangleFieldMaster :  public AFieldMaster
{
 protected:
  map<int,double> m_data; // id,value

  virtual void writeAsDX();
  virtual void writeAsRAW();
  virtual void writeAsSUM(){};
  virtual void writeAsMAX(){};
  virtual void writeAsRAW_SERIES(){};

  void collectFull();
  void collectFullDX();
  
 public:
  ScalarTriangleFieldMaster(TML_Comm*,const string&,const string&,const string&,const string&,int,int,int);
  ~ScalarTriangleFieldMaster();
  
  virtual void collect();
  virtual void write();
};

#endif //__SCALARTRIANGLEFIELDMASTER_H
