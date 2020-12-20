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

#ifndef __VECTORINTERACTIONFIELDMASTER_H
#define __VECTORINTERACTIONFIELDMASTER_H

//--- project includes ---
#include "FieldMaster.h"
#include "Foundation/vec3.h"
#include "Foundation/quintuple.h"

//--- STL includes ---
#include <map>
#include <vector>
#include <iostream>

using std::map;
using std::vector;

/*!
  \class VectorInteractionFieldMaster
  \brief Class for master part of a vector field which is defined on all particles

  \author Steffen Abe
  $Revision$
  $Date$
*/
class VectorInteractionFieldMaster : public AFieldMaster
{
 public:
  typedef std::pair<esys::lsm::quintuple<Vec3,double,Vec3,double,Vec3>, Vec3> IVecData2;
  typedef std::pair<esys::lsm::triplet<int,int,Vec3>, Vec3> DataWithID;
  typedef std::pair<esys::lsm::quintuple<int,int,Vec3,Vec3,Vec3>, Vec3> DataWithPosID;


 protected:
  vector<IVecData2> m_data2; // vector of <pos1,radius1,pos2,radius2,ipos,value> groups
  vector<DataWithID> m_data_with_id;
  vector<DataWithPosID> m_data_with_pos_id;
  vector<pair<Vec3,Vec3> >  m_data; // vector of <position,value> pairs
  vector<Vec3> m_sum_vec;

  virtual void writeAsDX();
  virtual void writeAsSUM();
  virtual void writeAsMAX(){};
  virtual void writeAsRAW_SERIES(){};
  virtual void writeAsRAW2();
  virtual void writeAsRawWithID();
  virtual void writeAsRawWithPosID();

  void collectFull();
  void collectSum();
  void collectMax();
  void collectFull2();
  void collectFullWithID();
  void collectFullWithPosID();

 public:
  VectorInteractionFieldMaster(TML_Comm*,const string&,const string&,const string&,const string&,const string&,int,int,int,bool);
  virtual ~VectorInteractionFieldMaster(){}; 

  virtual void collect();
};

#endif // __VECTORINTERACTIONFIELDMASTER_H
