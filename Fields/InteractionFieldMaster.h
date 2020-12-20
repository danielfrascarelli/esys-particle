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

#ifndef __INTERACTIONFIELDMASTER_H
#define __INTERACTIONFIELDMASTER_H

//--- project includes ---
#include "FieldMaster.h"
#include "vec3.h"
#include "Foundation/quintuple.h"
#include "Foundation/triplet.h"

//--- STL includes ---
#include <vector>

using std::vector;

//class TML_Comm;

/*!
  \class ScalarInteractionFieldMaster
  \brief Class for master part of a scalar field which is defined on all particles

  \author Steffen Abe
  $Revision$
  $Date$
*/
class ScalarInteractionFieldMaster : public AFieldMaster
{
 public:
  typedef std::pair<esys::lsm::quintuple<Vec3,double,Vec3,double,Vec3>, double> IVecData2;
  typedef std::pair<esys::lsm::triplet<int,int,Vec3>, double> DataWithID;
  typedef std::pair<esys::lsm::quintuple<int,int,Vec3,Vec3,Vec3>, double> DataWithPosID;

 protected:
  vector<IVecData2> m_data2; // vector of <pos1,radius1,pos2,radius2,ipos,value> groups
  vector<DataWithID> m_data_with_id;
  vector<DataWithPosID> m_data_with_pos_id;
  vector<pair<Vec3,double> >  m_data; // vector of <position,value> pairs
  vector<double> m_sum_vec;

  virtual void writeAsDX();
  virtual void writeAsSUM();
  virtual void writeAsMAX();
  virtual void writeAsRAW_SERIES();
  virtual void writeAsRAW2();
  virtual void writeAsRAW();
  virtual void writeAsRawWithID();
  virtual void writeAsRawWithPosID();

  void collectFull();
  void collectSum();
  void collectMax();
  void collectFull2();
  void collectFullWithID();
  void collectFullWithPosID();

 public:
  ScalarInteractionFieldMaster(TML_Comm*,const string&,const string&,const string&,const string&,const string&,int,int,int,bool);
  ScalarInteractionFieldMaster(TML_Comm*,const string&,const string&,const string&,const string&,const string&,int,int,int,int,int,bool);
  virtual ~ScalarInteractionFieldMaster(){};

  void collect();
 };

#endif //__INTERACTIONFIELDMASTER_H
