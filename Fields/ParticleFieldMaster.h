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

#ifndef __PARTICLEFIELDMASTER_H
#define __PARTICLEFIELDMASTER_H

//--- project includes ---
#include "FieldMaster.h"
#include "vec3.h"

//--- STL includes ---
#include <map>

using std::map;

class TML_Comm;

/*!
  \class ScalarParticleFieldMaster
  \brief Class for master part of a scalar field which is defined on all particles

  \author Steffen Abe
  $Revision$
  $Date$
*/
class ScalarParticleFieldMaster : public AFieldMaster
{
 protected:
  map<int,double>  m_save_map;
  map<int,double>  m_rad_map;
  map<int,Vec3>  m_pos_map;
  virtual void writeAsDX();
  virtual void writeAsPOV();
  virtual void writeAsSILO();
  virtual void writeAsSUM();
  virtual void writeAsMAX();
  virtual void writeAsRAW_SERIES();
  virtual void writeAsRawWithPosID();

  void collectFull();
  void collectSum();

 public:
  ScalarParticleFieldMaster(TML_Comm*,const string&,const string&,const string&,int,int,int);
  ScalarParticleFieldMaster(TML_Comm*,const string&,const string&,const string&,int,int,int,int,int);
  virtual ~ScalarParticleFieldMaster(){}; 

  virtual void collect();
 };

/*!
  \class VectorParticleFieldMaster
  \brief Class for master part of a vector field which is defined on all particles

  \author Steffen Abe
  $Revision$
  $Date$
*/
class VectorParticleFieldMaster : public AFieldMaster
{
 protected:
  map<int,Vec3>  m_save_map;
  map<int,Vec3>  m_pos_map;
  virtual void writeAsDX();
  virtual void writeAsPOV();
  virtual void writeAsSILO();
  virtual void writeAsSUM();
  virtual void writeAsMAX();
  virtual void writeAsRAW_SERIES();
  virtual void writeAsRAW2();
  virtual void writeAsRawWithID();

 public:
  VectorParticleFieldMaster(TML_Comm*,const string&,const string&,const string&,int,int,int);
  VectorParticleFieldMaster(TML_Comm*,const string&,const string&,const string&,int,int,int,int,int);
  virtual ~VectorParticleFieldMaster(){}; 

  void collect();
};

#endif //__PARTICLEFIELDMASTER_H
