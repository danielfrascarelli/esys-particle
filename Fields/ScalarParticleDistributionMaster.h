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

#ifndef __SCALARPARTICLEDISTRIBUTIONMASTER_H
#define __SCALARPARTICLEDISTRIBUTIONMASTER_H

//--- project includes ---
#include "ParticleFieldMaster.h"
#include "realdist.h"

class TML_Comm;

/*!
  \class ScalarParticleDistributionMaster
  \brief Class for master part of the distribution/histogram of a  
  scalar field which is defined on particles

  \author Steffen Abe
  $Revision$
  $Date$
*/
class ScalarParticleDistributionMaster : public ScalarParticleFieldMaster
{
 protected:
  RealDist* m_dist;
  int m_dt_write;
  bool m_is_global;
  bool m_is_writing_time;

 public:
  ScalarParticleDistributionMaster(TML_Comm*,const string&,const string&,const string&,int,int,int,int,double,double,int);
  ScalarParticleDistributionMaster(TML_Comm*,const string&,const string&,const string&,int,int,int,int,double,double,int,int,int);
  ~ScalarParticleDistributionMaster();

  virtual bool needSave(int);
  virtual void collect();
  virtual void write();
};

#endif //__SCALARPARTICLEDISTRIBUTIONMASTER_H
