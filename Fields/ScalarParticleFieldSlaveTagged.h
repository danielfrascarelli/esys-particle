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

#ifndef __SCALARPARTICLEFIELDSLAVETAGGED_H
#define __SCALARPARTICLEFIELDSLAVETAGGED_H

// -- project includes --
#include "ScalarParticleFieldSlave.h"

class TML_Comm;

template <class T> class ParallelParticleArray;

/*!
  \class ScalarParticleFieldSlaveTagged
  \brief class for slave part of scalar field defined on tagged particles

  \author Steffen Abe
  $Revision$
  $Date$
*/
template <typename T>
class ScalarParticleFieldSlaveTagged : public ScalarParticleFieldSlave<T>
{
 private:
  int m_tag,m_mask;

  virtual void SendDataFull();
  virtual void SendDataSum();
  virtual void SendDataMax();

 protected: 
 public:
  ScalarParticleFieldSlaveTagged(TML_Comm*,ParallelParticleArray<T>*,typename T::ScalarFieldFunction,int,int);
};

#include "ScalarParticleFieldSlaveTagged.hpp"

#endif //__SCALARPARTICLEFIELDSLAVETAGGED_H
