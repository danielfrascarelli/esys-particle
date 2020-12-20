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

#ifndef __SCALAR_PARTICLE_FIELD_SLAVE_H
#define __SCALAR_PARTICLE_FIELD_SLAVE_H

// -- project includes --
#include "FieldSlave.h"

template <class T> class ParallelParticleArray;
class TML_Comm;

/*!
  \class ScalarParticleFieldSlave
  \brief class for slave part of scalar field defined on the particles

  \author Steffen Abe
  $Revision$
  $Date$
*/
template <typename T>
class ScalarParticleFieldSlave : public AFieldSlave
{
 private:
  virtual void SendDataFull();
  virtual void SendDataSum();
  virtual void SendDataMax();

 protected: 
  typename T::ScalarFieldFunction m_rdf;
  ParallelParticleArray<T>* m_ppa;

 public:
  ScalarParticleFieldSlave(TML_Comm*,ParallelParticleArray<T>*,typename T::ScalarFieldFunction);

  virtual void sendData();
};

#include "ScalarParticleFieldSlave.hpp"

#endif //__SCALAR_PARTICLE_FIELD_SLAVE_H
