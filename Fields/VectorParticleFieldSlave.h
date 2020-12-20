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

#ifndef __VECTOR_PARTICLE_FIELD_SLAVE_H
#define __VECTOR_PARTICLE_FIELD_SLAVE_H

// -- project includes --
#include "FieldSlave.h"

template <class T> class ParallelParticleArray;
class TML_Comm;

/*!
  \class VectorParticleFieldSlave
  \brief class for slave part of scalar field defined on the particles

  \author Steffen Abe
  $Revision$
  $Date$
*/
template <typename T>
class VectorParticleFieldSlave : public AFieldSlave
{
 private:

 protected: 
  typename T::VectorFieldFunction m_rdf;
  ParallelParticleArray<T>* m_ppa;

 public:
  VectorParticleFieldSlave(TML_Comm*,ParallelParticleArray<T>*,typename T::VectorFieldFunction);
  virtual void sendData();
};

#include "VectorParticleFieldSlave.hpp"

#endif //__SCALAR_PARTICLE_FIELD_SLAVE_H
