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

#ifndef __VECTOR_PARTICLE_FIELD_SLAVE_TAGGED_H
#define __VECTOR_PARTICLE_FIELD_SLAVE_TAGGED__H

// -- project includes --
#include "VectorParticleFieldSlave.h"

class TML_Comm;

template <class T> class ParallelParticleArray;

/*!
  \class VectorParticleFieldSlaveTagged
  \brief class for slave part of scalar field defined on the particles

  \author Steffen Abe
  $Revision$
  $Date$
*/
template <typename T>
class VectorParticleFieldSlaveTagged : public VectorParticleFieldSlave<T>
{
 private:
  int m_tag,m_mask;

 protected: 
 public:
  VectorParticleFieldSlaveTagged(TML_Comm*,ParallelParticleArray<T>*,typename T::VectorFieldFunction,int,int);
  virtual void sendData();
};

#include "VectorParticleFieldSlaveTagged.hpp"

#endif //__SCALAR_PARTICLE_FIELD_SLAVE_TAGGED__H
