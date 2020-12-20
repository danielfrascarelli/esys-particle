/////////////////////////////////////////////////////////////
//                                                         //
// Copyright (c) 2003-2014 by The University of Queensland //
// Centre for Geoscience Computing                         //
// http://earth.uq.edu.au/centre-geoscience-computing      //
//                                                         //
// Primary Business: Brisbane, Queensland, Australia       //
// Licensed under the Open Software License version 3.0    //
// http://www.apache.org/licenses/LICENSE-2.0          //
//                                                         //
/////////////////////////////////////////////////////////////

#ifndef __VECTOR_FLUID_FIELD_SLAVE_H
#define __VECTOR_FLUID_FIELD_SLAVE_H

// -- project includes --
#include "FieldSlave.h"
#include "Model/FluidCell.h"

template <class T> class ParallelParticleArray;
class TML_Comm;

/*!
  \class VectorFluidFieldSlave
  \brief class for slave part of vector field defined on the fluid cells

  \author Qi Shao
  $Revision$
  $Date$
*/
template <typename T>
class VectorFluidFieldSlave : public AFieldSlave
{
 private:

 protected: 
  CFluidCell::VectorFieldFunction m_rdf;
  ParallelParticleArray<T>* m_ppa;

 public:
  VectorFluidFieldSlave(TML_Comm*,ParallelParticleArray<T>*,CFluidCell::VectorFieldFunction);
  virtual void sendData();
};

#include "VectorFluidFieldSlave.hpp"

#endif //__SCALAR_FLUID_FIELD_SLAVE_H
