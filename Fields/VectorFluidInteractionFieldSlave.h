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

#ifndef __VECTOR_FLUID_INTERACTION_FIELD_SLAVE_H
#define __VECTOR_FLUID_INTERACTION_FIELD_SLAVE_H

// -- project includes --
#include "Fields/FieldSlave.h"
#include "Model/FluidInteraction.h"
#include "pis/pi_storage_f.h"

template <typename T> class ParallelInteractionStorage_F;
class TML_Comm;

/*!
  \class VectorFluidInteractionFieldSlave
  \brief class for slave part of vector field defined on the fluid interactions

  \author Qi Shao
  $Revision$
  $Date$
*/

template <typename T>
class VectorFluidInteractionFieldSlave : public AFieldSlave
{
protected:
  virtual void SendDataWithPos();
  virtual void SendDataWithID();
  virtual void SendDataWithIDPos();
  virtual void SendDataWithParticle();
  virtual void SendDataSum();
  virtual void SendDataMax();
  virtual void sendData();

  CFluidInteraction::VectorFieldFunction m_rdf;
  ParallelInteractionStorage_F<T>* m_pis;

 public:
  VectorFluidInteractionFieldSlave(TML_Comm*,ParallelInteractionStorage_F<T>*,CFluidInteraction::VectorFieldFunction);
};

#include "Fields/VectorFluidInteractionFieldSlave.hpp"

#endif //__VECTOR_FLUID_INTERACTION_FIELD_SLAVE_H
