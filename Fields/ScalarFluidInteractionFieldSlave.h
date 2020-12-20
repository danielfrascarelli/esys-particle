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

#ifndef __SCALAR_FLUID_INTERACTION_FIELD_SLAVE_H
#define __SCALAR_FLUID_INTERACTION_FIELD_SLAVE_H

// -- project includes --
#include "Fields/FieldSlave.h"
#include "Model/FluidInteraction.h"
#include "pis/pi_storage_f.h"

//template <typename I> class TParallelInteractionStorage;
template <typename T> class ParallelInteractionStorage_F;
class TML_Comm;

/*!
  \class ScalarFluidInteractionFieldSlave
  \brief class for slave part of scalar field defined on the fluid interactions

  \author Qi Shao
  $Revision$
  $Date$
*/

template <typename T>
class ScalarFluidInteractionFieldSlave : public AFieldSlave
{
protected:
  virtual void SendDataWithPos();
  virtual void SendDataWithID();
  virtual void SendDataWithIDPos();
  virtual void SendDataWithParticle();
  virtual void SendDataSum();
  virtual void SendDataMax();
  virtual void sendData();

  CFluidInteraction::ScalarFieldFunction m_rdf;
  ParallelInteractionStorage_F<T>* m_pis;

 public:
  ScalarFluidInteractionFieldSlave(TML_Comm*,ParallelInteractionStorage_F<T>*,CFluidInteraction::ScalarFieldFunction);
};

#include "Fields/ScalarFluidInteractionFieldSlave.hpp"

#endif //__SCALAR_FLUID_INTERACTION_FIELD_SLAVE_H
