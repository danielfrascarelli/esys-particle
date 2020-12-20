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

#ifndef __SCALAR_INTERACTION_FIELD_SLAVE_H
#define __SCALAR_INTERACTION_FIELD_SLAVE_H

// -- project includes --
#include "Fields/InteractionFieldSlave.h"

template <typename I> class TParallelInteractionStorage;
class TML_Comm;

/*!
  \class ScalarInteractionFieldSlave
  \brief class for slave part of scalar field defined on the interactions

  \author Steffen Abe
  $Revision$
  $Date$
*/
template <typename T>
class ScalarInteractionFieldSlave : public InteractionFieldSlave<T>
{
protected:
  virtual void SendDataFull();
  virtual void SendDataFull2();
  virtual void SendDataWithID();
  virtual void SendDataWithPosID();
  virtual void SendDataSum();
  virtual void SendDataMax();

 protected:
  typename T::ScalarFieldFunction m_rdf;

 public:
  ScalarInteractionFieldSlave(TML_Comm*,TParallelInteractionStorage<T>*,typename T::ScalarFieldFunction);
};

#include "Fields/ScalarInteractionFieldSlave.hpp"

#endif //__SCALAR_INTERACTION_FIELD_SLAVE_H
