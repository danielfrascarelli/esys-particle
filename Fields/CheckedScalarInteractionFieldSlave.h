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

#ifndef __CHECKED_SCALAR_INTERACTION_FIELD_SLAVE_H
#define __CHECKED_SCALAR_INTERACTION_FIELD_SLAVE_H

// -- project includes --
#include "InteractionFieldSlave.h"

template <typename I> class TParallelInteractionStorage;
class TML_Comm;

/*!
  \class ScalarInteractionFieldSlave
  \brief class for slave part of scalar field defined on the particles

  \author Steffen Abe
  $Revision$
  $Date$
*/
template <typename T>
class CheckedScalarInteractionFieldSlave : public InteractionFieldSlave<T>
{
 private:
  virtual void SendDataFull();
  virtual void SendDataFull2();
  virtual void SendDataWithID();
  virtual void SendDataSum();
  virtual void SendDataMax();

 protected:
  typename T::CheckedScalarFieldFunction m_rdf;

 public:
  CheckedScalarInteractionFieldSlave(TML_Comm*,TParallelInteractionStorage<T>*,typename T::CheckedScalarFieldFunction);
};

#include "CheckedScalarInteractionFieldSlave.hpp"

#endif //__CHECKED_SCALAR_INTERACTION_FIELD_SLAVE_H
