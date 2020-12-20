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

#ifndef __VECTOR_INTERACTION_FIELD_SLAVE_H
#define __VECTOR_INTERACTION_FIELD_SLAVE_H

// -- project includes --
#include "Fields/InteractionFieldSlave.h"

template <typename I> class TParallelInteractionStorage;
class TML_Comm;

/*!
  \class VectorInteractionFieldSlave
  \brief class for slave part of vector field defined on the interactions

  \author Steffen Abe
  $Revision$
  $Date$
*/
template <typename T>
class VectorInteractionFieldSlave : public InteractionFieldSlave<T>
{
 private:
  virtual void SendDataFull();
  virtual void SendDataFull2();
  virtual void SendDataWithID();
  virtual void SendDataWithPosID();
  virtual void SendDataSum();
  virtual void SendDataMax();

 protected:
  typename T::VectorFieldFunction m_rdf;

 public:
  VectorInteractionFieldSlave(TML_Comm*,TParallelInteractionStorage<T>*,typename T::VectorFieldFunction);
};

#include "VectorInteractionFieldSlave.hpp"

#endif // __VECTOR_INTERACTION_FIELD_SLAVE_H
