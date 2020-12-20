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

#ifndef __CHECKED_VECTOR_INTERACTION_FIELD_SLAVE_H
#define __CHECKED_VECTOR_INTERACTION_FIELD_SLAVE_H

// -- project includes --
#include "InteractionFieldSlave.h"

template <typename I> class TParallelInteractionStorage;
class TML_Comm;

/*!
  \class ScalarInteractionFieldSlave
  \brief class for slave part of vector field defined on the particles

  \author Steffen Abe
  $Revision: 450 $
  $Date: 2004-10-04 06:51:59 +0100 (Mon, 04 Oct 2004) $
*/
template <typename T>
class CheckedVectorInteractionFieldSlave : public InteractionFieldSlave<T>
{
 private:
  virtual void SendDataFull();
  virtual void SendDataFull2();
  virtual void SendDataSum();
  virtual void SendDataMax();

 protected:
  typename T::CheckedVectorFieldFunction m_rdf;

 public:
  CheckedVectorInteractionFieldSlave(TML_Comm*,TParallelInteractionStorage<T>*,typename T::CheckedVectorFieldFunction);
};

#include "CheckedVectorInteractionFieldSlave.hpp"

#endif //__CHECKED_SCALAR_INTERACTION_FIELD_SLAVE_H
