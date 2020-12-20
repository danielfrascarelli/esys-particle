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

#ifndef __SCALAR_INTERACTION_FIELD_SLAVE_TAGGED_H
#define __SCALAR_INTERACTION_FIELD_SLAVE_TAGGED_H

// -- project includes --
#include "ScalarInteractionFieldSlave.h"

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
class ScalarInteractionFieldSlaveTagged : public ScalarInteractionFieldSlave<T>
{
 private:
  virtual void SendDataFull();
  virtual void SendDataSum();
  virtual void SendDataMax();

 protected: 
  int m_tag;
  int m_mask;

 public:
  ScalarInteractionFieldSlaveTagged(TML_Comm*,TParallelInteractionStorage<T>*,typename T::ScalarFieldFunction,int,int);
};

#include "ScalarInteractionFieldSlaveTagged.hpp"

#endif //__SCALAR_INTERACTION_FIELD_SLAVE_TAGGED_H
