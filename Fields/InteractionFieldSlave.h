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

#ifndef __INTERACTIONFIELDSLAVE_H
#define __INTERACTIONFIELDSLAVE_H

// -- project includes --
#include "FieldSlave.h"

template <typename I> class TParallelInteractionStorage;
class TML_Comm;

/*!
  \class InteractionFieldSlave
  \brief abstract base class for slave part of scalar field defined on the interactions
 
  \author Steffen Abe
  $Revision$
  $Date$
*/
template <typename T>
class InteractionFieldSlave : public AFieldSlave
{
private:

protected:
  TParallelInteractionStorage<T>* m_pis;
  virtual void SendDataFull()=0;
  virtual void SendDataFull2()=0;
  virtual void SendDataWithID()=0;
  virtual void SendDataWithPosID();
  virtual void SendDataSum()=0;
  virtual void SendDataMax()=0;
	
public:
  InteractionFieldSlave(TML_Comm*,TParallelInteractionStorage<T>*);
  virtual void sendData();
};

#include "InteractionFieldSlave.hpp"

#endif // __INTERACTIONFIELDSLAVE_H

