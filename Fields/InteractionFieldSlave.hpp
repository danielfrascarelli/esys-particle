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

#include "field_const.h"


template <typename T>
InteractionFieldSlave<T>::InteractionFieldSlave(TML_Comm* comm,TParallelInteractionStorage<T>* pis):AFieldSlave(comm)
{
  m_pis=pis;
}

/*!
  Dummy implementation for SendDataWithPosID - does nothing except writing an error message to 
stderr.
 */
template <typename T>
void InteractionFieldSlave<T>::SendDataWithPosID()
{
  std::cerr << "SendDataWithPosID is not implemented for this type of InteractionFieldSaver" << std::endl;
}

/*!
  Send data back to master. Determine the type of data (full/sum) to send back from
  the received coll_type and call the sendDataFull or sendDataSum.
*/
template <typename T>
void InteractionFieldSlave<T>::sendData()
{
  // debug output 
  console.XDebug() << "InteractionFieldSlave<T>::sendData()\n";

  int coll_type;
  m_comm->recv_broadcast(coll_type,0);

  // debug output 
  console.XDebug() << "received coll_type=" << coll_type << "\n"; 

  switch(coll_type){
  case COLL_TYPE_FULL : SendDataFull();break;
  case COLL_TYPE_SUM : SendDataSum();break;
  case COLL_TYPE_MAX : SendDataMax();break;
  case COLL_TYPE_FULL2 : SendDataFull2();break;	
  case COLL_TYPE_FULL_WITH_ID : SendDataWithID();break;	
  case COLL_TYPE_FULL_WITH_POS_ID : SendDataWithPosID();break;	

  default: std::cerr << "unknown collection type" << std::endl;
  }
}
