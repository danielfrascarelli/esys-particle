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

//-- STL includes --
#include <vector>
#include <utility>

using std::vector;
using std::pair;

// -- IO includes --
#include <iostream>

using std::cout;
using std::endl;

#include "Foundation/triplet.h"
#include "pis/pi_storage_f.h"

/*!
  constructor

  \param comm the TML communicator used for sending the data back to the master
  \param pis a pointer to the fluid interaction storage
  \param rdf the fluid interaction member function to access the data
*/
template <typename T>
VectorFluidInteractionFieldSlave<T>::VectorFluidInteractionFieldSlave(TML_Comm* comm,ParallelInteractionStorage_F<T>* pis,CFluidInteraction::VectorFieldFunction rdf):AFieldSlave(comm)
{
  m_rdf=rdf;
  m_pis=pis;
}

/*!
  Send data back to master. Determine the type of data (full/sum) to send back from
  the received coll_type and call the according send function.
*/
template <typename T>
void VectorFluidInteractionFieldSlave<T>::sendData()
{
  // debug output
  console.XDebug() << "InteractionFieldSlave<T>::sendData()\n";

  int coll_type;
  m_comm->recv_broadcast(coll_type,0);

  // debug output
  console.XDebug() << "received coll_type=" << coll_type << "\n";

  switch(coll_type){
  case COLL_TYPE_FULL_WITH_POS : SendDataWithPos();break;
  case COLL_TYPE_FULL_WITH_ID : SendDataWithID();break;
  case COLL_TYPE_FULL_WITH_ID_POS : SendDataWithIDPos();break;
  case COLL_TYPE_FULL_WITH_PARTICLE : SendDataWithParticle();break;
  case COLL_TYPE_SUM : SendDataSum();break;
  case COLL_TYPE_MAX : SendDataMax();break;

  default: std::cerr << "unknown collection type" << std::endl;
  }
}


/*!
  send full field data and position of the interaction
*/
template <typename T>
void VectorFluidInteractionFieldSlave<T>::SendDataWithPos()
{
  vector<pair<Vec3,Vec3> > data;

  data=this->m_pis->forAllInnerInteractionsGetDataWithPos(m_rdf);

  // send data to master
  this->m_comm->send_gather(data,0);
}

/*!
  send full field data and id of the particle
*/
template <typename T>
void VectorFluidInteractionFieldSlave<T>::SendDataWithID()
{
  vector<pair<int,Vec3> > data;

  data=this->m_pis->forAllInnerInteractionsGetDataWithID(m_rdf);

  // send data to master
  this->m_comm->send_gather(data,0);
}


/*!
  send full field data and id, position of the particle
*/
template <typename T>
void VectorFluidInteractionFieldSlave<T>::SendDataWithIDPos()
{
  vector<pair<pair<int,Vec3>,Vec3> > data;

  data=this->m_pis->forAllInnerInteractionsGetDataWithIDPos(m_rdf);

  // send data to master
  this->m_comm->send_gather(data,0);
}

/*!
  send full field data and id, position and radius of the particle
*/
template <typename T>
void VectorFluidInteractionFieldSlave<T>::SendDataWithParticle()
{
  vector<pair<esys::lsm::triplet<int,Vec3,double>,Vec3> > data;

  data=this->m_pis->forAllInnerInteractionsGetDataWithParticle(m_rdf);

  // send data to master
  this->m_comm->send_gather(data,0);
}

/*!
  send sum only
*/
template <typename T>
void VectorFluidInteractionFieldSlave<T>::SendDataSum()
{
  vector<Vec3> data_vec;

  // get data from interactions
  this->m_pis->forAllInnerInteractionsGet(data_vec,m_rdf);

  // sum data
  Vec3 sum=Vec3(0.0,0.0,0.0);
  for(vector<Vec3>::iterator iter=data_vec.begin();
      iter!=data_vec.end();
      iter++){
    sum+=*iter;
  }

  vector<Vec3> sum_vec;
  sum_vec.push_back(sum);
  this->m_comm->send_gather(sum_vec,0);
}


/*!
  send maximum only
*/
template <typename T>
void VectorFluidInteractionFieldSlave<T>::SendDataMax()
{
}
