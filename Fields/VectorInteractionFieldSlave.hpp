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


#include <Foundation/quadtuple.h>

/*!
  constructor 

  \param comm the TML communicator used for sending the data back to the master
  \param pis a pointer to the interaction storage
  \param rdf the particle member function to access the data
*/
template <typename T>
VectorInteractionFieldSlave<T>::VectorInteractionFieldSlave(TML_Comm* comm,TParallelInteractionStorage<T>* pis,typename T::VectorFieldFunction rdf)
  :InteractionFieldSlave<T>(comm,pis)
{
  m_rdf=rdf;
}

/*!
  send full field data and position of the interaction 
*/
template <typename T>
void VectorInteractionFieldSlave<T>::SendDataFull()
{
  vector<pair<Vec3,Vec3> > data;

  data=this->m_pis->forAllInnerInteractionsGetWithPos(m_rdf);

  // send data to master
  this->m_comm->send_gather(data,0);
}

/*!
  send full field data and position of the interaction 
*/
template <typename T>
void VectorInteractionFieldSlave<T>::SendDataFull2()
{
  vector<pair<esys::lsm::quintuple<Vec3,double,Vec3,double,Vec3>,Vec3> > data;

  data=this->m_pis->forAllInnerInteractionsGetRaw2(m_rdf);

  // send data to master
  this->m_comm->send_gather(data,0);
}

/*!
  send full field data and position of the interaction 
*/
template <typename T>
void VectorInteractionFieldSlave<T>::SendDataWithID()
{
  vector<pair<esys::lsm::triplet<int,int,Vec3>, Vec3> > data;

  // debug output 
  console.XDebug() << "VectorInteractionFieldSlave<T>::SendDataWithID()\n";

  data=this->m_pis->forAllInnerInteractionsGetDataWithID(m_rdf);

  // debug output 
  console.XDebug() << "sending " << data.size() << " data\n"; 

  // send data to master
  this->m_comm->send_gather(data,0);
}

/*!
  send full field data, position of the interaction and pos. and Id of the particles  
*/
template <typename T>
void VectorInteractionFieldSlave<T>::SendDataWithPosID()
{
  vector<pair<esys::lsm::quintuple<int,int,Vec3,Vec3,Vec3>, Vec3> > data;

  // debug output 
  console.XDebug() << "VectorInteractionFieldSlave<T>::SendDataWithPosID()\n";

  data=this->m_pis->forAllInnerInteractionsGetDataWithPosID(m_rdf);

  // debug output 
  console.XDebug() << "sending " << data.size() << " data\n"; 

  // send data to master
  this->m_comm->send_gather(data,0);
}

/*!
  send sum only
*/
template <typename T>
void VectorInteractionFieldSlave<T>::SendDataSum()
{
  vector<Vec3> data_vec;

  // get data from interactions
  this->m_pis->forAllInnerInteractionsGet(data_vec,m_rdf);

  // sum data
  Vec3 sum=Vec3(0.0,0.0,0.0);
  for(vector<Vec3>::iterator iter=data_vec.begin();
      iter!=data_vec.end();
      iter++){
    sum=sum+(*iter);
  }

  vector<Vec3> sum_vec;
  sum_vec.push_back(sum);
  this->m_comm->send_gather(sum_vec,0);
}

/*!
  send max of the field
*/
template <typename T>
void VectorInteractionFieldSlave<T>::SendDataMax()
{}
