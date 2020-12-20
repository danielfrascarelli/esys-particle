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

//-- STL includes --
#include <vector>
#include <utility>

using std::vector;
using std::pair;

// -- IO includes --
#include <iostream>

using std::cout;
using std::endl;

/*!
  constructor 

  \param comm the TML communicator used for sending the data back to the master
  \param ppa a pointer to the particle array
  \param rdf the particle member function to access the data
*/
template <typename T>
ScalarParticleFieldSlave<T>::ScalarParticleFieldSlave(TML_Comm* comm,ParallelParticleArray<T>* ppa,typename T::ScalarFieldFunction rdf):AFieldSlave(comm)
{
  m_ppa=ppa;
  m_rdf=rdf;
} 



/*!
  send full field date, position and size of the particles
*/
template <typename T>
void ScalarParticleFieldSlave<T>::SendDataFull()
{
  vector<pair<int,double> > data_vec;
  vector<pair<int,double> > rad_vec;
  vector<pair<int,Vec3> > pos_vec;

  data_vec=m_ppa->forAllInnerParticlesGetIndexed(m_rdf);
  pos_vec=m_ppa->forAllInnerParticlesGetIndexed(typename T::VectorFieldFunction(&T::getPos));
  rad_vec=m_ppa->forAllInnerParticlesGetIndexed(typename T::ScalarFieldFunction(&T::getRad));

  // send data to master
  m_comm->send_gather(data_vec,0);
  m_comm->send_gather(pos_vec,0);
  m_comm->send_gather(rad_vec,0);
}

/*!
  send sum only
*/
template <typename T>
void ScalarParticleFieldSlave<T>::SendDataSum()
{
  vector<double> data_vec;
  
  // get data from particles
  m_ppa->forAllInnerParticlesGet(data_vec,m_rdf);

  // sum data
  double sum=0.0;
  for(vector<double>::iterator iter=data_vec.begin();
      iter!=data_vec.end();
      iter++){
    sum+=*iter;
  }
  
  vector<double> sum_vec;
  sum_vec.push_back(sum);
  m_comm->send_gather(sum_vec,0);
}

/*!
  send maximum only
*/
template <typename T>
void ScalarParticleFieldSlave<T>::SendDataMax()
{
  vector<double> data_vec;
  
  // get data from particles
  m_ppa->forAllInnerParticlesGet(data_vec,m_rdf);

  // sum data
  double max=*(data_vec.begin());
  for(vector<double>::iterator iter=data_vec.begin();
      iter!=data_vec.end();
      iter++){
    max=(*iter > max) ? *iter : max;
  }
  
  vector<double> max_vec;
  max_vec.push_back(max);
  m_comm->send_gather(max_vec,0);
}

/*!
  send data back to master
*/
template <typename T>
void ScalarParticleFieldSlave<T>::sendData()
{
  int coll_type;
  m_comm->recv_broadcast(coll_type,0);

  switch(coll_type){
  case 1: SendDataFull();break;
  case 2: SendDataSum();break;
  case 3: SendDataMax();break;
  default: std::cerr << "unknown collection type" << std::endl;
  }
}
