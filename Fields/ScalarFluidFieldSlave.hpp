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

/*!
  constructor

  \param comm the TML communicator used for sending the data back to the master
  \param ppa a pointer to the particle array
  \param rdf the fluid cell member function to access the data
*/
template <typename T>
ScalarFluidFieldSlave<T>::ScalarFluidFieldSlave(TML_Comm* comm,ParallelParticleArray<T>* ppa,CFluidCell::ScalarFieldFunction rdf):AFieldSlave(comm)
{
  m_ppa=ppa;
  m_rdf=rdf;
}



/*!
  send full field data and position of the fluid cells
*/
template <typename T>
void ScalarFluidFieldSlave<T>::SendDataFull()
{
  vector<pair<Vec3,double> > data_vec;

  data_vec=m_ppa->forAllInnerCellsGet(m_rdf);

  // send data to master
  m_comm->send_gather(data_vec,0);
}

/*!
  send sum only
*/
template <typename T>
void ScalarFluidFieldSlave<T>::SendDataSum()
{
  vector<double> data_vec;

  // get data from particles
  data_vec=m_ppa->forAllInnerCellsGetSum(m_rdf);

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
void ScalarFluidFieldSlave<T>::SendDataMax()
{
  vector<double> data_vec;

  // get data from particles
  data_vec=m_ppa->forAllInnerCellsGetSum(m_rdf);

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
void ScalarFluidFieldSlave<T>::sendData()
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
