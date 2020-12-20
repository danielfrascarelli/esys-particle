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
  \param tag the tag
  \param mask the mask
*/
template <typename T>
ScalarParticleFieldSlaveTagged<T>::ScalarParticleFieldSlaveTagged(TML_Comm* comm,ParallelParticleArray<T>* ppa,typename T::ScalarFieldFunction rdf,int tag,int mask):ScalarParticleFieldSlave<T>(comm,ppa,rdf)
{
  m_tag=tag;
  m_mask=mask;
} 

/*!
  send full field date, position and size of the particles
*/
template <typename T>
void ScalarParticleFieldSlaveTagged<T>::SendDataFull()
{
  vector<pair<int,double> > data_vec;
  vector<pair<int,double> > rad_vec;
  vector<pair<int,Vec3> > pos_vec;

  data_vec =
    this->m_ppa->forAllInnerTaggedParticlesGetIndexed(
      this->m_rdf,
      m_tag,
      m_mask
  );
  pos_vec =
    this->m_ppa->forAllInnerTaggedParticlesGetIndexed(
      typename T::VectorFieldFunction(&T::getPos),
      m_tag,
      m_mask
    );
  rad_vec =
    this->m_ppa->forAllInnerTaggedParticlesGetIndexed(
      typename T::ScalarFieldFunction(&T::getRad),
      m_tag,
      m_mask
    );

  // send data to master
  this->m_comm->send_gather(data_vec,0);
  this->m_comm->send_gather(pos_vec,0);
  this->m_comm->send_gather(rad_vec,0);
}

/*!
  send sum only
*/
template <typename T>
void ScalarParticleFieldSlaveTagged<T>::SendDataSum()
{
  vector<double> data_vec;

  // get data from particles
  this->m_ppa->forAllTaggedInnerParticlesGet(data_vec,this->m_rdf,m_tag,m_mask);

  // sum data
  double sum=0.0;
  for(vector<double>::iterator iter=data_vec.begin();
      iter!=data_vec.end();
      iter++){
    sum+=*iter;
  }

  vector<double> sum_vec;
  sum_vec.push_back(sum);
  this->m_comm->send_gather(sum_vec,0);
}


/*!
  send maximum only
*/
template <typename T>
void ScalarParticleFieldSlaveTagged<T>::SendDataMax()
{
  vector<double> data_vec;

  // get data from particles
  this->m_ppa->forAllTaggedInnerParticlesGet(
    data_vec,
    this->m_rdf,
    m_tag,
    m_mask
  );

  // sum data
  double max=*(data_vec.begin());
  for(vector<double>::iterator iter=data_vec.begin();
      iter!=data_vec.end();
      iter++){
    max=(*iter > max) ? *iter : max;
  }

  vector<double> max_vec;
  max_vec.push_back(max);
  this->m_comm->send_gather(max_vec,0);
}
