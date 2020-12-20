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
VectorParticleFieldSlaveTagged<T>::VectorParticleFieldSlaveTagged(TML_Comm* comm,ParallelParticleArray<T>* ppa,typename T::VectorFieldFunction rdf,int tag,int mask):VectorParticleFieldSlave<T>(comm,ppa,rdf)
{
  m_tag=tag;
  m_mask=mask;
}

template <typename T>
void VectorParticleFieldSlaveTagged<T>::sendData()
{
  vector<pair<int,Vec3> > data_vec;
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

  // send data to master
  this->m_comm->send_gather(data_vec,0);
  this->m_comm->send_gather(pos_vec,0);
}
