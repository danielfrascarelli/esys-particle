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
VectorParticleFieldSlave<T>::VectorParticleFieldSlave(TML_Comm* comm,ParallelParticleArray<T>* ppa,typename T::VectorFieldFunction rdf):AFieldSlave(comm)
{
  m_ppa=ppa;
  m_rdf=rdf;
}

/*!
  send data back to master
*/
template <typename T>
void VectorParticleFieldSlave<T>::sendData()
{ 
  vector<pair<int,Vec3> > data_vec;
  vector<pair<int,Vec3> > pos_vec;

  data_vec=m_ppa->forAllInnerParticlesGetIndexed(m_rdf);
  pos_vec=m_ppa->forAllInnerParticlesGetIndexed(typename T::VectorFieldFunction(&T::getPos));
 
  // send data to master
  m_comm->send_gather(data_vec,0);
  m_comm->send_gather(pos_vec,0);
}
