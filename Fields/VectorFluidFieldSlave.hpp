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
VectorFluidFieldSlave<T>::VectorFluidFieldSlave(TML_Comm* comm,ParallelParticleArray<T>* ppa,CFluidCell::VectorFieldFunction rdf):AFieldSlave(comm)
{
  m_ppa=ppa;
  m_rdf=rdf;
}

/*!
  send data back to master
*/
template <typename T>
void VectorFluidFieldSlave<T>::sendData()
{ 
  vector<pair<Vec3,Vec3> > data_vec;

  data_vec=m_ppa->forAllInnerCellsGet(m_rdf);
 
  // send data to master
  m_comm->send_gather(data_vec,0);
}


