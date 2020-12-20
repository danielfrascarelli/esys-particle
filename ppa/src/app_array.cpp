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

#include "ppa/src/pp_array.h"

/*!
  Constructor, assuming no circular boundary conditions. Just sets up the cartesian communicator 

  \param comm the communicator
  \param dims the dimensions of the process space
*/
AParallelParticleArray::AParallelParticleArray(TML_Comm *comm, const std::vector<unsigned int> &dims):m_comm(comm, dims, vector<bool>(dims.size(),false))
{
}

/*!
  Constructor with given circular boundary conditions. Just sets up the cartesian communicator 

  \param comm the communicator
  \param dims the dimensions of the process space
  \param circ the circular boundary conditions
*/
AParallelParticleArray::AParallelParticleArray(TML_Comm *comm, const std::vector<unsigned int> &dims, const std::vector<bool> &circ):m_comm(comm,dims,circ)
{}
