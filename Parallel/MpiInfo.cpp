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


#include "Parallel/MpiInfo.h"
#include "Parallel/mpibuf.h"

MpiInfo::MpiInfo(MPI_Comm comm) : m_comm(comm)
{
}

MpiInfo::~MpiInfo()
{
}

int MpiInfo::rank() const
{
  int rank;
  MPI_Comm_rank(m_comm, &rank);
  
  return rank;
}

int MpiInfo::size() const
{
  int size;
  MPI_Comm_size(m_comm, &size);
  
  return size;
}
