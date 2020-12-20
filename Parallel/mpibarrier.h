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

#ifndef __MPIBARRIER_H
#define __MPIBARRIER_H

// -- system includes --
#include <mpi.h>

// -- project includes --

/*!
  \class CMPIBarrier
  \brief A convenience class encapsulating an MPI barrier. Includes timing of the wait and a debug message ( via console.XDebug() )

  \author Steffen Abe
  $Revision$
  $Date$

*/
class CMPIBarrier
{
private:
  MPI_Comm m_comm; //!< the MPI Communicator used 
  int m_id;         //!< the rank of the process
  double m_time; 


public:
  CMPIBarrier(MPI_Comm comm=MPI_COMM_WORLD);

  void wait(const char*);
};

#endif // __MPIBARRIER_H
