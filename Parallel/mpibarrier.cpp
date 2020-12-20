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

// -- project includes --
#include "mpibarrier.h"

// -- system includes --
#include <console.h>


CMPIBarrier::CMPIBarrier(MPI_Comm comm)
{
  m_comm=comm;
  MPI_Comm_rank(m_comm,&m_id);
}
  
/*!
  Wait on the barrier. The message and the time waited are output to console.XDebug()

  \param msg the message

*/
void CMPIBarrier::wait(const char* msg)
{
  m_time=MPI_Wtime();
  if(m_id==0){
    console.XDebug()<< "Master waiting on Barrier ( " << msg << " )\n";
  } else {
    console.XDebug()<< "Worker " << m_id << " waiting on Barrier ( " << msg << " )\n";
  }
    
  MPI_Barrier(m_comm);
  double p_time=MPI_Wtime();
  if(m_id==0){
    console.XDebug()<< "Master past Barrier ( " << msg << " ) after " << p_time-m_time << " sec \n";
  } else {
    console.XDebug()<< "Worker " << m_id << " past Barrier ( " << msg << " ) after " << p_time-m_time << " sec \n";
  }
}
