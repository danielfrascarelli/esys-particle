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


#ifndef MPIINFO_H
#define MPIINFO_H

#include <mpi.h>

/**
 * Convenience class for obtaining rank of current process
 * (and size of a specified MPI_Comm).
 */
class MpiInfo
{
public:
  /**
   * Rank and size for specified MPI_Comm.
   */
  MpiInfo(MPI_Comm comm = MPI_COMM_WORLD);
  
  virtual ~MpiInfo();
  
  /**
   *
   */
  int rank() const;
  
  /**
   *
   */
  int size() const;
  
private:
  MPI_Comm m_comm;
};

#endif
