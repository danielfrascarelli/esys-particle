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


#ifndef CHECKPOINTER_H
#define CHECKPOINTER_H

// --- MPI includes ---
#include <mpi.h>

//--- TML includes ---
#include "tml/comm/comm.h"

#include <string>
#include <iostream>


namespace esys
{
  namespace lsm
  {
    class CheckPointable;
  }
}

/**
 * Saves the state of a model.
*/
class CheckPointer
{
public:
  /**
  */
  CheckPointer(esys::lsm::CheckPointable &checkPointable, MPI_Comm mpiComm=MPI_COMM_WORLD);

  /**
  */
  virtual ~CheckPointer();

  /**
   * Saves the state of a model.
   */
  virtual void saveRestartable();
  virtual void saveDump();
  virtual void saveThroughMaster(TML_Comm&);
  virtual void loadCheckPoint();

  /**
   * Saves the state of a model to specified stream.
   */
  virtual void save(std::ostream &oStream);

  MPI_Comm getMpiComm() const;
  
  void setMpiComm(MPI_Comm mpiComm);

private:
  esys::lsm::CheckPointable *m_pCheckPointable;
  MPI_Comm m_mpiComm;
};

#endif
