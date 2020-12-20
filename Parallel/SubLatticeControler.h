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

#ifndef __SUBLATTICECONTROLER_H
#define __SUBLATTICECONTROLER_H

//--- MPI --- 
#include <mpi.h>

//--- TML includes ---
#include "tml/comm/comm.h"
#include "tml/comm/comm_world.h"

#include "Parallel/SubLattice.h"
#include "Parallel/CheckPointable.h"
#include "Foundation/Timer.h"

#include <boost/shared_ptr.hpp>

class CheckPointer;

/*!
  \class CSubLatticeControler
  \brief class for control of a SubLattice

  Does initialisation and control of a TSubLattice and comunicates with the TLatticeMaster
 
  \author Steffen Abe
  $Revision$
  $Date$
*/
class CSubLatticeControler : esys::lsm::CheckPointable
{
 private:
  int m_global_rank;
  int m_global_size;
  int m_local_rank;
  int m_local_size;

  MPI_Comm m_global_comm; // global MPI communicator
  MPI_Comm m_local_comm; //  MPI communicator of the spawned workers
  MPI_Group m_global_group, m_local_group; // MPI groups 

  TML_Comm m_tml_global_comm;
  TML_Comm m_tml_local_comm;
  
  ASubLattice*  m_lattice;
  CheckPointer* m_pCheckPointer; // for restatable checkpoints
  CheckPointer* m_pSnapShooter; // for visualization dumps
  std::string   m_timingFileName;
  typedef boost::shared_ptr<MpiWTimers> MpiWTimersPtr;
  MpiWTimersPtr m_timersPtr;

 public:
  CSubLatticeControler();
  ~CSubLatticeControler();

  void initMPI();
  void makeLattice();
  void initLattice();
  void initLatticeCirc();
  void init2DTriangularLocal();  
  void init3DTriangularLocal();
  void searchNeighbors();
  void performTiming();
  void saveTimingData();
  void getIdParticleData();
  void setTimeStepSize();
  void setTimingFileName(const std::string &timingFileName) {m_timingFileName=timingFileName;}
  const std::string &getTimingFileName() const {return m_timingFileName;}
  void do2dCalculations();
  void getNumParticles();
  void findParticleNearestToPoint();
  void getParticlePosn();
  void moveSingleParticle();
//  void getBondGroupIdPairs();
  void translateMeshBy();
  void rotateMeshBy();
  void run();

  void setVerbosity();
  void initializeConsole();
  void setConsoleFilename();
  void setConsoleBuffered();
  
  virtual void saveCheckPointData(std::ostream &oStream);
  virtual void saveSnapShotData(std::ostream &oStream);

  virtual void loadCheckPointData(std::istream &iStream);
  
};

#endif //__SUBLATTICECONTROLER_H
