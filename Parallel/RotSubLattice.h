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

#ifndef __ROTSUBLATTICE_H
#define __ROTSUBLATTICE_H

// -- project includes --
#include "Parallel/SubLattice.h"

/*!
  \class TRotSubLattice
  \brief class of a SubLattice of rotational particles

  \author Steffen Abe
  $Revision$
  $Date$
*/
template <typename T>
class TRotSubLattice : public TSubLattice<T>
{
 protected:
  
  // functions doing the actual work adding interaction groups
  virtual bool doAddPIG(const string&,const string&,CVarMPIBuffer&,bool tagged=false);
  virtual bool doAddDamping(const string&,CVarMPIBuffer&);

 public:
  TRotSubLattice(const esys::lsm::CLatticeParam &prm, int rank, MPI_Comm comm, MPI_Comm worker_comm);
  virtual ~TRotSubLattice();
  virtual void setParticleAngularVelocity();
  virtual void addRotBondedIG();
  virtual void addBrittleBeamSCIG();
  virtual void addBrittleBeamDZCIG();
  virtual void addRotThermBondedIG();
};

#include "Parallel/RotSubLattice.hpp"

#endif //__ROTSUBLATTICE_H
