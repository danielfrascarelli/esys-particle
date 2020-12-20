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

#include "ASubLattice.h"
#include "console.h"
#include "mpibarrier.h"

ASubLattice::~ASubLattice()
{}

/*!
  Adjust the size of the neighbor table

  \param size the size
*/
void ASubLattice::setNTSize(int size)
{
  //m_ntable->resize(size);
}

/*!
  Default implementation for complex init. Just prints an error message.
*/
void ASubLattice::initComplex()
{
  CMPIBarrier barrier(MPI_COMM_WORLD);
  CVarMPIBuffer buffer(MPI_COMM_WORLD);
  console.XDebug()<<"ASubLattice::initComplex()\n";

  console.Critical() << "Complex init not implemented for this type of lattice!\n";
  buffer.receiveBroadcast(0);
  barrier.wait("ASubLattice::initComplex()");
}
