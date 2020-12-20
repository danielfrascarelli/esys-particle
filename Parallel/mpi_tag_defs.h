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

#ifndef __MPI_TAG_DEFS_H
#define __MPI_TAG_DEFS_H

const int SUBLATTICE_INIT_TAG=1; 
const int SUBLATTICE_PARTICLE_TAG=2; 
const int SUBLATTICE_NSP_TAG=3;
const int SUBLATTICE_SHARED_PARTICLE_TAG=4;
const int SUBLATTICE_NBI_TAG=5;
const int SUBLATTICE_BOND_TAG=6;
const int SUBLATTICE_DT_TAG=7;
const int SUBLATTICE_XCHG_TAG=8;
const int SUBLATTICE_NBP_TAG=9;
const int SUBLATTICE_BOUNDARY_PARTICLE_TAG=10;
const int SUBLATTICE_SPACE_TAG=11; 
const int SUBLATTICE_LPARAM_TAG=12;

const int NEIGHBOR_XCHG_TAG=256;
const int PARTICLE_XCHG_TAG=257;
const int PARTICLE_XCHG2_TAG=258;

const int FIELD_DESC_TAG=512;

#endif
