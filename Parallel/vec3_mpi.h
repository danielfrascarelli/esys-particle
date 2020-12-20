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

#ifndef __VEC3_MPI_H
#define __VEC3_MPI_H

//--- MPI includes ---
#include <mpi.h>

//--- project includes ---
#include "vec3.h"

struct svec3_mpi{
  MPI_Datatype type;
  bool initialized;
};

static svec3_mpi vec3_mpi={MPI_DATATYPE_NULL,false};

#endif //__VEC3_MPI_H
