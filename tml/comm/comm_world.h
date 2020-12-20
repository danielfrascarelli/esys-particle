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

#ifndef __COMM_WORLD_H
#define __COMM_WORLD_H

//--- MPI ---
#include <mpi.h>

//--- TML ---
#include "tml/comm/comm.h"

/*!
  \class TML_CommWorld
  \brief abstract base class for communicator

  \author Steffen Abe
  $Revision$
  $Date$
*/
class TML_CommWorld : public TML_Comm
{
public:
  //! constructor
  TML_CommWorld();

};

#endif //__COMM_WORLD_H
