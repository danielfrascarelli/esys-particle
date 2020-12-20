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

#ifndef ESYS_LSMPARTICLEIDPAIRSETPY_H
#define ESYS_LSMPARTICLEIDPAIRSETPY_H

#include <Python/esys/lsm/util/SetPy.h>
#include <Python/esys/lsm/ParticleIdPairPy.h>

namespace esys
{
  namespace lsm
  {
    typedef SetPy<ParticleIdPairPy> ParticleIdPairSetPy;

    void exportParticleIdPairSet();
  }
}

#endif
