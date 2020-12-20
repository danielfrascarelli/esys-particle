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


#ifndef ESYS_LSMGRAINCOLLECTIONPY_H
#define ESYS_LSMGRAINCOLLECTIONPY_H

#include "Foundation/console.h"

#include "Geometry/GrainCollection.h"
#include "Python/esys/lsm/geometry/GrainPy.h"
#include "Python/esys/lsm/geometry/ParticleCollectionPy.h"

namespace esys
{
  namespace lsm
  {
    class GrainCollectionPy : public GrainCollection<GrainPy>
    {
    public:
      typedef GrainCollection<GrainPy>              Inherited;
      typedef IteratorPy<Inherited::GrainIterator>  GrainIteratorPy;
      typedef Inherited::ParticlePoolPtr            ParticlePoolPtr;
      typedef Inherited::GrainPoolPtr               GrainPoolPtr;

      GrainCollectionPy();

      GrainCollectionPy(ParticlePoolPtr particlePoolPtr);

      GrainCollectionPy(
        ParticlePoolPtr particlePoolPtr,
        GrainPoolPtr grainPoolPtr
      );

      GrainIteratorPy getGrainIteratorPy();

      GrainPy &createGrainPy();

      GrainPy &createGrainWithIdPy(GrainPy::Id id);

    private:
    };

    void exportGrainCollection();
  }
}

#endif
