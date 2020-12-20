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


#include "Geometry/RegularBlockGenerator.h"
#include "Geometry/GridIterator.h"
#include "Geometry/SimpleParticle.h"

namespace esys
{
  namespace lsm
  {
    RegularBlockGenerator::RegularBlockGenerator(
      NTable            &nTable,
      ParticlePool      &particlePool,
      const BoundingBox &bBox,
      const BoolVector  &periodicDimensions,
      double            tolerance,
      double            sphereRadius
    ) : BlockGenerator(nTable, particlePool, bBox, periodicDimensions, tolerance),
        m_radius(sphereRadius)
    {
    }

    RegularBlockGenerator::~RegularBlockGenerator()
    {
    }
    
    double RegularBlockGenerator::getRadius() const
    {
      return m_radius;
    }

    double RegularBlockGenerator::getGridRadius() const
    {
      return m_radius;
    }

    void RegularBlockGenerator::generate()
    {
      generateSeedParticles();
    }
  }
}
