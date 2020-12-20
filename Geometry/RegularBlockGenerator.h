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


#ifndef ESYS_LSMREGULARBLOCKGENERATOR_H
#define ESYS_LSMREGULARBLOCKGENERATOR_H

#include <Geometry/BlockGenerator.h>

namespace esys
{
  namespace lsm
  {
    typedef std::vector<bool> BoolVector;
    /**
     *
     */
    class RegularBlockGenerator : public BlockGenerator
    {
    public:
      RegularBlockGenerator(
        NTable            &nTable,
        ParticlePool      &particlePool,
        const BoundingBox &bBox,
        const BoolVector  &periodicDimensions,
        double            tolerance,
        double            sphereRadius
      );

      virtual ~RegularBlockGenerator();

      virtual double getRadius() const;
      
      virtual double getGridRadius() const;

      virtual void generate();

    private:
      double      m_radius;
    };
  }
}

#endif
