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


#include "Geometry/ClosePackIterator.h"

namespace esys
{
  namespace lsm
  {
    const double ClosePackIterator::SQRT_1_OVER_3 = 1.0/sqrt(3.0);
    const double ClosePackIterator::SQRT_8_OVER_3 = sqrt(8.0/3.0);
    const double ClosePackIterator::SQRT_3        = sqrt(3.0);

    Vec3L ClosePackIterator::s_orientationDimMap[NUM_ORIENTATIONS] =
      {
        Vec3L(0, 2, 1),
        Vec3L(0, 1, 2),
        Vec3L(0, 2, 1),
        Vec3L(1, 0, 2),
        Vec3L(1, 2, 0),
        Vec3L(2, 0, 1),
        Vec3L(2, 1, 0)
      };
  }
}
