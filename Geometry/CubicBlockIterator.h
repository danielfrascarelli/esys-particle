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


#ifndef ESYS_LSMCUBICBLOCKITERATOR_H
#define ESYS_LSMCUBICBLOCKITERATOR_H

#include "Foundation/BoundingBox.h"
#include "Foundation/vec3.h"
#include "Geometry/ClosePackIterator.h"
#include "Geometry/Vec3L.h"

namespace esys
{
  namespace lsm
  {
    
    /**
     * Class for iterating over the centre-points of spheres arranged
     * in a cubic-close-packed (face-centred-close packed).
     */
    class CubicBlockIterator : public ClosePackIterator
    {
    public:

      /**
       * Creates default empty iterator.
       */
      inline CubicBlockIterator();

      /**
       * Creates an iterator which will iterate over numI*numJ*numK
       * centre points of spheres with radius sphereRadius.
       * @param numI number of spheres in the i direction.
       * @param numJ number of spheres in the j direction.
       * @param numK number of spheres in the k direction.
       * @param sphereRadius radius of spheres in the packing.
       * @param orientation specifies the axis alignment of layers.
       */
      inline CubicBlockIterator(
        int numI,
        int numJ,
        int numK,
        double sphereRadius,
        ClosePackOrientation orientation = DEFAULT_ORIENT
      );
    };
  }
}

#include "Geometry/CubicBlockIterator.hpp"

#endif
