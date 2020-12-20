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


#ifndef ESYS_LSMHEXAGBLOCKITERATOR_HPP
#define ESYS_LSMHEXAGBLOCKITERATOR_HPP

namespace esys
{
  namespace lsm
  {
    HexagBlockIterator::HexagBlockIterator()
      : ClosePackIterator()
    {
    }

    HexagBlockIterator::HexagBlockIterator(
      int numI,
      int numJ,
      int numK,
      double sphereRadius,
      ClosePackOrientation orientation
    )
      : ClosePackIterator(numI, numJ, numK, sphereRadius, orientation)
    {
      setDimRepeat(Vec3L(2,2,2));
      
      OffsetMatrix offsetMatrix;
      offsetMatrix(0,0,0) = 0.0;
      offsetMatrix(0,0,1) = getRadius();
      offsetMatrix(0,0,2) = 0.0;

      offsetMatrix(0,1,0) = getRadius();
      offsetMatrix(0,1,1) = 0.0;
      offsetMatrix(0,1,2) = 0.0;

      offsetMatrix(0,2,0) = 0.0;
      offsetMatrix(0,2,1) = 0.0;
      offsetMatrix(0,2,2) = 0.0;

      offsetMatrix(1,0,0) = 0.0;
      offsetMatrix(1,0,1) = 0.0;
      offsetMatrix(1,0,2) = 0.0;

      offsetMatrix(1,1,0) = SQRT_1_OVER_3*getRadius();
      offsetMatrix(1,1,1) = SQRT_1_OVER_3*getRadius();
      offsetMatrix(1,1,2) = SQRT_1_OVER_3*getRadius();

      offsetMatrix(1,2,0) = 0.0;
      offsetMatrix(1,2,1) = 0.0;
      offsetMatrix(1,2,2) = 0.0;

      offsetMatrix(2,0,0) = 0.0;
      offsetMatrix(2,0,1) = 0.0;
      offsetMatrix(2,0,2) = 0.0;

      offsetMatrix(2,1,0) = 0.0;
      offsetMatrix(2,1,1) = 0.0;
      offsetMatrix(2,1,2) = 0.0;

      offsetMatrix(2,2,0) = 0.0;
      offsetMatrix(2,2,1) = 0.0;
      offsetMatrix(2,2,2) = 0.0;

      setOffsetMatrix(offsetMatrix);
    }
  }
}

#endif
