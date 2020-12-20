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


#ifndef ESYS_LSMCUBICBLOCKITERATOR_HPP
#define ESYS_LSMCUBICBLOCKITERATOR_HPP

namespace esys
{
  namespace lsm
  {
    CubicBlockIterator::CubicBlockIterator()
      : ClosePackIterator()
    {
    }

    CubicBlockIterator::CubicBlockIterator(
      int numI,
      int numJ,
      int numK,
      double sphereRadius,
      ClosePackOrientation orientation
    )
      : ClosePackIterator(numI, numJ, numK, sphereRadius, orientation)
    {
      setDimRepeat(Vec3L(6,3,3));
      
      OffsetMatrix offsetMatrix;
      offsetMatrix(0,0,0) = 0.0;
      offsetMatrix(0,0,1) = 0.0;
      offsetMatrix(0,0,2) = getRadius();
      offsetMatrix(0,0,3) = 0.0;
      offsetMatrix(0,0,4) = 0.0;
      offsetMatrix(0,0,5) = getRadius();

      offsetMatrix(0,1,0) = getRadius();
      offsetMatrix(0,1,1) = getRadius();
      offsetMatrix(0,1,2) = 0.0;
      offsetMatrix(0,1,3) = getRadius();
      offsetMatrix(0,1,4) = getRadius();
      offsetMatrix(0,1,5) = 0.0;

      offsetMatrix(0,2,0) = 0.0;
      offsetMatrix(0,2,1) = 0.0;
      offsetMatrix(0,2,2) = getRadius();
      offsetMatrix(0,2,3) = 0.0;
      offsetMatrix(0,2,4) = 0.0;
      offsetMatrix(0,2,5) = getRadius();

      offsetMatrix(0,3,0) = getRadius();
      offsetMatrix(0,3,1) = getRadius();
      offsetMatrix(0,3,2) = 0.0;
      offsetMatrix(0,3,3) = getRadius();
      offsetMatrix(0,3,4) = getRadius();
      offsetMatrix(0,3,5) = 0.0;

      offsetMatrix(0,4,0) = 0.0;
      offsetMatrix(0,4,1) = 0.0;
      offsetMatrix(0,4,2) = getRadius();
      offsetMatrix(0,4,3) = 0.0;
      offsetMatrix(0,4,4) = 0.0;
      offsetMatrix(0,4,5) = getRadius();

      offsetMatrix(0,5,0) = getRadius();
      offsetMatrix(0,5,1) = getRadius();
      offsetMatrix(0,5,2) = 0.0;
      offsetMatrix(0,5,3) = getRadius();
      offsetMatrix(0,5,4) = getRadius();
      offsetMatrix(0,5,5) = 0.0;

      offsetMatrix(1,0,0) = 0.0;
      offsetMatrix(1,0,1) = 0.0;
      offsetMatrix(1,0,2) = 0.0;

      offsetMatrix(1,1,0) = 2.0*SQRT_1_OVER_3*getRadius();
      offsetMatrix(1,1,1) = 2.0*SQRT_1_OVER_3*getRadius();
      offsetMatrix(1,1,2) = 2.0*SQRT_1_OVER_3*getRadius();

      offsetMatrix(1,2,0) = SQRT_1_OVER_3*getRadius();
      offsetMatrix(1,2,1) = SQRT_1_OVER_3*getRadius();
      offsetMatrix(1,2,2) = SQRT_1_OVER_3*getRadius();

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
