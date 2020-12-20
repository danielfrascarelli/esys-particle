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


#ifndef ESYS_LSMCUBICBLOCK_HPP
#define ESYS_LSMCUBICBLOCK_HPP

namespace esys
{
  namespace lsm
  {
    template <typename TmplParticle>
    CubicBlock<TmplParticle>::CubicBlock(
      unsigned int numX,
      unsigned int numY,
      unsigned int numZ,
      double radius,
      ClosePackOrientation orientation
    )
      : Inherited(numX, numY, numZ, radius, orientation)
    {
    }
  }
}

#endif
