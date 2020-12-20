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


#ifndef ESYS_LSMHEXAGBLOCK_H
#define ESYS_LSMHEXAGBLOCK_H

#include "Geometry/SimpleParticle.h"
#include "Geometry/HexagBlockIterator.h"
#include "Geometry/ClosePackBlock.h"

namespace esys
{
  namespace lsm
  {
    template <typename TmplParticle = SimpleParticle>
    class HexagBlock : public ClosePackBlock<HexagBlockIterator,TmplParticle>
    {
    public:
      typedef ClosePackBlock<HexagBlockIterator,TmplParticle> Inherited;
      HexagBlock(
          unsigned int numX,
          unsigned int numY,
          unsigned int numZ,
          double radius = 0.5,
          ClosePackOrientation orientation = DEFAULT_ORIENT
      );
    };
  }
}

#include "Geometry/HexagBlock.hpp"

#endif
