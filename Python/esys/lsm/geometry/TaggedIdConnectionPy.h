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


#ifndef ESYS_LSMTAGGEDIDCONNECTIONPY_H
#define ESYS_LSMTAGGEDIDCONNECTIONPY_H

#include "Geometry/BasicInteraction.h"

namespace esys
{
  namespace lsm
  {
    class TaggedIdConnectionPy : public BasicInteraction
    {
    public:
      TaggedIdConnectionPy(Id id1, Id id2, Tag tag=0);

      TaggedIdConnectionPy(const TaggedIdConnectionPy &idConn);
    };
    void exportTaggedIdConnection();
  }
}

#endif
