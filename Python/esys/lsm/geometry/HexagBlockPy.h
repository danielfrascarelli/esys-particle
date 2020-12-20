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


#ifndef ESYS_LSMHEXAGBLOCKPY_H
#define ESYS_LSMHEXAGBLOCKPY_H

#include "Geometry/ClosePackOrientation.h"
#include "Python/esys/lsm/geometry/SimpleSphereCollectionPy.h"

namespace boost
{
  namespace python
  {
    class list;
  }
}

namespace esys
{
  namespace lsm
  {
    class HexagBlockPy : public SimpleSphereCollectionPy
    {
    public:
      HexagBlockPy();

      HexagBlockPy(
        const boost::python::list &dimCount,
        double radius,
        const ClosePackOrientation &orientation
      );
    private:
    };

    void exportHexagBlock();
  }
}

#endif
