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


#ifndef ESYS_LSMSIMPLEBLOCKPY_H
#define ESYS_LSMSIMPLEBLOCKPY_H

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
    class SimpleBlockPy : public SimpleSphereCollectionPy
    {
    public:
      SimpleBlockPy();
      SimpleBlockPy(const boost::python::list &dimCount, double radius=0.5);
    private:
    };

    void exportSimpleBlock();
  }
}

#endif
