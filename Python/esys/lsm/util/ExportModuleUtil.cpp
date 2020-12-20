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


#include <boost/python.hpp>
#include "Python/esys/lsm/util/Vec3Py.h"
#include "Python/esys/lsm/util/BoundingBoxPy.h"
#include "Python/esys/lsm/util/BoundingSpherePy.h"
#include "Python/esys/lsm/util/RngPy.h"
#include "Python/esys/lsm/util/QuaternionPy.h"

BOOST_PYTHON_MODULE(FoundationPy)
{
  esys::lsm::exportVec3();
  esys::lsm::exportBoundingBox();
  esys::lsm::exportBoundingSphere();
  esys::lsm::exportRng();
  esys::lsm::exportQuaternion();
}
