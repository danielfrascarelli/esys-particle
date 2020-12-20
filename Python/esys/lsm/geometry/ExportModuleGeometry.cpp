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
#include "Python/esys/lsm/geometry/SimpleSpherePy.h"
#include "Python/esys/lsm/geometry/SimpleSphereCollectionPy.h"
#include "Python/esys/lsm/geometry/ParticleCollectionPy.h"
#include "Python/esys/lsm/geometry/GrainPy.h"
#include "Python/esys/lsm/geometry/GrainCollectionPy.h"
#include "Python/esys/lsm/geometry/PackerPy.h"
#include "Python/esys/lsm/geometry/SimpleBlockPy.h"
#include "Python/esys/lsm/geometry/CubicBlockPy.h"
#include "Python/esys/lsm/geometry/HexagBlockPy.h"
#include "Python/esys/lsm/geometry/TaggedIdConnectionPy.h"
#include "Python/esys/lsm/geometry/DistConnectionsPy.h"
#include "Python/esys/lsm/geometry/GougeConfigPy.h"
#include "Python/esys/lsm/geometry/GougeConfigPrmsPy.h"
#include "Python/esys/lsm/geometry/MiscPy.h"
#include "Python/esys/lsm/geometry/SimpleSphereNeighboursPy.h"

BOOST_PYTHON_MODULE(GeometryPy)
{
  esys::lsm::exportOrientation();
  esys::lsm::exportGougeConfigPrms();
  esys::lsm::exportGougeConfig();
  esys::lsm::exportSimpleSphere();
  esys::lsm::exportSimpleSphereCollection();
  esys::lsm::exportParticleCollection();
  esys::lsm::exportGrain();
  esys::lsm::exportGrainCollection();
  esys::lsm::exportPacker();
  esys::lsm::exportCubicBlock();
  esys::lsm::exportHexagBlock();
  esys::lsm::exportSimpleBlock();
  esys::lsm::exportTaggedIdConnection();
  esys::lsm::exportDistConnections();
  esys::lsm::exportMisc();
  esys::lsm::exportSimpleSphereNeighbours();
}
