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


#ifndef ESYS_LSMBOUNDINGBOXPY_H
#define ESYS_LSMBOUNDINGBOXPY_H

#include <boost/python.hpp>
#include "Foundation/BoundingBox.h"
#include "Foundation/StringUtil.h"
#include "Python/BoostPythonUtil/Util.h"
#include "Python/esys/lsm/util/Vec3Py.h"

#include <sstream>

namespace esys
{
  namespace lsm
  {
    class BoundingBoxPy : public BoundingBox
    {
    public:
      BoundingBoxPy();

      BoundingBoxPy(const Vec3Py &minPt, const Vec3Py &maxPt);
      
      BoundingBoxPy(const Vec3 &minPt, const Vec3 &maxPt);

      BoundingBoxPy(const BoundingBoxPy &v);

      BoundingBoxPy(const BoundingBox &v);

      BoundingBoxPy(const boost::python::object &pyMin, const boost::python::object &pyMax);
      
      bool operator==(const BoundingBoxPy &bBox) const;
      
      Vec3Py getMinPtPy() const;
      
      Vec3Py getMaxPtPy() const;
      
      Vec3Py getSizePy() const;

      Vec3Py getCentrePy() const;

      bool intersectsWithVec3Py(const Vec3Py &pt) const;

    };

    void exportBoundingBox();    
  }
}

std::ostream &operator<<(std::ostream &oStream, const esys::lsm::BoundingBoxPy &vec);

#endif
