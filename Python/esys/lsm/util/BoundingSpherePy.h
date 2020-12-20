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


#ifndef ESYS_LSMBOUNDINGSPHEREPY_H
#define ESYS_LSMBOUNDINGSPHEREPY_H

#include <boost/python.hpp>
#include "Foundation/BoundingSphere.h"
#include "Foundation/StringUtil.h"
#include "Python/BoostPythonUtil/Util.h"
#include "Python/esys/lsm/util/Vec3Py.h"

#include <sstream>

namespace esys
{
  namespace lsm
  {
    class BoundingSpherePy : public BoundingSphere
    {
    public:
      BoundingSpherePy();

      BoundingSpherePy(const Vec3Py &centrePt, double radius);

      BoundingSpherePy(const Vec3 &centrePt, double radius);

      BoundingSpherePy(const BoundingSpherePy &v);

      BoundingSpherePy(const BoundingSphere &v);

      BoundingSpherePy(const boost::python::object &pyCentre, double radius);

      bool operator==(const BoundingSpherePy &bBox) const;

      Vec3Py getCentrePy() const;
    };

    void exportBoundingSphere();
  }
}

std::ostream &operator<<(std::ostream &oStream, const esys::lsm::BoundingSpherePy &vec);

#endif
