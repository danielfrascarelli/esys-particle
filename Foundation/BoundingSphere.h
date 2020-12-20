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


#ifndef ESYS_LSMBOUNDINGSPHERE_H
#define ESYS_LSMBOUNDINGSPHERE_H

#include "Foundation/vec3.h"
#include "Foundation/BoundingBox.h"

namespace esys
{
  namespace lsm
  {
    /**
     *
     */
    class BoundingSphere
    {
    public:
      inline BoundingSphere();
      
      inline BoundingSphere(const Vec3 &centre, double radius);

      inline virtual ~BoundingSphere();

      inline const Vec3 &getCentre() const;

      inline double getRadius() const;

      inline BoundingBox getBBox() const;

      inline BoundingBox get2dBBox() const;
      
      inline bool operator==(const BoundingSphere &bSphere) const;
      
      inline bool contains(const Vec3 &pt, double tolerance = 0.0) const;

      inline bool contains(
        const BoundingSphere &bSphere,
        double tolerance = 0.0
      ) const;

    private:
      Vec3   m_centre;
      double m_radius;
    };
    inline std::ostream &operator<<(std::ostream &oStream, const BoundingSphere &bSphere);
  }
}

#include "Foundation/BoundingSphere.hpp"

#endif
