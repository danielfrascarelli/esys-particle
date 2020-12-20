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


#include "Foundation/BoundingSphere.h"

namespace esys
{

  namespace lsm
  {
    BoundingSphere::BoundingSphere()
      : m_centre(Vec3::ZERO),
        m_radius(0.0)
    {
    }

    BoundingSphere::BoundingSphere(const Vec3 &centre, double radius)
      : m_centre(centre),
        m_radius(radius)
    {
    }
  
    BoundingSphere::~BoundingSphere()
    {
    }

    const Vec3 &BoundingSphere::getCentre() const
    {
      return m_centre;
    }
 
    double BoundingSphere::getRadius() const
    {
      return m_radius;
    }

    BoundingBox BoundingSphere::getBBox() const
    {
      return BoundingBox(getCentre()-getRadius(), getCentre()+getRadius());
    }

    BoundingBox BoundingSphere::get2dBBox() const
    {
      return
        BoundingBox(
          Vec3(getCentre()[0]-getRadius(), getCentre()[1]-getRadius(), 0.0),
          Vec3(getCentre()[0]+getRadius(), getCentre()[1]+getRadius(), 0.0)
        );
    }

    bool BoundingSphere::operator==(const BoundingSphere &bSphere) const
    {
      return
        (
          (getCentre() == bSphere.getCentre())
          &&
          (getRadius() == bSphere.getRadius())
        );
    }
    
    bool BoundingSphere::contains(const Vec3 &pt, double tolerance) const
    {
      const double r = (getRadius() + tolerance);
      return ((getCentre()-pt).norm2() <= (r*r));
    }

    bool BoundingSphere::contains(
      const BoundingSphere &bSphere,
      double tolerance
    ) const
    {
      const double r = (getRadius()-bSphere.getRadius() + tolerance);
      return ((getCentre()-bSphere.getCentre()).norm2() <= (r*r));
    }
    
    std::ostream &operator<<(std::ostream &oStream, const BoundingSphere &bSphere)
    {
      oStream << bSphere.getCentre() << " " << bSphere.getRadius();
      return oStream;
    }
  }
}
