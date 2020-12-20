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


#ifndef ESYS_LSMSPHEREBOXVOLCALCULATOR_H
#define ESYS_LSMSPHEREBOXVOLCALCULATOR_H

#include "Foundation/vec3.h"
#include "Geometry/IntersectionVolCalculator.h"

namespace esys
{
  namespace lsm
  {
    /**
     * Calculates the volume of intersection between a box and
     * a sphere.
     */
    class SphereBoxVolCalculator : public impl::IntersectionVolCalculator<3, Vec3>
    {
    public:
      typedef impl::IntersectionVolCalculator<3, Vec3> Inherited;
      typedef Inherited::BasicSphere                   Sphere;
      typedef Inherited::BasicBox                      Box;
      SphereBoxVolCalculator(
        const Box    &box
      ) : Inherited(box)
      {
      }
      
      double getVolume(const Sphere &sphere)
      {
        return Inherited::getVolume(sphere);
      }

    private:
    };
  }
}
#endif
