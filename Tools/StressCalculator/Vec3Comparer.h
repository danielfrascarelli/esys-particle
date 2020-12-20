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


#ifndef ESYS_LSMVEC3COMPARE_H
#define ESYS_LSMVEC3COMPARE_H

#include "Foundation/vec3.h"

namespace esys
{
  namespace lsm
  {
    class Vec3XyzComparer
    {
    public:
      bool operator()(const Vec3 &v1, const Vec3 &v2) const
      {
        return
          (
            (v1.X() < v2.X())
            ||
            (
              (v1.X() == v2.X())
              &&
              (
                (v1.Y() < v2.Y())
                ||
                (
                  (v1.Y() == v2.Y())
                  &&
                  (v1.Z() < v2.Z())
                )
              )
            )
          );
      }
    };
    
    class Vec3ZyxComparer
    {
    public:
      bool operator()(const Vec3 &v1, const Vec3 &v2) const
      {
        return
          (
            (v1.Z() < v2.Z())
            ||
            (
              (v1.Z() == v2.Z())
              &&
              (
                (v1.Y() < v2.Y())
                ||
                (
                  (v1.Y() == v2.Y())
                  &&
                  (v1.X() < v2.X())
                )
              )
            )
          );
      }
    };
  }
}
#endif
