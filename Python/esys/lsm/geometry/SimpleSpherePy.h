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


#ifndef ESYS_LSM_SIMPLEPARTICLEPY_H
#define ESYS_LSM_SIMPLEPARTICLEPY_H

#include <Geometry/SimpleParticle.h>
#include "Python/esys/lsm/util/Vec3Py.h"

namespace esys
{
  namespace lsm
  {
    class SimpleSpherePy : public SimpleParticle
    {
    public:
      SimpleSpherePy(
          const Vec3Py &centre,
          double radius,
          int id = 0,
          int tag = 0,
          double mass = -1.0 // indicates use of default mass
      );

      SimpleSpherePy(
        int id,
        const Vec3Py &posn,
        double radius,
        double mass
      );

      SimpleSpherePy(const SimpleSpherePy &particle);

      bool operator==(const SimpleSpherePy &p) const;

      Vec3Py getPosnPy() const;

      void setPosnPy(const Vec3Py posn);

      void translateByPy(const Vec3Py &translation);

      void rigidRotatePy(const Vec3Py &axis, const Vec3Py &pt);

      static SimpleSpherePy INVALID;
    };

    void exportSimpleSphere();
  }
}

#endif
