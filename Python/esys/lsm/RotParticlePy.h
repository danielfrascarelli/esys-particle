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

#ifndef ESYS_LSMROTPARTICLEPY_H
#define ESYS_LSMROTPARTICLEPY_H

#include "Model/RotParticle.h"
#include "Python/esys/lsm/util/Vec3Py.h"
#include "Python/esys/lsm/util/QuaternionPy.h"

namespace esys
{
  namespace lsm
  {
    class RotParticlePy : public CRotParticle
    {
    public:
      RotParticlePy();

      RotParticlePy(const RotParticlePy &p);
      
      RotParticlePy(const CRotParticle &p);

      RotParticlePy(int id, const Vec3Py &posn, double radius, double mass);

      Vec3Py getInitialPosn() const;

      Vec3Py getPosn() const;

      void setPosn(const Vec3Py& posn);

      QuaternionPy getOrientation() const;

      void setOrientation(const QuaternionPy& quat);

      Vec3Py getVelocity() const;

      Vec3Py getLinearVelocity() const;

      void setLinearVelocity(const Vec3Py& vel);

      Vec3Py getAngularVelocity() const;

      void setAngularVelocity(const Vec3Py& vel);

      Vec3Py getAcceleration() const;

      Vec3Py getLinearAcceleration() const;

      void setLinearAcceleration(const Vec3Py& accel);

      Vec3Py getAngularAcceleration() const;

      void setAngularAcceleration(const Vec3Py& accel);

      Vec3Py getLinearForce() const;

      void setLinearForce(const Vec3Py& force);

      Vec3Py getAngularForce() const;

      void setAngularForce(const Vec3Py& force);
    };

    void exportRotParticle();
  }
}
#endif
