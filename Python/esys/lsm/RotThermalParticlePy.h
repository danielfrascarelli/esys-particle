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

#ifndef ESYS_LSMTHERMALROTPARTICLEPY_H
#define ESYS_LSMTHERMALROTPARTICLEPY_H

#include "Model/RotThermParticle.h"
#include "Python/esys/lsm/util/Vec3Py.h"
#include "Python/esys/lsm/util/QuaternionPy.h"

namespace esys
{
  namespace lsm
  {
    class RotThermalParticlePy : public CRotThermParticle
    {
    public:
      RotThermalParticlePy();

      RotThermalParticlePy(const RotThermalParticlePy &p);
      
      RotThermalParticlePy(const CRotThermParticle &p);

      RotThermalParticlePy(int id, const Vec3Py &posn, double radius, double mass);

      Vec3Py getPosn() const;

      void setPosn(const Vec3Py &posn);

      Vec3Py getInitialPosn() const;

      Vec3Py getLinearVelocity() const;

      void setLinearVelocity(const Vec3Py &vel);

      Vec3Py getLinearForce() const;

      void setLinearForce(const Vec3Py &force);

      Vec3Py getLinearAcceleration() const;

      void setLinearAcceleration(const Vec3Py &accel);

      Vec3Py getAngularVelocity() const;

      void setAngularVelocity(const Vec3Py &vel);

      Vec3Py getAngularForce() const;

      void setAngularForce(const Vec3Py &force);

      Vec3Py getAngularAcceleration() const;

      void setAngularAcceleration(const Vec3Py &accel);

      QuaternionPy getOrientation() const;

      void setOrientation(const QuaternionPy &quat);

      class PickleSuite : public boost::python::pickle_suite
      {
      public:
        static
        boost::python::tuple
        getstate(boost::python::object pcObj);

        static
        void
        setstate(boost::python::object pcObj, boost::python::tuple state);

        static bool getstate_manages_dict();
      };
    };

    void exportRotThermalParticle();
  }
}
#endif
