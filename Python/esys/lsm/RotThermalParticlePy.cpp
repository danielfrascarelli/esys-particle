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

#include <mpi.h>
#include <boost/version.hpp>
#include <boost/python.hpp>
#include "Python/esys/lsm/RotThermalParticlePy.h"
#include <stdexcept>
namespace esys
{
  namespace lsm
  {
    RotThermalParticlePy::RotThermalParticlePy() : CRotThermParticle()
    {
    }

    RotThermalParticlePy::RotThermalParticlePy(int id, const Vec3Py &posn, double radius, double mass)
      : CRotThermParticle(radius, mass, posn, Vec3(), Vec3(), id, true)
    {
    }
    
    RotThermalParticlePy::RotThermalParticlePy(const CRotThermParticle &p) : CRotThermParticle(p)
    {
    }

    RotThermalParticlePy::RotThermalParticlePy(const RotThermalParticlePy &p) : CRotThermParticle(p)
    {
    }

    Vec3Py RotThermalParticlePy::getPosn() const
    {
      return Vec3Py(getPos());
    }

    void RotThermalParticlePy::setPosn(const Vec3Py &posn)
    {
      setPos(posn);
    }

    Vec3Py RotThermalParticlePy::getInitialPosn() const
    {
      return Vec3Py(getInitPos());
    }

    Vec3Py RotThermalParticlePy::getLinearVelocity() const
    {
      return Vec3Py(getVel());
    }

    void RotThermalParticlePy::setLinearVelocity(const Vec3Py &vel)
    {
      return setVel(vel);
    }

    Vec3Py RotThermalParticlePy::getLinearForce() const
    {
      return Vec3Py(m_force);
    }

    void RotThermalParticlePy::setLinearForce(const Vec3Py &force)
    {
      setForce(force);
    }

    Vec3Py RotThermalParticlePy::getLinearAcceleration() const
    {
      return Vec3Py(getForce()*getInvMass());
    }

    void RotThermalParticlePy::setLinearAcceleration(const Vec3Py &accel)
    {
      setForce(accel*getMass());
    }

    Vec3Py RotThermalParticlePy::getAngularVelocity() const
    {
      return Vec3Py(getAngVel());
    }

    void RotThermalParticlePy::setAngularVelocity(const Vec3Py &vel)
    {
      setAngVel(vel);
    }

    Vec3Py RotThermalParticlePy::getAngularAcceleration() const
    {
      return getMoment()*getInvInertRot();
    }

    void RotThermalParticlePy::setAngularAcceleration(const Vec3Py &accel)
    {
      setMoment(accel*getInertRot());
    }

    Vec3Py RotThermalParticlePy::getAngularForce() const
    {
      return Vec3Py(m_moment);
    }

    void RotThermalParticlePy::setAngularForce(const Vec3Py &moment)
    {
      setMoment(moment);
    }

    QuaternionPy RotThermalParticlePy::getOrientation() const
    {
      return QuaternionPy(getQuat());
    }

    void RotThermalParticlePy::setOrientation(const QuaternionPy &quat)
    {
      setQuat(quat);
    }

    using boost::python::extract;

    boost::python::tuple
    RotThermalParticlePy::PickleSuite::getstate(boost::python::object pcObj)
    {
      const RotThermalParticlePy &p = extract<const RotThermalParticlePy &>(pcObj);
      boost::python::list l;

      l.append(p.getTag());
      l.append(p.getID());
      l.append(p.getPosn());
      l.append(p.getInitialPosn());
      l.append(p.getLinearVelocity());
      l.append(p.getLinearForce());
      l.append(p.getRad());
      l.append(p.getMass());
      l.append(p.getInertRot());
      l.append(p.getAngularVelocity());
      l.append(p.getAngularForce());
      l.append(p.getOrientation());

      l.append(p.getEquilibRadius());
      l.append(p.getEquilibTemperature());
      l.append(p.getCp());
      l.append(p.getThermExpansion0());
      l.append(p.getThermExpansion1());
      l.append(p.getThermExpansion2());
      l.append(p.getTemperature());

      return boost::python::make_tuple(pcObj.attr("__dict__"), l);
    }

    void
    RotThermalParticlePy::PickleSuite::setstate(
      boost::python::object pcObj,
      boost::python::tuple state
    )
    {
      // restore the object's __dict__
      boost::python::dict d =
        boost::python::extract<boost::python::dict>(
          pcObj.attr("__dict__")
        )();
      d.update(state[0]);

      RotThermalParticlePy &p = extract<RotThermalParticlePy &>(pcObj);
      boost::python::list l = extract<boost::python::list>(state[1]);
      int i = -1;
      p.setTag(extract<int>(l[++i]));
      p.setID(extract<int>(l[++i]));
      p.setPosn(extract<const Vec3Py &>(l[++i]));
      p.setInitPos(extract<const Vec3Py &>(l[++i]));
      p.setLinearVelocity(extract<const Vec3Py &>(l[++i]));
      p.setLinearForce(extract<const Vec3Py &>(l[++i]));
      p.setRad(extract<double>(l[++i]));
      p.setMass(extract<double>(l[++i]));
      p.setInertRot(extract<double>(l[++i]));
      p.setAngularVelocity(extract<const Vec3Py &>(l[++i]));
      p.setAngularForce(extract<const Vec3Py &>(l[++i]));
      p.setOrientation(extract<const QuaternionPy &>(l[++i]));

      p.setEquilibRadius(extract<double>(l[++i]));
      p.setEquilibTemperature(extract<double>(l[++i]));
      p.setCp(extract<double>(l[++i]));
      p.setThermExpansion0(extract<double>(l[++i]));
      p.setThermExpansion1(extract<double>(l[++i]));
      p.setThermExpansion2(extract<double>(l[++i]));
      p.setTemperature(extract<double>(l[++i]));
    }

    bool RotThermalParticlePy::PickleSuite::getstate_manages_dict()
    {
      return true;
    }

    using boost::python::arg;
    void exportRotThermalParticle()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<RotThermalParticlePy>(
        "RotThermalSphere",
        "EXPERIMENTAL Rotational sphere with additional thermal properties.\n"
      )
        .def(boost::python::init<>())
        .def(boost::python::init<const RotThermalParticlePy &>())
        .def(boost::python::init<int,const Vec3Py &, double, double>(
          (
            arg("id"),
            arg("posn"),
            arg("radius"),
            arg("mass")
          ),
          "Construct a rotational spherical particle.\n"
          "@type id: int\n"
          "@kwarg id: Unique identifier for particle.\n"
          "@type posn: L{Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
          "@kwarg posn: Initial position of particle, centre-point of sphere.\n"
          "@type radius: float\n"
          "@kwarg radius: The radius of the sphere.\n"
          "@type mass: float\n"
          "@kwarg mass: Mass of particle.\n"
        ))
        .def("getId",                   &RotThermalParticlePy::getID, "Returns the unique ID of the particle\n")
        .def("getTag",                  &RotThermalParticlePy::getTag, "Returns the non-unique tag of the particle\n")
        .def("setTag",                  &RotThermalParticlePy::setTag, "Specifies the tag of the particle\n")
        .def("getInitialPosn",          &RotThermalParticlePy::getInitialPosn, "Returns the initial position of the particle\n")
        .def("getInitialPosition",      &RotThermalParticlePy::getInitialPosn, "Returns the initial position of the particle\n")
        .def("getPosn",                 &RotThermalParticlePy::getPosn, "Returns the current position of the particle\n")
        .def("getPosition",             &RotThermalParticlePy::getPosn, "Returns the current position of the particle\n")
        .def("setPosn",                 &RotThermalParticlePy::setPosn, "Specifies the position of the particle\n")
        .def("setPosition",             &RotThermalParticlePy::setPosn, "Specifies the position of the particle\n")
        .def("getOrientation",          &RotThermalParticlePy::getOrientation, "Returns the current orientation of the particle\n")
        .def("setOrientation",          &RotThermalParticlePy::setOrientation, "Specifies the orientation of the particle\n")
        .def("getLinearVelocity",       &RotThermalParticlePy::getLinearVelocity, "Returns the current translational velocity of the particle\n")
        .def("setLinearVelocity",       &RotThermalParticlePy::setLinearVelocity, "Specifies the current translational velocity of the particle\n")
        .def("getAngularVelocity",      &RotThermalParticlePy::getAngularVelocity, "Returns the current angular velocity of the particle\n")
        .def("setAngularVelocity",      &RotThermalParticlePy::setAngularVelocity, "Specifies the angular velocity of the particle\n")
        .def("getLinearAcceleration",   &RotThermalParticlePy::getLinearAcceleration, "Returns the current translational acceleration of the particle\n")
        .def("setLinearAcceleration",   &RotThermalParticlePy::setLinearAcceleration, "Specifies the translational acceleration of the particle\n")
        .def("getAngularAcceleration",  &RotThermalParticlePy::getAngularAcceleration, "Returns the current angular acceleration of the particle\n")
        .def("setAngularAcceleration",  &RotThermalParticlePy::setAngularAcceleration, "Specifies the angular acceleration of the particle\n")
        .def("getRad",                  &RotThermalParticlePy::getRad, "Returns the radius of the particle\n")
        .def("getRadius",               &RotThermalParticlePy::getRad, "Returns the radius of the particle\n")
        .def("getCentre",               &RotThermalParticlePy::getPosn, "Returns the current position of the particle\n")
        .def("getCenter",               &RotThermalParticlePy::getPosn, "Returns the current position of the particle\n")
        .def("getMass",                 &RotThermalParticlePy::getMass, "Returns the mass of the particle\n")
        .def("getEquilibRadius",        &RotThermalParticlePy::getEquilibRadius, "Returns the thermal equilibrium radius of the particle\n")
        .def("setEquilibRadius",        &RotThermalParticlePy::setEquilibRadius, "Specifies the thermal equilibrium radius of the particle\n")
        .def("getEquilibTemperature",   &RotThermalParticlePy::getEquilibTemperature, "Returns the thermal equilibrium temperature of the particle\n")
        .def("setEquilibTemperature",   &RotThermalParticlePy::setEquilibTemperature, "Specifies the thermal equilibrium temperature of the particle\n")
        .def("getTemperature",          &RotThermalParticlePy::getTemperature, "Returns the current temperature of the particle\n")
        .def("setTemperature",          &RotThermalParticlePy::setTemperature, "Specifies the temperature of the particle\n")
        .def("getCp",                   &RotThermalParticlePy::getCp, "Returns the thermal conductivity of the particle\n")
        .def("setCp",                   &RotThermalParticlePy::setCp, "Specifies the thermal conductivity of the particle\n")
        .def("getExpansionCoeff0",      &RotThermalParticlePy::getThermExpansion0, "Returns the X-component of the thermal expansion coefficient of the particle\n")
        .def("setExpansionCoeff0",      &RotThermalParticlePy::setThermExpansion0, "Specifies the X-component of the thermal expansion coefficient of the particle\n")
        .def("getExpansionCoeff1",      &RotThermalParticlePy::getThermExpansion1, "Returns the Y-component of the thermal expansion coefficient of the particle\n")
        .def("setExpansionCoeff1",      &RotThermalParticlePy::setThermExpansion1, "Specifies the Y-component of the thermal expansion coefficient of the particle\n")
        .def("getExpansionCoeff2",      &RotThermalParticlePy::getThermExpansion2, "Returns the Z-component of the thermal expansion coefficient of the particle\n")
        .def("setExpansionCoeff2",      &RotThermalParticlePy::setThermExpansion2, "Specifies the Z-component of the thermal expansion coefficient of the particle\n")
        ;
    }
  }
}
