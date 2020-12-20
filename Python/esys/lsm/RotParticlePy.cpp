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
#include <stdexcept>
#include <boost/python.hpp>
#include "Python/esys/lsm/RotParticlePy.h"

namespace esys
{
  namespace lsm
  {
    RotParticlePy::RotParticlePy() : CRotParticle()
    {
    }

    RotParticlePy::RotParticlePy(const RotParticlePy &p) : CRotParticle(p)
    {
    }

    RotParticlePy::RotParticlePy(const CRotParticle &p) : CRotParticle(p)
    {
    }

    // dynamic & rotational default to true
    RotParticlePy::RotParticlePy(int id, const Vec3Py &posn, double radius, double mass)
      : CRotParticle(radius, mass, posn, Vec3(), Vec3(), id,true,true)
    {
    }

    Vec3Py RotParticlePy::getInitialPosn() const
    {
      return Vec3Py(getInitPos());
    }

    Vec3Py RotParticlePy::getPosn() const
    {
      return Vec3Py(getPos());
    }

    void RotParticlePy::setPosn(const Vec3Py& posn)
    {
      setPos(posn);
    }

    QuaternionPy RotParticlePy::getOrientation() const
    {
      return QuaternionPy(getQuat());
    }

    void RotParticlePy::setOrientation(const QuaternionPy& quat)
    {
      setQuat(quat);
    }

    Vec3Py RotParticlePy::getVelocity() const
    {
      return Vec3Py(getVel());
    }

    Vec3Py RotParticlePy::getLinearVelocity() const
    {
      return Vec3Py(getVel());
    }

    void RotParticlePy::setLinearVelocity(const Vec3Py& vel)
    {
      return setVel(vel);
    }

    Vec3Py RotParticlePy::getAngularVelocity() const
    {
      return Vec3Py(getAngVel());
    }

    void RotParticlePy::setAngularVelocity(const Vec3Py& vel)
    {
      setAngVel(vel);
    }

    Vec3Py RotParticlePy::getAcceleration() const
    {
      return Vec3Py(getForce()*getInvMass());
    }

    Vec3Py RotParticlePy::getLinearAcceleration() const
    {
      return Vec3Py(getForce()*getInvMass());
    }

    void RotParticlePy::setLinearAcceleration(const Vec3Py& accel)
    {
      setForce(accel*getMass());
    }

    Vec3Py RotParticlePy::getAngularAcceleration() const
    {
      return getMoment()*getInvInertRot();
    }

    void RotParticlePy::setAngularAcceleration(const Vec3Py& accel)
    {
      setMoment(accel*getInertRot());
    }

    Vec3Py RotParticlePy::getLinearForce() const
    {
      return Vec3Py(m_force);
    }

    void RotParticlePy::setLinearForce(const Vec3Py& force)
    {
      setForce(force);
    }

    Vec3Py RotParticlePy::getAngularForce() const
    {
      return Vec3Py(m_moment);
    }

    void RotParticlePy::setAngularForce(const Vec3Py& moment)
    {
      setMoment(moment);
    }

    using boost::python::arg;
    void exportRotParticle()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<RotParticlePy>("RotSphere","Class defining the properties of spheres with translational and rotational degrees of freedom.\n")
        .def(boost::python::init<>())
        .def(boost::python::init<const RotParticlePy &>())
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
        .def("getId",                   &RotParticlePy::getID, "Returns the unique ID of the particle\n")
        .def("getTag",                  &RotParticlePy::getTag, "Returns the non-unique tag of the particle\n")
        .def("setTag",                  &RotParticlePy::setTag, "Specifies the tag of the particle\n")
        .def("getInitialPosn",          &RotParticlePy::getInitialPosn, "Returns the initial position of the particle\n")
        .def("getInitialPosition",      &RotParticlePy::getInitialPosn, "Returns the initial position of the particle\n")
        .def("getPosn",                 &RotParticlePy::getPosn, "Returns the current position of the particle\n")
        .def("getPosition",             &RotParticlePy::getPosn, "Returns the current position of the particle\n")
        .def("setPosn",                 &RotParticlePy::setPosn, "Specifies the position of the particle\n")
        .def("setPosition",             &RotParticlePy::setPosn, "Specifies the position of the particle\n")
        .def("getOrientation",          &RotParticlePy::getOrientation, "Returns the current orientation of the particle\n")
        .def("setOrientation",          &RotParticlePy::setOrientation, "Specifies the orientation of the particle\n")
        .def("getVelocity",             &RotParticlePy::getVelocity, "Returns the current translational velocity of the particle\n")
        .def("getLinearVelocity",       &RotParticlePy::getLinearVelocity, "Returns the current translational velocity of the particle\n")
        .def("setLinearVelocity",       &RotParticlePy::setLinearVelocity, "Specifies the current translational velocity of the particle\n")
        .def("getAngularVelocity",      &RotParticlePy::getAngularVelocity, "Returns the current angular velocity of the particle\n")
        .def("setAngularVelocity",      &RotParticlePy::setAngularVelocity, "Specifies the angular velocity of the particle\n")
        .def("getAcceleration",         &RotParticlePy::getAcceleration, "Returns the current translational acceleration of the particle\n")
        .def("getLinearAcceleration",   &RotParticlePy::getLinearAcceleration, "Returns the current translational acceleration of the particle\n")
        .def("setLinearAcceleration",   &RotParticlePy::setLinearAcceleration, "Specifies the translational acceleration of the particle\n")
        .def("getAngularAcceleration",  &RotParticlePy::getAngularAcceleration, "Returns the current angular acceleration of the particle\n")
        .def("setAngularAcceleration",  &RotParticlePy::setAngularAcceleration, "Specifies the angular acceleration of the particle\n")
        .def("getRad",                  &RotParticlePy::getRad, "Returns the radius of the particle\n")
        .def("getRadius",               &RotParticlePy::getRad, "Returns the radius of the particle\n")
        .def("getCentre",               &RotParticlePy::getPosn, "Returns the current position of the particle\n")
        .def("getCenter",               &RotParticlePy::getPosn, "Returns the current position of the particle\n")
        .def("getMass",                 &RotParticlePy::getMass, "Returns the mass of the particle\n")
      ;
    }
  }
}
