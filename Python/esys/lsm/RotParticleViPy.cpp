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
#include "Python/esys/lsm/RotParticleViPy.h"
#include <stdexcept>
namespace esys
{
  namespace lsm
  {
    RotParticleViPy::RotParticleViPy() : CRotParticleVi()
    {
    }

    RotParticleViPy::RotParticleViPy(int id, const Vec3Py &posn, double radius, double mass)
      : CRotParticleVi(radius, mass, posn, Vec3(), Vec3(), id, true)
    {
    }
    
    RotParticleViPy::RotParticleViPy(const CRotParticleVi &p) : CRotParticleVi(p)
    {
    }

    RotParticleViPy::RotParticleViPy(const RotParticleViPy &p) : CRotParticleVi(p)
    {
    }

    Vec3Py RotParticleViPy::getPosn() const
    {
      return Vec3Py(getPos());
    }

    void RotParticleViPy::setPosn(const Vec3Py &posn)
    {
      setPos(posn);
    }

    Vec3Py RotParticleViPy::getInitialPosn() const
    {
      return Vec3Py(getInitPos());
    }

    Vec3Py RotParticleViPy::getLinearVelocity() const
    {
      return Vec3Py(getVel());
    }

    Vec3Py RotParticleViPy::getVelocity() const
    {
      return Vec3Py(getVel());
    }

    void RotParticleViPy::setLinearVelocity(const Vec3Py &vel)
    {
      return setVel(vel);
    }

    Vec3Py RotParticleViPy::getLinearForce() const
    {
      return Vec3Py(m_force);
    }

    void RotParticleViPy::setLinearForce(const Vec3Py &force)
    {
      setForce(force);
    }

    Vec3Py RotParticleViPy::getLinearAcceleration() const
    {
      return Vec3Py(getForce()*getInvMass());
    }

    Vec3Py RotParticleViPy::getAcceleration() const
    {
      return Vec3Py(getForce()*getInvMass());
    }

    void RotParticleViPy::setLinearAcceleration(const Vec3Py &accel)
    {
      setForce(accel*getMass());
    }

    Vec3Py RotParticleViPy::getAngularVelocity() const
    {
      return Vec3Py(getAngVel());
    }

    void RotParticleViPy::setAngularVelocity(const Vec3Py &vel)
    {
      setAngVel(vel);
    }

    Vec3Py RotParticleViPy::getAngularVelocityT() const
    {
      return Vec3Py(getAngVel_t());
    }

    void RotParticleViPy::setAngularVelocityT(const Vec3Py &vel)
    {
      setAngVel_t(vel);
    }

    Vec3Py RotParticleViPy::getAngularAcceleration() const
    {
      return getMoment()*getInvInertRot();
    }

    void RotParticleViPy::setAngularAcceleration(const Vec3Py &accel)
    {
      setMoment(accel*getInertRot());
    }

    Vec3Py RotParticleViPy::getAngularForce() const
    {
      return Vec3Py(m_moment);
    }

    void RotParticleViPy::setAngularForce(const Vec3Py &moment)
    {
      setMoment(moment);
    }

    QuaternionPy RotParticleViPy::getOrientation() const
    {
      return QuaternionPy(getQuat());
    }

    void RotParticleViPy::setOrientation(const QuaternionPy &quat)
    {
      setQuat(quat);
    }

    using boost::python::extract;

    boost::python::tuple
    RotParticleViPy::PickleSuite::getstate(boost::python::object pcObj)
    {
      const RotParticleViPy &p = extract<const RotParticleViPy &>(pcObj);
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

      return boost::python::make_tuple(pcObj.attr("__dict__"), l);
    }

    void
    RotParticleViPy::PickleSuite::setstate(
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

      RotParticleViPy &p = extract<RotParticleViPy &>(pcObj);
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
    }

    bool RotParticleViPy::PickleSuite::getstate_manages_dict()
    {
      return true;
    }

    using boost::python::arg;
    void exportRotParticleVi()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<RotParticleViPy>(
        "RotSphereVi",
        "EXPERIMENTAL Rotational sphere using Verlet 2nd order time integration.\n"
      )
        .def(boost::python::init<>())
        .def(boost::python::init<const RotParticleViPy &>())
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
        .def("getId",                   &RotParticleViPy::getID, "Returns the unique ID of the particle\n")
        .def("getTag",                  &RotParticleViPy::getTag, "Returns the non-unique tag of the particle\n")
        .def("setTag",                  &RotParticleViPy::setTag, "Specifies the tag of the particle\n")
        .def("getInitialPosn",          &RotParticleViPy::getInitialPosn, "Returns the initial position of the particle\n")
        .def("getInitialPosition",      &RotParticleViPy::getInitialPosn, "Returns the initial position of the particle\n")
        .def("getPosn",                 &RotParticleViPy::getPosn, "Returns the current position of the particle\n")
        .def("getPosition",             &RotParticleViPy::getPosn, "Returns the current position of the particle\n")
        .def("setPosn",                 &RotParticleViPy::setPosn, "Specifies the position of the particle\n")
        .def("setPosition",             &RotParticleViPy::setPosn, "Specifies the position of the particle\n")
        .def("getOrientation",          &RotParticleViPy::getOrientation, "Returns the current orientation of the particle\n")
        .def("setOrientation",          &RotParticleViPy::setOrientation, "Specifies the orientation of the particle\n")
        .def("getVelocity",             &RotParticleViPy::getVelocity, "Returns the current translational velocity of the particle\n")
        .def("getLinearVelocity",       &RotParticleViPy::getLinearVelocity, "Returns the current translational velocity of the particle\n")
        .def("setLinearVelocity",       &RotParticleViPy::setLinearVelocity, "Specifies the current translational velocity of the particle\n")
        .def("getAngularVelocity",      &RotParticleViPy::getAngularVelocity, "Returns the angular velocity of the particle at time t - 0.5*dt\n")
        .def("setAngularVelocity",      &RotParticleViPy::setAngularVelocity, "Specifies the angular velocity of the particle at time t - 0.5*dt\n")
        .def("getAngularVelocityT",   &RotParticleViPy::getAngularVelocityT, "Returns the angular velocity of the particle at time t\n")
        .def("setAngularVelocityT",   &RotParticleViPy::setAngularVelocityT, "Specifies the angular velocity of the particle at time t\n")
        .def("getAcceleration",         &RotParticleViPy::getAcceleration, "Returns the current translational acceleration of the particle\n")
        .def("getLinearAcceleration",   &RotParticleViPy::getLinearAcceleration, "Returns the current translational acceleration of the particle\n")
        .def("setLinearAcceleration",   &RotParticleViPy::setLinearAcceleration, "Specifies the translational acceleration of the particle\n")
        .def("getAngularAcceleration",  &RotParticleViPy::getAngularAcceleration, "Returns the current angular acceleration of the particle\n")
        .def("setAngularAcceleration",  &RotParticleViPy::setAngularAcceleration, "Specifies the angular acceleration of the particle\n")
        .def("getRad",                  &RotParticleViPy::getRad, "Returns the radius of the particle\n")
        .def("getRadius",               &RotParticleViPy::getRad, "Returns the radius of the particle\n")
        .def("getCentre",               &RotParticleViPy::getPosn, "Returns the current position of the particle\n")
        .def("getCenter",               &RotParticleViPy::getPosn, "Returns the current position of the particle\n")
        .def("getMass",                 &RotParticleViPy::getMass, "Returns the mass of the particle\n")
        ;
    }
  }
}
