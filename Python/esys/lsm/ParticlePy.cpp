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

#include <boost/version.hpp>
#include <boost/python.hpp>
#include "Python/esys/lsm/ParticlePy.h"

namespace esys
{
  namespace lsm
  {
    ParticlePy::ParticlePy() : CParticle()
    {
    }

    ParticlePy::ParticlePy(int id, const Vec3Py &posn, double radius, double mass)
      : CParticle(radius, mass, posn, Vec3(), Vec3(), id, true)
    {
    }

    ParticlePy::ParticlePy(const CParticle &p) : CParticle(p)
    {
    }

    ParticlePy::ParticlePy(const ParticlePy &p) : CParticle(p)
    {
    }

    Vec3Py ParticlePy::getPosn() const
    {
      return Vec3Py(getPos());
    }

    void ParticlePy::setPosn(const Vec3Py &posn)
    {
      setPos(posn);
    }

    Vec3Py ParticlePy::getVelocity() const
    {
      return Vec3Py(getVel());
    }

    void ParticlePy::setVelocity(const Vec3Py &vel)
    {
      setVel(vel);
    }

    Vec3Py ParticlePy::getInitialPosn() const
    {
      return Vec3Py(getInitPos());
    }

    Vec3Py ParticlePy::getAcceleration() const
    {
      return Vec3Py(getForce()*getInvMass());
    }

    void ParticlePy::setAcceleration(const Vec3Py &accel)
    {
      setForce(accel*getMass());
    }

    Vec3Py ParticlePy::getLForce() const
    {
      return Vec3Py(getForce());
    }

    void ParticlePy::setForce(const Vec3Py &force)
    {
      setForce(force);
    }

    using boost::python::arg;
    void exportParticle()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<ParticlePy>("NRotSphere","Class defining the properties of non-rotational spheres.\n")
        
        .def(boost::python::init<>())
        .def(boost::python::init<const ParticlePy &>())
        .def(
          boost::python::init<int,const Vec3Py &, double, double>(
            (
              arg("id"),
              arg("posn"),
              arg("radius"),
              arg("mass")
            ),
            "Construct a non-rotational spherical particle.\n"
            "@type id: int\n"
            "@kwarg id: Unique identifier for particle.\n"
            "@type posn: L{Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
            "@kwarg posn: Initial position of particle, centre-point of sphere."
            "@type radius: float\n"
            "@kwarg radius: The radius of the sphere.\n"
            "@type mass: float\n"
            "@kwarg mass: Mass of particle."
          )
        )
        .def("getId",           	&ParticlePy::getID, "Returns the unique ID of the particle\n")
        .def("getTag",          	&ParticlePy::getTag, "Returns the non-unique tag of the particle\n")
        .def("setTag",          	&ParticlePy::setTag, "Specifies the particle's tag\n")
        .def("getPosn",         	&ParticlePy::getPosn, "Returns the current position of the particle's centre of mass\n")
        .def("getPosition",     	&ParticlePy::getPosn, "Returns the current position of the particle's centre of mass\n")
        .def("setPosn",        		&ParticlePy::setPosn, "Specifies the position of the particle's centre of mass\n")
        .def("setPosition",     	&ParticlePy::setPosn, "Specifies the position of the particle's centre of mass\n")
        .def("getInitialPosn",		&ParticlePy::getInitialPosn, "Returns the initial position of the particle's centre of mass\n")
        .def("getInitialPosition",	&ParticlePy::getInitialPosn, "Returns the initial position of the particle's centre of mass\n")
        .def("getVelocity",     	&ParticlePy::getVelocity, "Returns the current velocity of the particle\n")
        .def("getLinearVelocity",     	&ParticlePy::getVelocity, "Returns the current velocity of the particle\n")
        .def("setVelocity",     	&ParticlePy::setVelocity, "Specifies the velocity of the particle\n")
        .def("setLinearVelocity",     	&ParticlePy::setVelocity, "Specifies the velocity of the particle\n")
        .def("getAcceleration", 	&ParticlePy::getAcceleration, "Returns the current acceleration of the particle\n")
        .def("getLinearAcceleration",	&ParticlePy::getAcceleration, "Returns the current acceleration of the particle\n")
        .def("setAcceleration",	 	&ParticlePy::setAcceleration, "Specifies the acceleration of the particle\n")
        .def("setLinearAcceleration",	&ParticlePy::setAcceleration, "Specifies the acceleration of the particle\n")
        .def("getRadius",       	&ParticlePy::getRad, "Returns the radius of the particle\n")
        .def("getRad",          	&ParticlePy::getRad, "Returns the radius of the particle\n")
        .def("getCenter",       	&ParticlePy::getPosn, "Returns the current position of the particle's centre of mass\n")
        .def("getCentre",       	&ParticlePy::getPosn, "Returns the current position of the particle's centre of mass\n")
        .def("getMass",         	&ParticlePy::getMass, "Returns the mass of the particle\n")
        .def("getForce",        	&ParticlePy::getLForce, "Returns the current net force acting on the particle\n")
        .def("getLinearForce",        	&ParticlePy::getLForce, "Returns the current net force acting on the particle\n")
        .def("setForce",        	&ParticlePy::setForce, "Specifies the net force acting on the particle\n")
        .def("setLinearForce",        	&ParticlePy::setForce, "Specifies the net force acting on the particle\n")
        ;
    }
  }
}
