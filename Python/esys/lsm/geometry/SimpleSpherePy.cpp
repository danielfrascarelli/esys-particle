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
#include <Geometry/SimpleParticle.h>
#include "Python/esys/lsm/geometry/SimpleSpherePy.h"

namespace esys
{
  namespace lsm
  {
    SimpleSpherePy SimpleSpherePy::INVALID = SimpleSpherePy(Vec3Py(Vec3::ZERO), 0.0, -1, -1, -1.0);

    SimpleSpherePy::SimpleSpherePy(
        const Vec3Py &center,
        double radius,
        int id,
        int tag,
        double mass
    ) : SimpleParticle(center, radius, id, tag)
    {
      if (mass >= 0.0)
      {
        setMass(mass);
      }
    }

    SimpleSpherePy::SimpleSpherePy(
      int id,
      const Vec3Py &posn,
      double radius,
      double mass
    )
      : SimpleParticle(posn, radius, id, 0)
    {
      if (mass >= 0.0)
      {
        setMass(mass);
      }
    }

    SimpleSpherePy::SimpleSpherePy(const SimpleSpherePy &particle)
      : SimpleParticle(particle)
    {
    }

    void SimpleSpherePy::setPosnPy(const Vec3Py posn)
    {
      setPos(posn);
    }

    bool SimpleSpherePy::operator==(const SimpleSpherePy &p) const
    {
      return (getID() == p.getID());
    }

    Vec3Py SimpleSpherePy::getPosnPy() const
    {
      return Vec3Py(getPos());
    }

    void SimpleSpherePy::translateByPy(const Vec3Py &translation)
    {
      translateBy(translation);
    }

    void SimpleSpherePy::rigidRotatePy(const Vec3Py &axis, const Vec3Py &pt)
    {
      rotate(axis, pt);
    }

    class SimpleSpherePyPickleSuite : public boost::python::pickle_suite
    {
    public:
      static
      boost::python::tuple
      getinitargs(SimpleSpherePy const& p)
      {
        return
          boost::python::make_tuple(
            p.getPosnPy(),
            p.getRadius(),
            p.getId(),
            p.getTag(),
            p.getMass()
          );
      }
    };

    using boost::python::arg;
    void exportSimpleSphere()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<SimpleSpherePy>(
        "SimpleSphere",
        boost::python::init<
          const Vec3Py &,
          double,
          int,
          int,
          double
        >(
          (
            arg("centre")=Vec3Py(Vec3::ZERO),
            arg("radius")=1.0,
            arg("id")=0,
            arg("tag")=0,
            arg("mass")=-1.0
          )
        )
      )
      .def(
        boost::python::init<
          const SimpleSpherePy &
        >(
          (
            arg("p")
          )
        )
      )
      .def(
        boost::python::init<
          int,
          const Vec3Py &,
          double,
          double
        >(
          (
            arg("id"),
            arg("posn"),
            arg("radius"),
            arg("mass")
          ),
          "Constructs a spherical particle.\n"
          "@type centre: L{esys.lsm.util.Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
          "@kwarg centre: Centre point of the spherical particle.\n"
          "@type radius: float\n"
          "@kwarg radius: Radius of the spherical particle.\n"
          "@type id: int\n"
          "@kwarg id: Unique identifier for the particle.\n"
          "@type tag: int\n"
          "@kwarg tag: Label for the particle.\n"
          "@type mass: float\n"
          "@kwarg mass: The mass of the spherical particle.\n"
        )
      )
      .def(
        "getRadius",
        &SimpleSpherePy::getRad,
        "Returns the radius of this sphere.\n"
        "@rtype: float\n"
        "@return: Radius of this sphere."
      )
      .def(
        "getRad",
        &SimpleSpherePy::getRad,
        "Returns the radius of this sphere.\n"
        "@rtype: float\n"
        "@return: Radius of this sphere."
      )
      .def(
        "getCenter",
        &SimpleSpherePy::getPosnPy,
        "Returns the centre point of this sphere.\n"
        "@rtype: L{esys.lsm.util.Vec3<esys.lsm.util.FoundationPy>}\n"
        "@return: Returns the centre point of this sphere.\n"
      )
      .def(
        "getCentre",
        &SimpleSpherePy::getPosnPy,
        "Returns the centre point of this sphere.\n"
        "@rtype: L{esys.lsm.util.Vec3<esys.lsm.util.FoundationPy>}\n"
        "@return: Returns the centre point of this sphere.\n"
      )
      .def(
        "getPosn",
        &SimpleSpherePy::getPosnPy,
        "Returns the centre point of this sphere.\n"
        "@rtype: L{esys.lsm.util.Vec3<esys.lsm.util.FoundationPy>}\n"
        "@return: Returns the centre point of this sphere.\n"
      )
      .def(
        "getPos",
        &SimpleSpherePy::getPosnPy,
        "Returns the centre point of this sphere.\n"
        "@rtype: L{esys.lsm.util.Vec3<esys.lsm.util.FoundationPy>}\n"
        "@return: Returns the centre point of this sphere.\n"
      )
      .def(
        "getId",
        &SimpleSpherePy::getID,
        "Returns the id associated with this particle.\n"
        "@rtype: int\n"
        "@return: Particle ID."
      )
      .def(
        "getTag",
        &SimpleSpherePy::getTag,
        "Returns the tag assigned to this particle.\n"
        "@rtype: int\n"
        "@return: tag value assigned to this particle."
      )
      .def(
        "getMass",
        &SimpleSpherePy::getMass,
        "Returns the mass of this particle.\n"
        "@rtype: float\n"
        "@return: mass"
      )
/*      .def(
        "getVolume",
        &SimpleSpherePy::getVolume,
        "Returns the volume of this particle.\n"
        "@rtype: float\n"
        "@return: volume"
      )
*/      .def(
        "setRadius",
        &SimpleSpherePy::setRad,
        (arg("radius")),
        "Sets the radius of this sphere.\n"
        "@type radius: float\n"
        "@kwarg radius: radius."
      )
      .def(
        "setRad",
        &SimpleSpherePy::setRad,
        (arg("radius")),
        "Sets the radius of this sphere.\n"
        "@type radius: float\n"
        "@kwarg radius: radius."
      )
      .def(
        "setCenter",
        &SimpleSpherePy::setPosnPy,
        (arg("position")),
        "Sets the centre-point of this sphere.\n"
        "@type position: L{esys.lsm.util.Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
        "@kwarg position: centre point."
      )
      .def(
        "setCentre",
        &SimpleSpherePy::setPosnPy,
        (arg("position")),
        "Sets the centre-point of this sphere.\n"
        "@type position: L{esys.lsm.util.Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
        "@kwarg position: centre point."
      )
      .def(
        "setPosn",
        &SimpleSpherePy::setPosnPy,
        (arg("position")),
        "Sets the centre-point of this sphere.\n"
        "@type position: L{esys.lsm.util.Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
        "@kwarg position: centre point."
      )
      .def(
        "setPos",
        &SimpleSpherePy::setPosnPy,
        (arg("position")),
        "Sets the centre-point of this sphere.\n"
        "@type position: L{esys.lsm.util.Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
        "@kwarg position: centre point."
      )
      .def(
        "translate",
        &SimpleSpherePy::translateByPy,
        (arg("translation")),
        "Moves the centre-point of this sphere by a specified translation.\n"
        "@type translation: L{esys.lsm.util.Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
        "@kwarg translation: Vector added to the current centre-point position."
      )
      .def(
        "rotate",
        &SimpleSpherePy::rigidRotatePy,
        (
          arg("axis"),
          arg("pt")
        ),
        "Performs rigid body rotation of the sphere about a specified axis.\n"
        "@type axis: L{esys.lsm.util.Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
        "@kwarg axis: The direction and magnitude of the axis of rotation.\n"
        "@type pt: L{esys.lsm.util.Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
        "@kwarg pt: Point through which the axis passes.\n"
      )
      .def(
        "setId",
        &SimpleSpherePy::setID,
        (arg("id")),
        "Sets the ID of this spherical particle.\n"
        "@type id: int\n"
        "@kwarg id: Identifier for this particle."
      )
      .def(
        "setTag",
        &SimpleSpherePy::setTag,
        (arg("tag")),
        "Sets the tag associated with this particle.\n"
        "@type tag: int\n"
        "@kwarg tag: Identifier for this particle."
      )
      .def(
        "setMass",
        &SimpleSpherePy::setMass,
        (arg("mass")),
        "Sets the mass of this particle.\n"
        "@type mass: float\n"
        "@kwarg mass: the new mass."
      )
      .def(
        "__eq__",
        &SimpleSpherePy::operator==,
        (arg("p")),
        "Equality operator for particles,"
        " simply compares particle ID values.\n"
        "@rtype: bool\n"
        "@return: self.getId() == p.getId()"
      )
      .def(
        "__hash__",
        &SimpleSpherePy::getID,
        "Returns a hash value for use in dictionaries and sets.\n"
        "@rtype: int\n"
        "@return: self.getId()"
      )
      .def_pickle(SimpleSpherePyPickleSuite())
      ;
    }
  }
}
