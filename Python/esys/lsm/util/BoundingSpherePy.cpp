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
#include <iostream>
#include <sstream>
#include "Python/esys/lsm/util/BoundingSpherePy.h"

using namespace boost::python;
using namespace esys::lsm;

namespace esys
{
  namespace lsm
  {
    BoundingSpherePy::BoundingSpherePy() : BoundingSphere()
    {
    }

    BoundingSpherePy::BoundingSpherePy(const Vec3Py &centrePt, double radius)
      : BoundingSphere(centrePt, radius)
    {
    }

    BoundingSpherePy::BoundingSpherePy(const Vec3 &centrePt, double radius)
      : BoundingSphere(centrePt, radius)
    {
    }

    BoundingSpherePy::BoundingSpherePy(const BoundingSpherePy &v)
      : BoundingSphere(v)
    {
    }

    BoundingSpherePy::BoundingSpherePy(const BoundingSphere &v)
      : BoundingSphere(v)
    {
    }

    BoundingSpherePy::BoundingSpherePy(const boost::python::object &pyCentre, double radius)
    {
      if (esys::lsm::bpu::len(pyCentre) == 3)
      {
        *this = BoundingSpherePy(Vec3Py(pyCentre), radius);
      }
      else
      {
        const std::string objectString =
          boost::python::extract<std::string>(boost::python::str(pyCentre))();
        std::stringstream msg;
        msg << "Could not extract (x,y,z) elements from: " << objectString;
        throw runtime_error(msg.str());
      }
    }

    bool BoundingSpherePy::operator==(const BoundingSpherePy &bSphere) const
    {
      return BoundingSphere::operator ==(bSphere);
    }

    Vec3Py BoundingSpherePy::getCentrePy() const
    {
      return getCentre();
    }

    class BoundingSpherePyPickleSuite : public boost::python::pickle_suite
    {
    public:
      static
      boost::python::tuple
      getinitargs(BoundingSpherePy const& bSphere)
      {
          return boost::python::make_tuple(bSphere.getCentrePy(), bSphere.getRadius());
      }
    };

    using boost::python::arg;
    void exportBoundingSphere()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      class_<esys::lsm::BoundingSpherePy>(
        "BoundingSphere",
        "A sphere - centre point and radius.",
        init<>()
      )
        .def(
          init<const BoundingSpherePy &>(
            ( arg("bSphere") )
          )
        )
        .def(
          init<const object &, double>(
            ( arg("seq"), arg("radius") )
          )
        )
        .def(
          init<const Vec3Py &, double>(
            ( arg("centrePt")=Vec3Py(0,0,0), arg("radius")=0.0 ),
            "Creates sphere with specified centre and radius.\n"
            "@type centrePt: L{Vec3}\n"
            "@kwarg centrePt: Centre of sphere.\n"
            "@type radius: float\n"
            "@kwarg radius: Radius of sphere."
          )
        )
        .def(
          "getCentre",
          &BoundingSpherePy::getCentrePy,
          "Returns centre point of this sphere.\n"
          "@rtype: L{Vec3}"
        )
        .def(
          "getCenter",
          &BoundingSpherePy::getCentrePy,
          "Returns centre point of this sphere.\n"
          "@rtype: L{Vec3}"
        )
        .def(
          "getRadius",
          &BoundingSpherePy::getRadius,
          "Returns radius of this sphere.\n"
          "@rtype: float"
        )
        .def(
          self == self//, "Tests equality of spheres - equality of centre point and radius."
        )
        .def(self_ns::str(self))
        .def_pickle(BoundingSpherePyPickleSuite())
      ;
    }
  }
}

std::ostream &operator<<(std::ostream &oStream, const esys::lsm::BoundingSpherePy &bSphere)
{
  oStream
    << bSphere.getCentre()[0]
    << " " << bSphere.getCentre()[1]
    << " " << bSphere.getCentre()[2]
    << " " << bSphere.getRadius();
  return oStream;
}
