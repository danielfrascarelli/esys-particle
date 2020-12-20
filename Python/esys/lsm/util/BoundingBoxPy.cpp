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
#include "Python/esys/lsm/util/BoundingBoxPy.h"

using namespace boost::python;
using namespace esys::lsm;

namespace esys
{
  namespace lsm
  {
    BoundingBoxPy::BoundingBoxPy() : BoundingBox()
    {
    }

    BoundingBoxPy::BoundingBoxPy(const Vec3Py &minPt, const Vec3Py &maxPt) : BoundingBox(minPt, maxPt)
    {
    }

    BoundingBoxPy::BoundingBoxPy(const Vec3 &minPt, const Vec3 &maxPt) : BoundingBox(minPt, maxPt)
    {
    }

    BoundingBoxPy::BoundingBoxPy(const BoundingBoxPy &v) : BoundingBox(v)
    {
    }

    BoundingBoxPy::BoundingBoxPy(const BoundingBox &v) : BoundingBox(v)
    {
    }

    BoundingBoxPy::BoundingBoxPy(const boost::python::object &pyMin, const boost::python::object &pyMax)
    {
      if ((esys::lsm::bpu::len(pyMin) == 3) && (esys::lsm::bpu::len(pyMax) == 3))
      {
        *this = BoundingBoxPy(Vec3Py(pyMin), Vec3Py(pyMax));
      }
      else
      {
        std::string objectString;
        if (esys::lsm::bpu::len(pyMin) != 3)
        {
          objectString = boost::python::extract<std::string>(boost::python::str(pyMin))();
        }
        else
        {
          objectString = boost::python::extract<std::string>(boost::python::str(pyMax))();
        }
        std::stringstream msg;
        msg << "Could not extract (x,y,z) elements from: " << objectString;
        throw runtime_error(msg.str());
      }
    }

    bool BoundingBoxPy::intersectsWithVec3Py(const Vec3Py &pt) const
    {
      return this->contains(pt);
    }

    bool BoundingBoxPy::operator==(const BoundingBoxPy &bBox) const
    {
      return
        (
          (getMinPt() == bBox.getMinPt())
          &&
          (getMaxPt() == bBox.getMaxPt())
        );
    }

    Vec3Py BoundingBoxPy::getMinPtPy() const
    {
      return Vec3Py(BoundingBox::getMinPt());
    }

    Vec3Py BoundingBoxPy::getMaxPtPy() const
    {
      return Vec3Py(BoundingBox::getMaxPt());
    }

    Vec3Py BoundingBoxPy::getSizePy() const
    {
      return Vec3Py(getMaxPt()-getMinPt());
    }

    Vec3Py BoundingBoxPy::getCentrePy() const
    {
      return Vec3Py((getMaxPt()+getMinPt())*0.5);
    }

    class BoundingBoxPyPickleSuite : public boost::python::pickle_suite
    {
    public:
      static
      boost::python::tuple
      getinitargs(BoundingBoxPy const& bBox)
      {
          return boost::python::make_tuple(bBox.getMinPtPy(), bBox.getMaxPtPy());
      }
    };

    using boost::python::arg;
    void exportBoundingBox()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      class_<esys::lsm::BoundingBoxPy>(
        "BoundingBox",
        "An axis aligned box.",
        init<>()
      )
      .def(init<const BoundingBoxPy &>())
      .def(init<const object &, const object &>())
      .def(
        init<const Vec3Py &, const Vec3Py &>(
          ( arg("minPt"), arg("maxPt") ),
          "Construct box by specifying lower left corner and\n"
          "upper right corner coordinates.\n"
          "@type minPt: L{Vec3}\n"
          "@kwarg minPt: lower left corner coordinate.\n"
          "@type maxPt: L{Vec3}\n"
          "@kwarg maxPt: upper right corner coordinate."
        )
      )
      .def(
        "getMinPt",
        &BoundingBoxPy::getMinPtPy,
        "Returns lower left corner coordinate.\n"
        "@rtype: L{Vec3}"
      )
      .def(
        "getMaxPt",
        &BoundingBoxPy::getMaxPtPy,
        "Returns upper right corner coordinate.\n"
        "@rtype: L{Vec3}"
      )
      .def(
        "getSize",
        &BoundingBoxPy::getSizePy,
        "Returns side-lengths for each coordinate component.\n"
        "@rtype: L{Vec3}"
      )
      .def(
        "getCentre",
        &BoundingBoxPy::getCentrePy,
        "Returns the centre point of this box.\n"
        "@rtype: L{Vec3}"
      )
      .def(
        "getCenter",
        &BoundingBoxPy::getCentrePy,
        "Returns the centre point of this box.\n"
        "@rtype: L{Vec3}"
      )
      .def(
        "intersectsWith",
        &BoundingBoxPy::intersectsWithVec3Py,
        "Returns C{True} if the specified point intersects"
        " with this box.\n"
        "@type pt: L{esys.lsm.util.Vec3<esys.lsm.util.FoundationPy.Vec3>}\n"
        "@kwarg pt: A point.\n"
        "@rtype: L{bool}"
      )
      .def(self == self)
      .def(self_ns::str(self))
      .def_pickle(BoundingBoxPyPickleSuite())
      ;
    }
  }
}

std::ostream &operator<<(std::ostream &oStream, const esys::lsm::BoundingBoxPy &bBox)
{
  oStream <<        bBox.getMinPt()[0] << " " << bBox.getMinPt()[1] << " " << bBox.getMinPt()[2];
  oStream << " " << bBox.getMaxPt()[0] << " " << bBox.getMaxPt()[1] << " " << bBox.getMaxPt()[2];
  return oStream;
}
