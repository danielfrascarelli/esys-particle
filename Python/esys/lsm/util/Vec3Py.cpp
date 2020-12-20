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
#include <iostream>
#include <sstream>
#include "Python/esys/lsm/util/Vec3Py.h"

using namespace boost::python;
using namespace esys::lsm;

namespace esys
{
  namespace lsm
  {
    Vec3Py::Vec3Py() : Vec3()
    {
    }

    Vec3Py::Vec3Py(double x, double y, double z) : Vec3(x,y,z)
    {
    }

    Vec3Py::Vec3Py(const Vec3Py &v) : Vec3(v)
    {
    }

    Vec3Py::Vec3Py(const Vec3 &v) : Vec3(v)
    {
    }

    Vec3Py::Vec3Py(const boost::python::object &pyOb)
    {
      if (esys::lsm::bpu::len(pyOb) == 3)
      {
        *this = 
          Vec3Py(
            boost::python::extract<double>(pyOb[0]),
            boost::python::extract<double>(pyOb[1]),
            boost::python::extract<double>(pyOb[2])
          );
      }
      else
      {
        std::stringstream msg;
        msg 
          << "Could not extract (x,y,z) elements from: "
          << boost::python::extract<std::string>(boost::python::str(pyOb))();
        throw runtime_error(msg.str());
      }
    }
    
    int Vec3Py::len() const
    {
      return 3;
    }

    int Vec3Py::getIndex(int i) const
    {
      if (i < 0)
          i += len();
      if ((i >= len()) || i < 0)
      {
          PyErr_SetString(PyExc_IndexError, "Index out of range");
          throw_error_already_set();
      }
      return i;
    }
    
    double Vec3Py::getItem(int i) const
    {
      return (*this)[getIndex(i)];
    }

    void Vec3Py::setItem(int i, double val)
    {
      (*this)[getIndex(i)] = val;
    }

    Vec3Py Vec3Py::operator-(const Vec3Py &v) const
    {
        return Vec3Py(Vec3::operator-(v));
    }

    Vec3Py Vec3Py::operator+(const Vec3Py &v) const
    {
        return Vec3Py(Vec3::operator+(v));
    }

    Vec3Py Vec3Py::operator+(double s) const
    {
        return Vec3Py(Vec3::operator+(s));
    }

    Vec3Py Vec3Py::operator-(double s) const
    {
        return Vec3Py(Vec3::operator-(s));
    }

    Vec3Py Vec3Py::operator*(double s) const
    {
        return Vec3Py(Vec3::operator*(s));
    }

	Vec3Py operator * (double s, const Vec3Py& v) 
	{
		return v*s;
	}

    Vec3Py Vec3Py::operator/(double s) const
    {
        return Vec3Py(Vec3::operator/(s));
    }

    double Vec3Py::norm() const
    {
      return Vec3::norm();
    }

    Vec3Py Vec3Py::rotatePy(const Vec3Py &axis, const Vec3Py &axisPt) const
    {
      return rotate(axis, axisPt);
    }

    double Vec3Py::dot(const Vec3Py &v) const
    {
      return ::dot(*this, v);
    }

    Vec3Py Vec3Py::cross(const Vec3Py &v) const
    {
      return ::cross(*this, v);
    }

    std::string Vec3Py::toString() const
    {
      return StringUtil::toString(*this);
    }

    boost::python::list Vec3Py::toList() const
    {
      boost::python::list l;
      l.append(getItem(0));
      l.append(getItem(1));
      l.append(getItem(2));
      return l;
    }

    boost::python::tuple Vec3Py::toTuple() const
    {
      return boost::python::tuple(toList());
    }

    class Vec3PyPickleSuite : public boost::python::pickle_suite
    {
    public:
      static
      boost::python::tuple
      getinitargs(Vec3Py const& v)
      {
          return v.toTuple();
      }
    };

    using boost::python::arg;
    void exportVec3()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      class_<esys::lsm::Vec3Py>(
        "Vec3",
        "A 3D coordinate.",
        init<>()
      )
      .def(init<const Vec3Py &>())
      .def(init<const object &>())
      .def(
        init<double,double,double>(
          ( arg("x"), arg("y"), arg("z") ),
          "Constructs coordinate with specifed component values.\n"
          "@type x: float\n"
          "@kwarg x: index 0\n"
          "@type y: float\n"
          "@kwarg y: index 1\n"
          "@type z: float\n"
          "@kwarg z: index 2\n"
        )
      )
      .def(self == self)
      .def(self - self)
      .def(self + self)
      .def(self + double(0.0))
      .def(self - double(0.0))
      .def(self * double(0.0))
      .def(double() * self)
      .def(self / double(0.0))
      .def(
        "dot",
        &Vec3Py::dot,
        ( arg("v") ),
        "Returns the dot product of this 3 element vector with\n"
        "the specified L{Vec3}.\n"
        "@type v: L{Vec3}\n"
        "@kwarg v: dot product with this\n"
        "@rtype: float"
      )
      .def(
        "rotate",
        &Vec3Py::rotatePy,
        ( arg("axis"), arg("axisPt") ),
        "Returns a point which is this point rotated about a specified axis.\n"
        "@type axis: L{Vec3}\n"
        "@kwarg axis: axis of rotation and counter-clockwise angle of rotation"
        " C{= axis.norm()} radians.\n"
        "@type axisPt: L{Vec3}\n"
        "@kwarg axisPt: The C{axis} of rotation is assumed to pass through"
        " this point.\n"
        "@rtype: L{Vec3}\n"
        "@return: The point which is the rotation of this (self) coordinate"
        " about the specified axis."
      )
      .def(
        "cross",
        &Vec3Py::cross,
        ( arg("v") ),
        "Returns the cross product of this 3 element vector with\n"
        "the specified L{Vec3}.\n"
        "@type v: L{Vec3}\n"
        "@kwarg v: cross product with this\n"
        "@rtype: L{Vec3}"
      )
      .def(
        "norm",
        &Vec3Py::norm,
        "Returns the magnitude of this 3 element vector.\n"
        "@rtype: float\n"
        "@return: math.sqrt(self.dot(self))."
      )
      .def(self_ns::str(self))
      .def("__len__", &Vec3Py::len)
      .def("__getitem__", &Vec3Py::getItem)
      .def("__setitem__", &Vec3Py::setItem)
      .def(
        "toList",
        &Vec3Py::toList,
        "Returns a python list of 3 elements.\n"
        "@rtype: list of three floats\n"
        "@return: [self[0],self[1],self[2]]"
      )
      .def(
        "toTuple",
        &Vec3Py::toTuple,
        "Returns a python tuple of 3 elements.\n"
        "@rtype: tuple of three floats\n"
        "@return: (self[0],self[1],self[2])"
      )
      .def_pickle(Vec3PyPickleSuite())
      ;
    }
  }
}

std::ostream &operator<<(std::ostream &oStream, const esys::lsm::Vec3Py &vec)
{
  oStream << vec[0] << " " << vec[1] << " " << vec[2];
  return oStream;
}
