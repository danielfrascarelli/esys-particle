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
#include "Python/esys/lsm/util/QuaternionPy.h"
#include "Python/esys/lsm/util/Vec3Py.h"
#include "Foundation/StringUtil.h"
#include "Python/BoostPythonUtil/Util.h"

using namespace boost::python;
using namespace esys::lsm;

namespace esys
{
  namespace lsm
  {
    QuaternionPy::QuaternionPy() : Quaternion()
    {
    }

    QuaternionPy::QuaternionPy(
      double q0,
      double q1,
      double q2,
      double q3
    ) : Quaternion(q0, Vec3(q1, q2, q3))
    {
    }

    QuaternionPy::QuaternionPy(const Vec3Py &v)
      : Quaternion(
          cos(v.norm()/2.0),
          v*(sin(v.norm()/2.0)/v.norm())
        )
    {
    }

    QuaternionPy::QuaternionPy(const QuaternionPy &q) : Quaternion(q)
    {
    }

    QuaternionPy::QuaternionPy(const Quaternion &q) : Quaternion(q)
    {
    }

    QuaternionPy::QuaternionPy(const boost::python::object &pyOb)
    {
      if (esys::lsm::bpu::len(pyOb) == 4)
      {
        *this = 
          QuaternionPy(
            boost::python::extract<double>(pyOb[0]),
            boost::python::extract<double>(pyOb[1]),
            boost::python::extract<double>(pyOb[2]),
            boost::python::extract<double>(pyOb[3])
          );
      }
      else
      {
        std::stringstream msg;
        msg 
          << "Could not extract (q0,q1,q2,q3) elements from: "
          << boost::python::extract<std::string>(boost::python::str(pyOb))();
        throw runtime_error(msg.str());
      }
    }
    
    int QuaternionPy::len() const
    {
      return 4;
    }

    int QuaternionPy::getIndex(int i) const
    {
      const int origI = i;
      if (i < 0)
      {
          i += len();
      }
      if ((i >= len()) || i < 0)
      {
        std::stringstream msg;
        msg << "Index " << origI << " out of range [0,4)";
        PyErr_SetString(PyExc_IndexError, msg.str().c_str());
        throw_error_already_set();
      }
      return i;
    }
    
    double QuaternionPy::getItem(int i) const
    {
      return ((i=getIndex(i)) == 0) ? return_sca() : return_vec()[i-1];
    }

    void QuaternionPy::setItem(int i, double val)
    {
      i = getIndex(i);
      if (i == 0)
      {
        set_scalar(val);
      }
      else
      {
        Vec3 v = return_vec();
        v[i-1] = val;
        set_vector(v);
      }
    }

    std::string QuaternionPy::toString() const
    {
      return StringUtil::toString(*this);
    }

    boost::python::list QuaternionPy::toList() const
    {
      boost::python::list l;
      l.append(getItem(0));
      l.append(getItem(1));
      l.append(getItem(2));
      l.append(getItem(3));
      return l;
    }

    boost::python::tuple QuaternionPy::toTuple() const
    {
      return boost::python::tuple(toList());
    }

    boost::python::tuple
    QuaternionPy::PickleSuite::getinitargs(QuaternionPy const& q)
    {
      return q.toTuple();
    }

    Vec3Py QuaternionPy::asAngleAxis() const
    {
      return Vec3Py(Quaternion::asAngleAxis());
    }

    boost::python::tuple QuaternionPy::asAngleAxisPair() const
    {
      Quaternion::AngleAxisPair p = Quaternion::asAngleAxisPair();
      return boost::python::make_tuple(p.first, Vec3Py(p.second));
    }

    using boost::python::arg;
    void exportQuaternion()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      class_<esys::lsm::QuaternionPy>(
        "Quaternion",
        "A quaternion.",
        init<>()
      )
      .def(init<const object &>())
      .def(init<const Vec3Py &>())
      .def(init<const QuaternionPy &>())
      .def(
        init<double,double,double,double>(
          (arg("q0"), arg("q1"), arg("q2"), arg("q3") ),
          "Constructs quaternion with specifed component values.\n"
          "@type q0: float\n"
          "@kwarg q0: Scalar part, index 0.\n"
          "@type q1: float\n"
          "@kwarg q1: index 1\n"
          "@type q2: float\n"
          "@kwarg q2: index 2\n"
          "@type q3: float\n"
          "@kwarg q3: index 3\n"
        )
      )
      .def(
        "normalise",
        &QuaternionPy::normalize,
        "Normalises this quaternion.\n"
      )
      .def(
        "normalize",
        &QuaternionPy::normalize,
        "Normalizes this quaternion.\n"
      )
      .def(
        "asAngleAxis",
        &QuaternionPy::asAngleAxis,
        "Returns angle and axis rotation representation, where"
        " angle is the magnitude of the returned "
        " L{Vec3<esys.lsm.util.Vec3>} object.\n"
        "@rtype: L{Vec3<esys.lsm.util.Vec3>}\n"
        "@return: angle and axis rotation representation."
      )
      .def(
        "asAngleAxisPair",
        &QuaternionPy::asAngleAxis,
        "Returns angle and axis rotation representation as a"
        " two element tuple, the first element is the"
        " angle and the second element is the"
        " L{Vec3<esys.lsm.util.Vec3>} axis.\n"
        "@rtype: (float, L{Vec3<esys.lsm.util.Vec3>})\n"
        "@return: angle and axis rotation tuple-pair."
      )
      .def(self_ns::str(self))
      .def(boost::python::self == boost::python::self)
      .def("__len__", &QuaternionPy::len)
      .def("__getitem__", &QuaternionPy::getItem)
      .def("__setitem__", &QuaternionPy::setItem)
      .def(
        "toList",
        &QuaternionPy::toList,
        "Returns a python list of 4 elements.\n"
        "@rtype: list of four floats\n"
        "@return: C{[self[0],self[1],self[2],self[3]]}"
      )
      .def(
        "toTuple",
        &QuaternionPy::toTuple,
        "Returns a python tuple of 4 elements.\n"
        "@rtype: tuple of four floats\n"
        "@return: C{(self[0],self[1],self[2],self[3])}"
      )
      .def_pickle(QuaternionPy::PickleSuite())
      ;
    }
  }
}

std::ostream &operator<<(std::ostream &oStream, const esys::lsm::QuaternionPy &quat)
{
  oStream
    << quat.return_sca() << " "
    << quat.return_vec()[0] << " "
    << quat.return_vec()[1] << " "
    << quat.return_vec()[2];
  return oStream;
}
