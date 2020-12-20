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


#ifndef ESYS_LSMVEC3PY_H
#define ESYS_LSMVEC3PY_H

#include <boost/python.hpp>
#include "Foundation/vec3.h"
#include "Foundation/StringUtil.h"
#include "Python/BoostPythonUtil/Util.h"

#include <sstream>

namespace esys
{
  namespace lsm
  {
    class Vec3Py : public Vec3
    {
    public:
      Vec3Py();

      Vec3Py(double x, double y, double z);

      Vec3Py(const Vec3Py &v);

      Vec3Py(const Vec3 &v);

      Vec3Py(const boost::python::object &pyOb);

      int len() const;

      double getItem(int i) const;

      void setItem(int i, double val);

      Vec3Py operator-(const Vec3Py &v) const;
      Vec3Py operator+(const Vec3Py &v) const;

      Vec3Py operator+(double s) const;
      Vec3Py operator-(double s) const;
      Vec3Py operator*(double s) const;
      Vec3Py operator/(double s) const;
	  friend Vec3Py operator * (double, const Vec3Py&);

      Vec3Py rotatePy(const Vec3Py &axis, const Vec3Py &axisPt) const;
      
      double norm() const;

      double dot(const Vec3Py &v) const;

      Vec3Py cross(const Vec3Py &v) const;

      std::string toString() const;

      boost::python::list toList() const;

      boost::python::tuple toTuple() const;

      int getIndex(int i) const;
    };

    void exportVec3();    
  }
}

std::ostream &operator<<(std::ostream &oStream, const esys::lsm::Vec3Py &vec);

#endif
