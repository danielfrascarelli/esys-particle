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

#ifndef ESYS_LSMQUATERNIONPY_H
#define ESYS_LSMQUATERNIONPY_H

#include <boost/python.hpp>
#include "Foundation/Quaternion.h"

#include <sstream>

namespace esys
{
  namespace lsm
  {
    class Vec3Py;
    
    class QuaternionPy : public Quaternion
    {
    public:
      QuaternionPy();

      QuaternionPy(double q0, double q1, double q2, double q3);

      QuaternionPy(const Vec3Py &v);
      
      QuaternionPy(const QuaternionPy &q);

      QuaternionPy(const Quaternion &q);

      QuaternionPy(const boost::python::object &pyOb);

      int len() const;

      double getItem(int i) const;

      void setItem(int i, double val);

      double norm() const;

      Vec3Py asAngleAxis() const;

      boost::python::tuple asAngleAxisPair() const;

      std::string toString() const;

      boost::python::list toList() const;

      boost::python::tuple toTuple() const;

      int getIndex(int i) const;

      class PickleSuite : public boost::python::pickle_suite
      {
      public:
        static
        boost::python::tuple
        getinitargs(QuaternionPy const& q);
      };
    };

    void exportQuaternion();
  }
}

std::ostream &operator<<(std::ostream &oStream, const esys::lsm::QuaternionPy &vec);

#endif //ESYS_LSMQUATERNIONPY_H
