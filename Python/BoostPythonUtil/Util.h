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


#ifndef ESYS_LSM_BPUUTIL_H
#define ESYS_LSM_BPUUTIL_H

#include <boost/python.hpp>

namespace esys
{
  namespace lsm
  {
    namespace bpu
    {
      inline boost::python::object iter(const boost::python::object &pyOb)
      {
        return pyOb.attr("__iter__")();
      }

      inline int len(const boost::python::object &pyOb)
      {
        return boost::python::extract<int>(pyOb.attr("__len__")());
      }
    }
  }
}

#endif
