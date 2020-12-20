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
#include "Parallel/LatticeMaster.h"
#include "Python/esys/lsm/sim/coupled/LsmMpiEscriptPy.h"
#include "Python/esys/lsm/LsmMpiPy.h"
#include "Python/esys/lsm/util/Vec3Py.h"

namespace esys
{
  namespace lsm
  {
    class LsmMpiEscriptPy : public LsmMpiPy
    {
    public:
      LsmMpiEscriptPy(int numWorkers, const boost::python::list &dimList)
        : LsmMpiPy(numWorkers, dimList)
      {
      }
    };

    void exportLsmMpiEscript()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<LsmMpiEscriptPy, boost::python::bases<LsmMpiPy> >(
        "LsmMpiEscript",
        boost::python::init<int, const boost::python::list &>()
      );
    }
  }
}
