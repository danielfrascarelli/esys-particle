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
#include "Python/esys/lsm/RunnablePy.h"

namespace esys
{
  namespace lsm
  {
    void RunnablePy::run()
    {
      boost::python::override overrideMethod = this->get_override("run");
      if (overrideMethod)
      {
        overrideMethod();
      }
      else
      {
        defaultRun();
      }
    }

    void RunnablePy::defaultRun()
    {
      PyErr_SetString(
        PyExc_NotImplementedError,
        "run method needs to be implemented in derived class"
      );
      boost::python::throw_error_already_set();
    }

    void exportRunnable()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<
        RunnablePy,
        boost::noncopyable
      >(
        "Runnable",
        "Base class with overridable C{run} method, that is subclassed to implement user-defined subroutines called each timestep of a simulation.\n"
      )
      .def(
        "run",
        &Runnable::run,
        "Method called each timestep, overridden by subclasses of the Runnable base class.",
        &RunnablePy::defaultRun
      )
      ;
    }
  }
}
