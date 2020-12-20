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

#include <mpi.h>
#include <boost/version.hpp>
#include <boost/python.hpp>
#include "Python/esys/lsm/InteractionGroupPy.h"
#include "Python/esys/lsm/LsmMpiPy.h"

namespace esys
{
  namespace lsm
  {
    InteractionGroupPy::InteractionGroupPy(
      LsmMpiPy &lsm,
      const std::string &name
    )
      : m_pLsm(&lsm), m_name(name)
    {
    }

    using boost::python::arg;
    
    void exportInteractionGroup()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<InteractionGroupPy>(
        "InteractionGroup",
        "Base class for interaction groups.",
        boost::python::no_init
      )
      .def(
        "getName",
        &InteractionGroupPy::getName,
        boost::python::return_value_policy<
          boost::python::copy_const_reference
        >(),
        "Returns the name of this interaction group.\n"
        "@rtype: str\n"
        "@return: Name of this interaction group."
      )
      ;
    }
  }
}
