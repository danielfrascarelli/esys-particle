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
#include <boost/version.hpp>
#include <boost/python.hpp>
#include "Foundation/Rng.h"

namespace esys
{
  namespace lsm
  {
    void seedDefaultRng(unsigned int seed)
    {
      rng::s_zeroOneUniform.seed(seed);
    }
    
    using boost::python::arg;
    void exportRng()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);
      boost::python::docstring_options doc_str_opt(true,false);

      boost::python::def(
        "seedDefaultRng",
        &seedDefaultRng,
        ( arg("seed")=1234567891 ),
        "Seeds the default C{[0,1]} uniform pseudo-random number generator.\n"
        "@type seed: int\n"
        "@kwarg seed: unsigned integer seed value."
      );
    }
  }
}
