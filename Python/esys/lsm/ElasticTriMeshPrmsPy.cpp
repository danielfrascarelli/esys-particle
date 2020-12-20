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
#include "Python/esys/lsm/ElasticTriMeshPrmsPy.h"

namespace esys
{
  namespace lsm
  {
    NRotElasticTriMeshPrmsPy::NRotElasticTriMeshPrmsPy(
      const string &interactionName,
      const string &meshName,
      double normalK
    )
      : ETriMeshIP(interactionName, meshName, normalK)
    {
    }

    using boost::python::arg;
    void exportElasticTriMeshPrms()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<NRotElasticTriMeshPrmsPy>(
        "NRotElasticTriMeshPrms",
        "Defines linear elastic contact interaction between particles and a"
        " triangular mesh surface.",
        boost::python::init<
          const std::string &,
          const std::string &,
          double
        >(
            (
              arg("name"),
              arg("meshName"),
              arg("normalK")
            ),
	    "Parameters defining elastic contacts between particles and a triangle mesh wall\n"
            "@type name: string\n"
            "@kwarg name: Name assigned to the group of interactions.\n"
            "@type meshName: string\n"
            "@kwarg meshName: Name of the mesh for which the interactions apply.\n"
            "@type normalK: float\n"
            "@kwarg normalK: spring constant for the linear elastic contact interaction."
        )
      );
    }
  }
}
