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
#include "Python/esys/lsm/InteractionParamsPy.h"
#include "Python/esys/lsm/SphereBodyPrmsPy.h"

namespace esys
{
  namespace lsm
  {

    /*!
      constructor of non-rotational elastic wall parameters
      
      \param name the name of the interaction
      \param wname the name of the wall
      \param normalSpringK the spring constant for the elastic interactions
    */
    NRotElasticSphereBodyPrmsPy::NRotElasticSphereBodyPrmsPy(
      const std::string& name,
      const std::string& wname,
      double normalSpringK
    )
      : CESphereBodyIGP(name,wname,normalSpringK)
    {
    }

    using boost::python::arg;
    void exportSphereBodyPrms()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);


      boost::python::class_<NRotElasticSphereBodyPrmsPy,boost::python::bases<InteractionPrmsPy> >(
        "NRotElasticSphereBodyPrms",
        "Parameters for linear elastic contact between particles and"
        " a sphere body.",
        boost::python::init<
          const std::string&,
          const std::string&,
          double
        >(
          (
            arg("name"),
            arg("sphereName"),
            arg("normalK")
          ),
	  "Parameters defining elastic contacts between particles and a planar wall\n"
          "@type name: string\n"
          "@kwarg name: Name assigned to the created interaction group.\n"
          "@type sphereName: string\n"
          "@kwarg sphereName: The name of an existing sphere body.\n"
          "@type normalK: float\n"
          "@kwarg normalK: spring constant for the linear elastic force"
          " calculation.\n"
        )
      )
      ;

      
    }
  }
}

