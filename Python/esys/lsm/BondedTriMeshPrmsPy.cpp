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
#include "Python/esys/lsm/BondedTriMeshPrmsPy.h"
#include "Python/esys/lsm/MeshBuildParamsPy.h"

namespace esys
{
  namespace lsm
  {
    NRotBondedTriMeshPrmsPy::NRotBondedTriMeshPrmsPy(
      const string &interactionName,
      const string &meshName,
      double normalK,
      double breakDistance,
      const MeshTagBuildPrmsPy &buildPrms)
      : BTriMeshIP(interactionName, meshName, normalK, breakDistance),
        m_tagPrmsPtr(TagBuildPrmsPtr(new MeshTagBuildPrmsPy(buildPrms))),
        m_gapPrmsPtr()
    {
    }
    
    NRotBondedTriMeshPrmsPy::NRotBondedTriMeshPrmsPy(
      const string &interactionName,
      const string &meshName,
      double normalK,
      double breakDistance,
      const MeshGapBuildPrmsPy &buildPrms
    )
      : BTriMeshIP(interactionName, meshName, normalK, breakDistance),
        m_tagPrmsPtr(),
        m_gapPrmsPtr(GapBuildPrmsPtr(new MeshGapBuildPrmsPy(buildPrms)))
    {
    }

    using boost::python::arg;
    void exportBondedTriMeshPrms()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<NRotBondedTriMeshPrmsPy,boost::python::bases<InteractionPrmsPy> >(
        "NRotBondedTriMeshPrms",
        "Parameters for specifying linear elastic bonds between particles\n"
        "and a triangular mesh surface.",
        boost::python::init<
          const std::string &,
          const std::string &,
          double,
          double,
          const MeshTagBuildPrmsPy &
        >(
          (
            arg("name"),
            arg("meshName"),
            arg("normalK"),
            arg("breakDistance"),
            arg("buildPrms")
          )
        )
      ).def(
        boost::python::init<
          const std::string &,
          const std::string &,
          double,
          double,
          const MeshGapBuildPrmsPy &
        >(
          (
            arg("name"),
            arg("meshName"),
            arg("normalK"),
            arg("breakDistance"),
            arg("buildPrms")
          ),
          "Parameters defining bonded elastic interactions between particles and triangle mesh walls\n"
          "@type name: string\n"
          "@kwarg name: name of interaction group.\n"
          "@type meshName: string\n"
          "@kwarg meshName: name of the triangular mesh for which elastic"
          " bonds will be created.\n"
          "@type normalK: float\n"
          "@kwarg normalK: spring constant for linear elastic bond force"
          " calculation.\n"
          "@type breakDistance: float\n"
          "@kwarg breakDistance: When distance between mesh and particle exceeds"
          " this distance, the bond breaks.\n"
          "@type buildPrms: L{MeshTagBuildPrms} or L{MeshGapBuildPrms}\n"
          "@kwarg buildPrms: Object which specifies the method of bond creation."
        )
      )
      .def(
        "getName",
        &NRotBondedTriMeshPrmsPy::getName,
        boost::python::return_value_policy<boost::python::copy_const_reference>(),
        "@rtype: string\n"
        "@return: Name of this interaction group."
      )
      .def(
        "getMeshName",
        &NRotBondedTriMeshPrmsPy::getMeshName,
        boost::python::return_value_policy<boost::python::copy_const_reference>(),
        "@rtype: string\n"
        "@return: name of the mesh for which the bonded interactions apply."
      )
     ;
    }
  }
}
