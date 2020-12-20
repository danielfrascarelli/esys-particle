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
#include "Python/esys/lsm/MeshBuildParamsPy.h"

namespace esys
{
  namespace lsm
  {
    /*! 
      constructor for MeshTagBuildPrmsPy
    
      \param tag
      \param mask
    */
    MeshTagBuildPrmsPy::MeshTagBuildPrmsPy(int tag, int mask) : MeshTagBuildPrms(tag, mask)
    {
    }
    
    /*!
      constructor for MeshGapBuildPrmsPy
    
      \param maxGap
    */
    MeshGapBuildPrmsPy::MeshGapBuildPrmsPy(double maxGap) : MeshGapBuildPrms(maxGap)
    {
    }

    using boost::python::arg;
    void exportMeshBuildPrms()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<MeshTagBuildPrmsPy>(
        "MeshTagBuildPrms",
        "Parameters for bonding particles to a mesh which uses\n"
        "particle-tags as the criterion for creating bonds.",
        boost::python::init<int, int>(
          (
            arg("tag"),
            arg("mask")
          ),
	  "Build parameters for particles with a specified tag\n"
          "@type tag: int\n"
          "@kwarg tag: particle tag, particles with this tag will be\n"
          " bonded to the specified mesh.\n"
          "@type mask: int\n"
          "@kwarg mask: mask applied (anded) to particle tag before"
          " being compared with C{tag}.\n"
        )
      );

      boost::python::class_<MeshGapBuildPrmsPy>(
        "MeshGapBuildPrms",
        "Parameters for bonding particles to a mesh which uses\n"
        "distance-to-mesh criterion for creating bonds.",
        boost::python::init<double>(
          (arg("maxDistance")),
	  "Build parameters for particles within a specified distance of a mesh wall\n"
          "@type maxDistance: float\n"
          "@kwarg maxDistance: Particles which are closer than C{maxDistance}"
          " to the mesh get bonded to the mesh.\n"
        )
      );
    }
  }
}
