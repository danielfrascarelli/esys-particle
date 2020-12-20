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
#include "Python/esys/lsm/WallPrmsPy.h"

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
    NRotElasticWallPrmsPy::NRotElasticWallPrmsPy(
      const std::string& name,
      const std::string& wname,
      double normalSpringK
    )
      : CEWallIGP(name,wname,normalSpringK)
    {
    }

    /*!
      constructor of non-rotational bonded wall parameters
   
      \param name the name of the interaction
      \param wname  the name of the wall
      \param normalSpringK the spring constant for the elastic interactions
      \param particleTag the tag of the particles to which the wall is bonded
      (if build via bond and not via distance)
    */
    NRotBondedWallPrmsPy::NRotBondedWallPrmsPy(
      const std::string& name,
      const std::string& wname,
      double normalSpringK,
      int particleTag
    )
      : CBWallIGP(name,wname,normalSpringK,particleTag,-1)
    {
    }

    /*!
      constructor of non-rotational bonded wall parameters
   
      \param name the name of the interaction
      \param wname  the name of the wall
      \param normalSpringK the spring constant for the elastic interactions
      \param particleTag the tag of the particles to which the wall is bonded
      \param tagMask the particle tag mask (which bits of the tag are significant)
    */
    NRotBondedWallPrmsPy::NRotBondedWallPrmsPy(
      const std::string& name,
      const std::string& wname,
      double normalSpringK,
      int particleTag,
      int tagMask
    )
      : CBWallIGP(name,wname,normalSpringK,particleTag,tagMask)
    {
    }

    /*!
      constructor of non-rotational bonded wall parameters with directional stiffness
   
      \param name the name of the interaction
      \param wname  the name of the wall
      \param SpringKx the spring constant for the elastic interactions
      \param SpringKy the spring constant for the elastic interactions
      \param SpringKz the spring constant for the elastic interactions
      \param particleTag the tag of the particles to which the wall is bonded
      (if build via bond and not via distance)
    */
    NRotSoftBondedWallPrmsPy::NRotSoftBondedWallPrmsPy(
      const std::string& name,
      const std::string& wname,
      double normalK,
      double shearK,
      int particleTag,
      int tagMask,
      bool scaling
    ) : CSoftBWallIGP(name,wname,normalK,shearK,particleTag,tagMask,scaling)
    {}

    /*!
      constructor of non-rotational bonded wall parameters with directional stiffness
      with scaling defaulting to true
   
      \param name the name of the interaction
      \param wname  the name of the wall
      \param SpringKx the spring constant for the elastic interactions
      \param SpringKy the spring constant for the elastic interactions
      \param SpringKz the spring constant for the elastic interactions
      \param particleTag the tag of the particles to which the wall is bonded
      (if build via bond and not via distance)
    */
    NRotSoftBondedWallPrmsPy::NRotSoftBondedWallPrmsPy(
      const std::string& name,
      const std::string& wname,
      double normalK,
      double shearK,
      int particleTag,
      int tagMask
    ) : CSoftBWallIGP(name,wname,normalK,shearK,particleTag,tagMask,true)
    {}

    using boost::python::arg;
    void exportWallPrms()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<NRotBondedWallPrmsPy,boost::python::bases<InteractionPrmsPy> >(
        "NRotBondedWallPrms",
        "Parameters for linear elastic bond between particles and"
        " a wall.",
	      boost::python::init<const std::string&, const std::string&, double, int, int>(
          (
            arg("name"),
            arg("wallName"),
            arg("normalK"),
            arg("particleTag"),
            arg("tagMask")
          ),
	  "Parameters defining bonded elastic interactions between particles and a planar wall\n"
          "@type name: string\n"
          "@kwarg name: Name assigned to the created interaction group.\n"
          "@type wallName: string\n"
          "@kwarg wallName: The name of an existing wall.\n"
          "@type normalK: float\n"
          "@kwarg normalK: spring constant for the linear elastic force"
          " calculation.\n"
          "@type particleTag: int\n"
          "@kwarg particleTag: Particles with this tag are bonded to the wall.\n"
          "@type tagMask: int\n"
          "@kwarg tagMask: the tag mask (default: -1) shows the significant"
          " bits of the tag.\n"
				)
			)
			.def(boost::python::init<const std::string&, const std::string&, double, int>(
        (
          arg("name"),
          arg("wallName"),
          arg("normalK"),
          arg("particleTag")
        )
			));
      
   

      boost::python::class_<NRotSoftBondedWallPrmsPy,boost::python::bases<InteractionPrmsPy> >(
        "NRotSoftBondedWallPrms",
        "Parameters for soft linear elastic bonds between particles and"
        " a wall. Soft bonds may have differing normal and shear stiffnesses.",
        boost::python::init<const std::string&, const std::string&, double, double,int,int,bool>(
          (
            arg("name"),
            arg("wallName"),
            arg("normalK"),
            arg("shearK"),
            arg("particleTag"),
            arg("tagMask"),
            arg("scaling")
          ),
	  "Parameters defining soft bonded interactions between particles and a planar wall\n"
          "@type name: string\n"
          "@kwarg name: Name assigned to the created interaction group.\n"
          "@type wallName: string\n"
          "@kwarg wallName: The name of an existing wall.\n"
          "@type normalK: float\n"
          "@kwarg normalK: spring constant for the normal elastic forces.\n"
          "@type shearK: float\n"
          "@kwarg shearK: spring constant for the shear elastic forces.\n"
          "@type particleTag: int\n"
          "@kwarg particleTag: Particles with this tag are bonded to the wall.\n"
          "@type tagMask: int\n"
          "@kwarg tagMask: the tag mask.\n"
          "@type scaling: bool\n"
          "@kwarg scaling: toggle whether to scale elastic constants by radius.\n"
        )
        )
        .def(boost::python::init<const std::string&, const std::string&, double, double, int, int>(
        (
          arg("name"),
          arg("wallName"),
          arg("normalK"),
          arg("shearK"),
          arg("particleTag"),
          arg("tagMask")
        )
        )
      );



      boost::python::class_<NRotElasticWallPrmsPy,boost::python::bases<InteractionPrmsPy> >(
        "NRotElasticWallPrms",
        "Parameters for linear elastic contact between particles and"
        " a wall.",
        boost::python::init<
          const std::string&,
          const std::string&,
          double
        >(
          (
            arg("name"),
            arg("wallName"),
            arg("normalK")
          ),
	  "Parameters defining elastic contacts between particles and a planar wall\n"
          "@type name: string\n"
          "@kwarg name: Name assigned to the created interaction group.\n"
          "@type wallName: string\n"
          "@kwarg wallName: The name of an existing wall.\n"
          "@type normalK: float\n"
          "@kwarg normalK: spring constant for the linear elastic force"
          " calculation.\n"
        )
      )
      ;

      
    }
  }
}

