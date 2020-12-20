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
#include "Python/esys/lsm/geometry/GougeConfigPrmsPy.h"
#include "Python/BoostPythonUtil/ListConverter.h"

namespace esys
{
  namespace lsm
  {
    ParticleRndPackPrmsPy::ParticleRndPackPrmsPy(double size, double minRadius, double maxRadius)
      : Inherited(size, minRadius, maxRadius)
    {
    }

    GrainRndPackPrmsPy::GrainRndPackPrmsPy(
      double             size,
      ParticleGrainGenPy &grainGen,
      int                connTag
    ) : Inherited(size, grainGen, connTag)
    {
    }

    GougeConfigPrmsPy::GougeConfigPrmsPy(
      const BoundingBoxPy       &bBox,
      double                    padRadius,
      const ParticleRndPackPrmsPy        &roughnessPrms,
      const GrainRndPackPrmsPy   &gougePrms,
      int                       maxInsertionFailures,
      const boost::python::list &periodicDimensionList,
      double                    tolerance,
      double                    connectionTolerance,
      int                       blockConnTag
    ) :
        Inherited(
          bBox,
          padRadius,
          XZ,
          roughnessPrms,
          gougePrms,
          bpu::listToVector<bool>(periodicDimensionList),
          maxInsertionFailures,
          tolerance,
          connectionTolerance,
          blockConnTag
        )
    {
    }

    using boost::python::arg;
    void exportGougeConfigPrms()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<ParticleRndPackPrmsPy>(
        "ParticleRndPackPrms",
        boost::python::init<
          double,
          double,
          double
        >(
          ( arg("size"), arg("minRadius"), arg("maxRadius") )
        )
      );

      boost::python::class_<GrainRndPackPrmsPy>(
        "GrainRndPackPrms",
        boost::python::init<
          double,
          GrainRndPackPrmsPy::ParticleGrainGenPy &,
          int
        >(
          ( arg("size"), arg("grainGenerator"), arg("connTag")=0 )
        )
      );

      boost::python::class_<GougeConfigPrmsPy>(
        "GougeConfigPrms",
        boost::python::init<
          const BoundingBoxPy &,
          double,
          const ParticleRndPackPrmsPy &,
          const GrainRndPackPrmsPy &,
          int,
          const boost::python::list &,
          double,
          double,
          int
        >(
          (
            arg("bBox"),
            arg("padRadius"),
            arg("rndBlockPrms"),
            arg("rndGougePrms"),
            arg("maxInsertFails"),
            arg("circDimList")=bpu::vectorToList(BoolVector(3,false)),
            arg("overlapTol")=0.001,
            arg("connectionTol")=0.01,
            arg("blockConnTag")=0
          )
        )
      )
      .def("getMinParticleRadius", &GougeConfigPrmsPy::getMinRadius)
      .def("getMaxParticleRadius", &GougeConfigPrmsPy::getMaxRadius)
      .def("is2d",                 &GougeConfigPrmsPy::is2d)
      ;
    }

  }
}
