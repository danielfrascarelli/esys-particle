/////////////////////////////////////////////////////////////
//                                                         //
// Copyright (c) 2003-2014 by The University of Queensland //
// Centre for Geoscience Computing                         //
// http://earth.uq.edu.au/centre-geoscience-computing      //
//                                                         //
// Primary Business: Brisbane, Queensland, Australia       //
// Licensed under the Open Software License version 3.0    //
// http://www.apache.org/licenses/LICENSE-2.0          //
//                                                         //
/////////////////////////////////////////////////////////////


#include <boost/version.hpp>
#include <boost/python.hpp>
#include "Python/esys/lsm/FluidInteractionFieldSaverPrmsPy.h"

namespace esys
{
  namespace lsm
  {
    FluidInteractionScalarFieldSaverPrmsPy::FluidInteractionScalarFieldSaverPrmsPy(
      const std::string &fieldName,
      const std::string &fileName,
      const std::string &fileFormat,
      int beginTimeStep,
      int endTimeStep,
      int timeStepIncr
    ) :
      FieldSaverPrmsPy(
        fieldName,
        fileName,
        fileFormat,
        beginTimeStep,
        endTimeStep,
        timeStepIncr
      )
    {
    }

    FluidInteractionVectorFieldSaverPrmsPy::FluidInteractionVectorFieldSaverPrmsPy(
      const std::string &fieldName,
      const std::string &fileName,
      const std::string &fileFormat,
      int beginTimeStep,
      int endTimeStep,
      int timeStepIncr
    ) :
      FieldSaverPrmsPy(
        fieldName,
        fileName,
        fileFormat,
        beginTimeStep,
        endTimeStep,
        timeStepIncr
      )
    {
    }


    using boost::python::arg;
    void exportFluidInteractionFieldSaverPrms()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<
        FluidInteractionScalarFieldSaverPrmsPy,
        boost::python::bases<FieldSaverPrmsPy>
      >(
        "FluidInteractionScalarFieldSaverPrms",
        "Parameters for saving scalar fluid-particle-interaction-data to file.",
        boost::python::init<
          const std::string &,
          const std::string &,
          const std::string &,
          int,
          int,
          int
        >(
          (
            arg("fieldName"),
            arg("fileName"),
            arg("fileFormat"),
            arg("beginTimeStep"),
            arg("endTimeStep"),
            arg("timeStepIncr")
          ),
	  "Defines fluid-particle-interaction scalar field saver parameters\n"
          "@type fieldName: string\n"
          "@kwarg fieldName: Name of the data field, e.g. 'Force_abs'\n"
          "@type fileName: string\n"
          "@kwarg fileName: Name of the file where data is saved\n"
          "@type fileFormat: string\n"
          "@kwarg fileFormat: Format of the data, e.g. 'SUM' or 'MAX'\n"
          "@type beginTimeStep: int\n"
          "@kwarg beginTimeStep: start saving data at this time step\n"
          "@type endTimeStep: int\n"
          "@kwarg endTimeStep:   finish saving data at this time step\n"
          "@type timeStepIncr: int\n"
          "@kwarg timeStepIncr:  save data every timeStepIncr time steps\n"
        )
      )
      ;

      boost::python::class_<
        FluidInteractionVectorFieldSaverPrmsPy,
        boost::python::bases<FieldSaverPrmsPy>
      >(
        "FluidInteractionVectorFieldSaverPrms",
        "Parameters for saving vector fluid-particle-interaction-data to file.",
        boost::python::init<
          const std::string &,
          const std::string &,
          const std::string &,
          int,
          int,
          int
        >(
          (
            arg("fieldName"),
            arg("fileName"),
            arg("fileFormat"),
            arg("beginTimeStep"),
            arg("endTimeStep"),
            arg("timeStepIncr")
          ),
	  "Defines fluidcell vector field saver parameters\n"
          "@type fieldName: string\n"
          "@kwarg fieldName: Name of the data field, e.g. 'drag_force'\n"
          "@type fileName: string\n"
          "@kwarg fileName: Name of the file where data is saved\n"
          "@type fileFormat: string\n"
          "@kwarg fileFormat: Format of the data, e.g. 'SUM' or 'MAX'\n"
          "@type beginTimeStep: int\n"
          "@kwarg beginTimeStep: start saving data at this time step\n"
          "@type endTimeStep: int\n"
          "@kwarg endTimeStep:   finish saving data at this time step\n"
          "@type timeStepIncr: int\n"
          "@kwarg timeStepIncr:  save data every timeStepIncr time steps\n"
        )
      )
      ;
    }
  }
}
