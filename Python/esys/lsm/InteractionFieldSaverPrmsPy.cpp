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
#include "Python/esys/lsm/InteractionFieldSaverPrmsPy.h"

namespace esys
{
  namespace lsm
  {
    InteractionFieldSaverPrmsPy::InteractionFieldSaverPrmsPy(
      const std::string &interactionName,
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
      ),
      m_interactionName(interactionName)
    {
    }

    InteractionScalarFieldSaverPrmsPy::InteractionScalarFieldSaverPrmsPy(
      const std::string &interactionName,
      const std::string &fieldName,
      const std::string &fileName,
      const std::string &fileFormat,
      int beginTimeStep,
      int endTimeStep,
      int timeStepIncr
    ) :
      InteractionFieldSaverPrmsPy(
        interactionName,
        fieldName,
        fileName,
        fileFormat,
        beginTimeStep,
        endTimeStep,
        timeStepIncr
      )
    {
    }

    CheckedInteractionScalarFieldSaverPrmsPy::CheckedInteractionScalarFieldSaverPrmsPy(
      const std::string &interactionName,
      const std::string &fieldName,
      const std::string &fileName,
      const std::string &fileFormat,
      int beginTimeStep,
      int endTimeStep,
      int timeStepIncr
    ) :
      InteractionFieldSaverPrmsPy(
        interactionName,
        fieldName,
        fileName,
        fileFormat,
        beginTimeStep,
        endTimeStep,
        timeStepIncr
      )
    {
    }

    TaggedInteractionScalarFieldSaverPrmsPy::TaggedInteractionScalarFieldSaverPrmsPy(
      const std::string &interactionName,
      const std::string &fieldName,
      const std::string &fileName,
      const std::string &fileFormat,
      int beginTimeStep,
      int endTimeStep,
      int timeStepIncr,
      int tag,
      int mask
    ) :
      InteractionScalarFieldSaverPrmsPy(
        interactionName,
        fieldName,
        fileName,
        fileFormat,
        beginTimeStep,
        endTimeStep,
        timeStepIncr
      )
    {
      m_tag=tag;
      m_mask=mask;
    }

    

    InteractionVectorFieldSaverPrmsPy::InteractionVectorFieldSaverPrmsPy(
      const std::string &interactionName,
      const std::string &fieldName,
      const std::string &fileName,
      const std::string &fileFormat,
      int beginTimeStep,
      int endTimeStep,
      int timeStepIncr
    ) :
      InteractionFieldSaverPrmsPy(
        interactionName,
        fieldName,
        fileName,
        fileFormat,
        beginTimeStep,
        endTimeStep,
        timeStepIncr
      )
    {
    }

    CheckedInteractionVectorFieldSaverPrmsPy::CheckedInteractionVectorFieldSaverPrmsPy(
      const std::string &interactionName,
      const std::string &fieldName,
      const std::string &fileName,
      const std::string &fileFormat,
      int beginTimeStep,
      int endTimeStep,
      int timeStepIncr
    ) :
      InteractionFieldSaverPrmsPy(
        interactionName,
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
    void exportInteractionFieldSaverPrms()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<
        InteractionFieldSaverPrmsPy,
        boost::python::bases<FieldSaverPrmsPy>
      >(
        "InteractionFieldSaverPrms",
        "Base class describing parameters for saving interaction-data to"
        " file.",
        boost::python::init<
          const std::string &,
          const std::string &,
          const std::string &,
          const std::string &,
          int,
          int,
          int
        >(
          (
            arg("interactionName"),
            arg("fieldName"),
            arg("fileName"),
            arg("fileFormat"),
            arg("beginTimeStep"),
            arg("endTimeStep"),
            arg("timeStepIncr")
          ),
	  "Base class for interaction field savers\n"
          "@type interactionName: string\n"
          "@kwarg interactionName: Name of the interaction group for"
          " which data are saved\n"
          "@type fieldName: string\n"
          "@kwarg fieldName: Name of the data field\n"
          "@type fileName: string\n"
          "@kwarg fileName: Name of the file where data are saved\n"
          "@type fileFormat: string\n"
          "@kwarg fileFormat: Format of the data\n"
          "@type beginTimeStep: int\n"
          "@kwarg beginTimeStep: start saving data at this time step\n"
          "@type endTimeStep: int\n"
          "@kwarg endTimeStep:   finish saving data at this time step\n"
          "@type timeStepIncr: int\n"
          "@kwarg timeStepIncr:  save data every timeStepIncr time steps\n"
        )
      )
        .def(
          "getInteractionName",
          &InteractionFieldSaverPrmsPy::getInteractionName,
          boost::python::return_value_policy<
            boost::python::copy_const_reference
          >(),
          "Returns the name of the interaction associated with these saver parameters.\n"
          "@rtype: string\n"
          "@return: Interaction name."
        )
      ;

      boost::python::class_<
        InteractionScalarFieldSaverPrmsPy,
        boost::python::bases<InteractionFieldSaverPrmsPy>
      >(
        "InteractionScalarFieldSaverPrms",
        "Parameters for saving scalar interaction-data to file.",
        boost::python::init<
          const std::string &,
          const std::string &,
          const std::string &,
          const std::string &,
          int,
          int,
          int
        >(
          (
            arg("interactionName"),
            arg("fieldName"),
            arg("fileName"),
            arg("fileFormat"),
            arg("beginTimeStep"),
            arg("endTimeStep"),
            arg("timeStepIncr")
          ),
          "Defines interaction scalar field saver parameters\n"
          "@type interactionName: string\n"
          "@kwarg interactionName: Name of the interaction group for"
          " which data are saved\n"
          "@type fieldName: string\n"
          "@kwarg fieldName: Name of the data field, e.g. 'e_pot_normal'\n"
          "@type fileName: string\n"
          "@kwarg fileName: Name of the file where data are saved\n"
          "@type fileFormat: string\n"
          "@kwarg fileFormat: Format of the data, e.g. 'SUM' or 'RAW2'\n"
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
        CheckedInteractionScalarFieldSaverPrmsPy,
        boost::python::bases<InteractionFieldSaverPrmsPy>
      >(
        "CheckedInteractionScalarFieldSaverPrms",
        "Parameters for saving checked scalar interaction-data to file.",
        boost::python::init<
          const std::string &,
          const std::string &,
          const std::string &,
          const std::string &,
          int,
          int,
          int
        >(
          (
            arg("interactionName"),
            arg("fieldName"),
            arg("fileName"),
            arg("fileFormat"),
            arg("beginTimeStep"),
            arg("endTimeStep"),
            arg("timeStepIncr")
          ),
          "Defines checked interaction scalar field saver parameters\n"
          "@type interactionName: string\n"
          "@kwarg interactionName: Name of the interaction group for"
          " which data are saved\n"
          "@type fieldName: string\n"
          "@kwarg fieldName: Name of the data field, e.g. 'e_pot_normal'\n"
          "@type fileName: string\n"
          "@kwarg fileName: Name of the file where data are saved\n"
          "@type fileFormat: string\n"
          "@kwarg fileFormat: Format of the data, e.g. 'SUM' or 'RAW2'\n"
          "@type beginTimeStep: int\n"
          "@kwarg beginTimeStep: start saving data at this time step\n"
          "@type endTimeStep: int\n"
          "@kwarg endTimeStep:   finish saving data at this time step\n"
          "@type timeStepIncr: int\n"
          "@kwarg timeStepIncr:  save data every timeStepIncr time steps\n"
        )
      );

      boost::python::class_<
        TaggedInteractionScalarFieldSaverPrmsPy,
        boost::python::bases<InteractionFieldSaverPrmsPy>
      >(
        "TaggedInteractionScalarFieldSaverPrms",
        "Parameters for saving scalar data on tagged interactions to file.",
        boost::python::init<
          const std::string &,
          const std::string &,
          const std::string &,
          const std::string &,
          int,
          int,
	  int,
	  int,
	  int
        >(
          (
            arg("interactionName"),
            arg("fieldName"),
            arg("fileName"),
            arg("fileFormat"),
            arg("beginTimeStep"),
            arg("endTimeStep"),
            arg("timeStepIncr"),
	    arg("tag"),
	    arg("mask")
          ),
          "Defines tagged interaction scalar field saver parameters\n"
	  "@type interactionName: string\n"
          "@kwarg interactionName: Name of the interaction group for"
          " which data are saved\n"
          "@type fieldName: string\n"
          "@kwarg fieldName: Name of the data field, e.g. 'e_pot_normal'\n"
          "@type fileName: string\n"
          "@kwarg fileName: Name of the file where data are saved\n"
          "@type fileFormat: string\n"
          "@kwarg fileFormat: Format of the data, e.g. 'SUM' or 'RAW2'\n"
          "@type beginTimeStep: int\n"
          "@kwarg beginTimeStep: start saving data at this time step\n"
          "@type endTimeStep: int\n"
          "@kwarg endTimeStep:   finish saving data at this time step\n"
          "@type timeStepIncr: int\n"
          "@kwarg timeStepIncr:  save data every timeStepIncr time steps\n"
	  "@type tag: int\n"
	  "@kwarg tag: the tag of the particles\n"
	  "@type mask: int\n"
	  "@kwarg mask: the tag mask\n"
        )
      )
      ;

      boost::python::class_<
        InteractionVectorFieldSaverPrmsPy,
        boost::python::bases<InteractionFieldSaverPrmsPy>
      >(
        "InteractionVectorFieldSaverPrms",
        "Parameters for saving vector interaction-data to file.",
        boost::python::init<
          const std::string &,
          const std::string &,
          const std::string &,
          const std::string &,
          int,
          int,
          int
        >(
          (
            arg("interactionName"),
            arg("fieldName"),
            arg("fileName"),
            arg("fileFormat"),
            arg("beginTimeStep"),
            arg("endTimeStep"),
            arg("timeStepIncr")
          ),
          "Defines interaction vector field saver parameters\n"
	  "@type interactionName: string\n"
          "@kwarg interactionName: Name of the interaction group for"
          " which data are saved\n"
          "@type fieldName: string\n"
          "@kwarg fieldName: Name of the data field, e.g. 'normal_force'\n"
          "@type fileName: string\n"
          "@kwarg fileName: Name of the file where data are saved\n"
          "@type fileFormat: string\n"
          "@kwarg fileFormat: Format of the data, e.g. 'SUM' or 'RAW2'\n"
          "@type beginTimeStep: int\n"
          "@kwarg beginTimeStep: start saving data at this time step\n"
          "@type endTimeStep: int\n"
          "@kwarg endTimeStep:   finish saving data at this time step\n"
          "@type timeStepIncr: int\n"
          "@kwarg timeStepIncr:  save data every timeStepIncr time steps\n"
        )
      );


      boost::python::class_<
        CheckedInteractionVectorFieldSaverPrmsPy,
        boost::python::bases<InteractionFieldSaverPrmsPy>
      >(
        "CheckedInteractionVectorFieldSaverPrms",
        "Parameters for saving checked vector interaction-data to file.",
        boost::python::init<
          const std::string &,
          const std::string &,
          const std::string &,
          const std::string &,
          int,
          int,
          int
        >(
          (
            arg("interactionName"),
            arg("fieldName"),
            arg("fileName"),
            arg("fileFormat"),
            arg("beginTimeStep"),
            arg("endTimeStep"),
            arg("timeStepIncr")
          ),
          "Defines checked interaction vector field saver parameters\n"
	  "@type interactionName: string\n"
          "@kwarg interactionName: Name of the interaction group for"
          " which data are saved\n"
          "@type fieldName: string\n"
          "@kwarg fieldName: Name of the data field, e.g. 'normal_force'\n"
          "@type fileName: string\n"
          "@kwarg fileName: Name of the file where data are saved\n"
          "@type fileFormat: string\n"
          "@kwarg fileFormat: Format of the data, e.g. 'SUM' or 'RAW2'\n"
          "@type beginTimeStep: int\n"
          "@kwarg beginTimeStep: start saving data at this time step\n"
          "@type endTimeStep: int\n"
          "@kwarg endTimeStep:   finish saving data at this time step\n"
          "@type timeStepIncr: int\n"
          "@kwarg timeStepIncr:  save data every timeStepIncr time steps\n"
        )
      );

    }
  }
}
