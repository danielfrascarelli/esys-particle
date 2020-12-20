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
#include "Python/esys/lsm/ParticleFieldSaverPrmsPy.h"

namespace esys
{
  namespace lsm
  {
    ParticleFieldSaverPrmsPy::ParticleFieldSaverPrmsPy(
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

    ParticleScalarFieldSaverPrmsPy::ParticleScalarFieldSaverPrmsPy(
      const std::string &fieldName,
      const std::string &fileName,
      const std::string &fileFormat,
      int beginTimeStep,
      int endTimeStep,
      int timeStepIncr
    ) :
      ParticleFieldSaverPrmsPy(
        fieldName,
        fileName,
        fileFormat,
        beginTimeStep,
        endTimeStep,
        timeStepIncr
      )
    {
    }

    ParticleVectorFieldSaverPrmsPy::ParticleVectorFieldSaverPrmsPy(
      const std::string &fieldName,
      const std::string &fileName,
      const std::string &fileFormat,
      int beginTimeStep,
      int endTimeStep,
      int timeStepIncr
    ) :
      ParticleFieldSaverPrmsPy(
        fieldName,
        fileName,
        fileFormat,
        beginTimeStep,
        endTimeStep,
        timeStepIncr
      )
    {
    }

    TaggedParticleScalarFieldSaverPrmsPy::TaggedParticleScalarFieldSaverPrmsPy(
      const std::string &fieldName,
      const std::string &fileName,
      const std::string &fileFormat,
      int beginTimeStep,
      int endTimeStep,
      int timeStepIncr,
      int tag,
      int mask
    ) : ParticleScalarFieldSaverPrmsPy (fieldName,
					fileName,
					fileFormat,
					beginTimeStep,
					endTimeStep,
					timeStepIncr)
    {
      m_tag=tag;
      m_mask=mask;
    }
					    
    TaggedParticleVectorFieldSaverPrmsPy::TaggedParticleVectorFieldSaverPrmsPy(
      const std::string &fieldName,
      const std::string &fileName,
      const std::string &fileFormat,
      int beginTimeStep,
      int endTimeStep,
      int timeStepIncr,
      int tag,
      int mask
    ) : ParticleVectorFieldSaverPrmsPy (fieldName,
					fileName,
					fileFormat,
					beginTimeStep,
					endTimeStep,
					timeStepIncr)
    {
      m_tag=tag;
      m_mask=mask;
    }
					    

    using boost::python::arg;
    void exportParticleFieldSaverPrms()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<
        ParticleFieldSaverPrmsPy,
        boost::python::bases<FieldSaverPrmsPy>
      >(
        "ParticleFieldSaverPrms",
        "Base class describing parameters for saving particle-data to file.",
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
	  "Base class for particle field savers\n"
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
      ;

      boost::python::class_<
        ParticleScalarFieldSaverPrmsPy,
        boost::python::bases<ParticleFieldSaverPrmsPy>
      >(
        "ParticleScalarFieldSaverPrms",
        "Parameters for saving scalar particle-data to file.",
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
	  "Defines particle scalar field saver parameters\n"
          "@type fieldName: string\n"
          "@kwarg fieldName: Name of the data field, e.g. 'e_kin'\n"
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
        ParticleVectorFieldSaverPrmsPy,
        boost::python::bases<ParticleFieldSaverPrmsPy>
      >(
        "ParticleVectorFieldSaverPrms",
        "Parameters for saving vector particle-data to file.",
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
	  "Defines particle vector field saver parameters\n"
          "@type fieldName: string\n"
          "@kwarg fieldName: Name of the data field, e.g. 'force'\n"
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
      );

      boost::python::class_<
        TaggedParticleScalarFieldSaverPrmsPy,
        boost::python::bases<ParticleScalarFieldSaverPrmsPy>
      >(
        "TaggedParticleScalarFieldSaverPrms",
        "Parameters for saving scalar data of tagged particles to file.",
        boost::python::init<
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
            arg("fieldName"),
            arg("fileName"),
            arg("fileFormat"),
            arg("beginTimeStep"),
            arg("endTimeStep"),
            arg("timeStepIncr"),
            arg("tag"),
            arg("mask")
          ),
	  "Defines tagged particle scalar field saver parameters\n"
          "@type fieldName: string\n"
          "@kwarg fieldName: Name of the data field, e.g. 'e_kin'\n"
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
          "@type tag: int\n"
          "@kwarg tag: the tag of the particles\n"
          "@type mask: int\n"
          "@kwarg mask: the tag mask\n"
        )
      )
      ;

      boost::python::class_<
        TaggedParticleVectorFieldSaverPrmsPy,
        boost::python::bases<ParticleVectorFieldSaverPrmsPy>
      >(
        "TaggedParticleVectorFieldSaverPrms",
        "Parameters for saving vector data of tagged particles to file.",
        boost::python::init<
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
            arg("fieldName"),
            arg("fileName"),
            arg("fileFormat"),
            arg("beginTimeStep"),
            arg("endTimeStep"),
            arg("timeStepIncr"),
            arg("tag"),
            arg("mask")
          ),
	  "Defines tagged particle vector field saver parameters\n"
          "@type fieldName: string\n"
          "@kwarg fieldName: Name of the data field, e.g. 'force'\n"
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
           "@type tag: int\n"
          "@kwarg tag: the tag of the particles\n"
          "@type mask: int\n"
          "@kwarg mask: the tag mask\n"
        )
      );
    }
  }
}
