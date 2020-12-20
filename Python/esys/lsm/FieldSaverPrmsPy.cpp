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
#include "Python/esys/lsm/FieldSaverPrmsPy.h"

namespace esys
{
  namespace lsm
  {
    FieldSaverPrmsPy::FieldSaverPrmsPy(
      const std::string &fieldName,
      const std::string &fileName,
      const std::string &fileFormat,
      int beginTimeStep,
      int endTimeStep,
      int timeStepIncr
    ) :
      m_fieldName(fieldName),
      m_fileName(fileName),
      m_fileFormat(fileFormat),
      m_beginTimeStep(beginTimeStep),
      m_endTimeStep(endTimeStep),
      m_timeStepIncr(timeStepIncr)
    {
    }

    using boost::python::arg;
    void exportFieldSaverPrms()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<FieldSaverPrmsPy>(
        "FieldSaverPrms",
        "Base class for parameters defining field savers to store model data.",
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
          "Base class for field saver parameters.\n"
          "@type fieldName: string\n"
          "@kwarg fieldName: Name of a data-field.\n"
          "@type fileName: string\n"
          "@kwarg fileName: Name of the file where the data is saved.\n"
          "@type fileFormat: string\n"
          "@kwarg fileFormat: Format of the data which is saved to file.\n"
          "@type beginTimeStep: int\n"
          "@kwarg beginTimeStep: Time step when first data is saved.\n"
          "@type endTimeStep: int\n"
          "@kwarg endTimeStep: Time step when last data is saved.\n"
          "@type timeStepIncr: int\n"
          "@kwarg timeStepIncr: Number of time steps between data saves.\n"
        )
      )
        .def(
          "getFieldName",
          &FieldSaverPrmsPy::getFieldName,
          boost::python::return_value_policy<
            boost::python::copy_const_reference
          >(),
          "Returns name of save field.\n"
          "@rtype: string\n"
          "@return: Field name\n"
        )
        .def(
          "getFileName",
          &FieldSaverPrmsPy::getFileName,
          boost::python::return_value_policy<
            boost::python::copy_const_reference
          >(),
          "Returns name of file where data is to be saved.\n"
          "@rtype: string\n"
          "@return: File name.\n"
        )
        .def("getFileFormat",
          &FieldSaverPrmsPy::getFileFormat,
          boost::python::return_value_policy<
            boost::python::copy_const_reference
          >(),
          "Returns format of saved data.\n"
          "@rtype: string\n"
          "@return: Name of saved data format."
        )
        .def(
          "getBeginTimeStep",
          &FieldSaverPrmsPy::getBeginTimeStep,
          "Returns time step when data is first saved.\n"
          "@rtype: int\n"
          "@return: Time step."
        )
        .def(
          "getEndTimeStep",
          &FieldSaverPrmsPy::getEndTimeStep,
          "Returns time step when data is last saved.\n"
          "@rtype: int\n"
          "@return: Time step."
        )
        .def(
          "getTimeStepIncr",
          &FieldSaverPrmsPy::getTimeStepIncr,
          "Returns number of time steps between data saves.\n"
          "@rtype: int\n"
          "@return: number of time steps."
        )
        ;
    }
  }
}
