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
#include "Foundation/StringUtil.h"
#include "Python/esys/lsm/CheckPointParamsPy.h"

namespace esys
{
  namespace lsm
  {
    CheckPointPrmsPy::CheckPointPrmsPy(
      const std::string &fileNamePrefix,
      int beginTimeStep,
      int endTimeStep,
      int timeStepIncr
    )
      : m_fileNamePrefix(fileNamePrefix),
        m_beginTimeStep(beginTimeStep),
        m_endTimeStep(endTimeStep),
        m_timeStepIncr(timeStepIncr)
    {
    }

    boost::python::list CheckPointPrmsPy::getFileNameList() const
    {
      boost::python::list fileNameList;
      for (int i = getBeginTimeStep(); i < getEndTimeStep(); i += getTimeStepIncr())
        {
          fileNameList.append(getFileName(i));
        }
      fileNameList.append(getFileName(getEndTimeStep()));
      return fileNameList;
    }

    std::string CheckPointPrmsPy::getFileName(int timeStep, int rank) const
    {
      return (getFileNamePrefix()+ "_t=" +
              StringUtil::toString(timeStep) + "_" +
              StringUtil::toString(rank) + ".txt");
    }

    /*!
      Construct RestartCheckPointPrmsPy with full set of arguments
    */
    RestartCheckPointPrmsPy::RestartCheckPointPrmsPy(
      const std::string &fileNamePrefix,
      int beginTimeStep,
      int endTimeStep,
      int timeStepIncr,
      int Precision)
      : CheckPointPrmsPy(fileNamePrefix,beginTimeStep,endTimeStep,timeStepIncr), m_Precision(Precision)
    {}

    /*!
      Construct RestartCheckPointPrmsPy with precision defaulting to 12 digits
    */
    RestartCheckPointPrmsPy::RestartCheckPointPrmsPy(
      const std::string &fileNamePrefix,
      int beginTimeStep,
      int endTimeStep,
      int timeStepIncr)
      : CheckPointPrmsPy(fileNamePrefix,beginTimeStep,endTimeStep,timeStepIncr), m_Precision(12)
    {}

    

    using boost::python::arg;
    void exportCheckPointPrms()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<CheckPointPrmsPy>(
        "CheckPointPrms",
        "Parameters for specifying check-pointing intervals.",
        boost::python::init<const std::string &, int, int, int>(
          (
            arg("fileNamePrefix"),
            arg("beginTimeStep"),
            arg("endTimeStep"),
            arg("timeStepIncr")
          ),
	  "Defines CheckPoint snapshot saver parameters\n"
          "@type fileNamePrefix: string\n"
          "@kwarg fileNamePrefix: prefix of files where checkpoint data"
          " is saved.\n"
          "@type beginTimeStep: int\n"
          "@kwarg beginTimeStep: time step when checkpoint saving begins.\n"
          "@type endTimeStep: int\n"
          "@kwarg endTimeStep: time step when checkpoint saving stops.\n"
          "@type timeStepIncr: int\n"
          "@kwarg timeStepIncr: a checkpoint is made every"
          "C{timeStepIncr} time-steps.\n"
        )
      )
        .def(
          "getFileNamePrefix",
          &CheckPointPrmsPy::getFileNamePrefix,
          "@rtype: string\n"
          "@return: prefix of files where checkpoint data is saved.\n"
        )
        .def(
          "getBeginTimeStep",
          &CheckPointPrmsPy::getBeginTimeStep,
          "@rtype: int"
	  "@return: Timestep when snapshot saving commences"
        )
        .def(
          "getEndTimeStep",
          &CheckPointPrmsPy::getEndTimeStep,
          "@rtype: int"
	  "@return: Timestep when snapshot saving ceases"
        )
        .def(
          "getTimeStepIncr",
          &CheckPointPrmsPy::getTimeStepIncr,
          "@rtype: int"
	  "@return: Number of timesteps between each snapshot"
        )
        .def(
          "getFileNameList",
          &CheckPointPrmsPy::getFileNameList,
          "@rtype: list\n"
          "@return: List of 'master' check-point file names."
        )
        ;

      boost::python::class_<RestartCheckPointPrmsPy,boost::python::bases<CheckPointPrmsPy> >(
        "RestartCheckPointPrms",
        "Parameters for specifying check-pointing intervals.",
        boost::python::init<const std::string &, int, int, int, int>(
          (
            arg("fileNamePrefix"),
            arg("beginTimeStep"),
            arg("endTimeStep"),
            arg("timeStepIncr"),
	    arg("Precision")
          ),
	  "Defines CheckPoint snapshot saver parameters with the output precision defaulting "
          "to 12 digits.  This function is provided for backward compatibility.\n"
          "@type fileNamePrefix: string\n"
          "@kwarg fileNamePrefix: prefix of files where checkpoint data"
          " is saved.\n"
          "@type beginTimeStep: int\n"
          "@kwarg beginTimeStep: time step when checkpoint saving begins.\n"
          "@type endTimeStep: int\n"
          "@kwarg endTimeStep: time step when checkpoint saving stops.\n"
          "@type timeStepIncr: int\n"
          "@kwarg timeStepIncr: a checkpoint is made every "
          "C{timeStepIncr} time-steps.\n"
	  "@type Precision: int\n"
	  "@kwarg Precision: the number of digits in the output data, defaulting to twelve "
          "when not specified.\n"
        )
      )
	.def(
	boost::python::init<const std::string &, int, int, int>(
          (
            arg("fileNamePrefix"),
            arg("beginTimeStep"),
            arg("endTimeStep"),
            arg("timeStepIncr")
	  )
        )
      )
        .def(
          "getPrecision",
          &RestartCheckPointPrmsPy::getPrecision,
          "@rtype: int\n"
          "@return: output precision of checkpoint files.\n"
        )
	;
      
    }
  }
}
