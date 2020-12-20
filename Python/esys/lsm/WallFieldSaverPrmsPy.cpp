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
#include "Python/esys/lsm/WallFieldSaverPrmsPy.h"
#include "Python/BoostPythonUtil/ListConverter.h"

namespace esys
{
  namespace lsm
  {
    WallFieldSaverPrmsPy::WallFieldSaverPrmsPy(
      const std::string &wallName,
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
      m_wallNameVector()
    {
      StringVector nameVec;
      nameVec.push_back(wallName);
      setWallNames(nameVec);
    }

    WallFieldSaverPrmsPy::WallFieldSaverPrmsPy(
      const boost::python::list &wallNameList,
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
      m_wallNameVector()
    {
      setWallNames(wallNameList);
    }

    WallFieldSaverPrmsPy::WallFieldSaverPrmsPy(
      const boost::python::tuple &wallNameTuple,
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
      m_wallNameVector()
    {
      setWallNames(wallNameTuple);
    }

    const WallFieldSaverPrmsPy::StringVector &
    WallFieldSaverPrmsPy::getWallNameVector() const
    {
      return m_wallNameVector;
    }

    boost::python::list WallFieldSaverPrmsPy::getWallNameList() const
    {
      return bpu::vectorToList(getWallNameVector());
    }

    void WallFieldSaverPrmsPy::setWallNames(
      const StringVector &wallNameVector
    )
    {
      m_wallNameVector = wallNameVector;
    }

    void WallFieldSaverPrmsPy::setWallNames(
      const boost::python::list &wallNameList
    )
    {
      setWallNames(bpu::listToVector<std::string>(wallNameList));
    }

    void WallFieldSaverPrmsPy::setWallNames(
      const boost::python::tuple &wallNameTuple
    )
    {
      setWallNames(bpu::tupleToVector<std::string>(wallNameTuple));
    }

    WallVectorFieldSaverPrmsPy::WallVectorFieldSaverPrmsPy(
      const std::string &wallName,
      const std::string &fieldName,
      const std::string &fileName,
      const std::string &fileFormat,
      int beginTimeStep,
      int endTimeStep,
      int timeStepIncr
    ) :
      WallFieldSaverPrmsPy(
        wallName,
        fieldName,
        fileName,
        fileFormat,
        beginTimeStep,
        endTimeStep,
        timeStepIncr
      )
    {
    }

    WallVectorFieldSaverPrmsPy::WallVectorFieldSaverPrmsPy(
      const boost::python::list &wallNameList,
      const std::string &fieldName,
      const std::string &fileName,
      const std::string &fileFormat,
      int beginTimeStep,
      int endTimeStep,
      int timeStepIncr
    ) :
      WallFieldSaverPrmsPy(
        wallNameList,
        fieldName,
        fileName,
        fileFormat,
        beginTimeStep,
        endTimeStep,
        timeStepIncr
      )
    {
    }

    WallVectorFieldSaverPrmsPy::WallVectorFieldSaverPrmsPy(
      const boost::python::tuple &wallNameTuple,
      const std::string &fieldName,
      const std::string &fileName,
      const std::string &fileFormat,
      int beginTimeStep,
      int endTimeStep,
      int timeStepIncr
    ) :
      WallFieldSaverPrmsPy(
        wallNameTuple,
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
    void exportWallFieldSaverPrms()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<
        WallFieldSaverPrmsPy,
        boost::python::bases<FieldSaverPrmsPy>
      >(
        "WallFieldSaverPrms",
        "Base class describing parameters for saving wall-data to"
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
            arg("wallName"),
            arg("fieldName"),
            arg("fileName"),
            arg("fileFormat"),
            arg("beginTimeStep"),
            arg("endTimeStep"),
            arg("timeStepIncr")
          )
        )
      )
        .def(
          boost::python::init<
            const boost::python::tuple &,
            const std::string &,
            const std::string &,
            const std::string &,
            int,
            int,
            int
          >(
            (
              arg("wallName"),
              arg("fieldName"),
              arg("fileName"),
              arg("fileFormat"),
              arg("beginTimeStep"),
              arg("endTimeStep"),
              arg("timeStepIncr")
            )
          )
        )
        .def(
          boost::python::init<
            const boost::python::list &,
            const std::string &,
            const std::string &,
            const std::string &,
            int,
            int,
            int
          >(
            (
              arg("wallName"),
              arg("fieldName"),
              arg("fileName"),
              arg("fileFormat"),
              arg("beginTimeStep"),
              arg("endTimeStep"),
              arg("timeStepIncr")
            ),
	  "Base class for wall field savers\n"
          "@type wallName: str or tuple or list of str\n"
          "@kwarg wallName: Name of the wall(s) for"
          " which data are saved.\n"
          "@type fieldName: string\n"
          "@kwarg fieldName: Name of the data field, e.g. 'Position'\n"
          "@type fileName: string\n"
          "@kwarg fileName: Name of the file where data are saved\n"
          "@type fileFormat: string\n"
          "@kwarg fileFormat: Format of the data, e.g. 'RAW_SERIES'\n"
          "@type beginTimeStep: int\n"
          "@kwarg beginTimeStep: start saving data at this time step\n"
          "@type endTimeStep: int\n"
          "@kwarg endTimeStep:   finish saving data at this time step\n"
          "@type timeStepIncr: int\n"
          "@kwarg timeStepIncr:  save data every timeStepIncr time steps"
          )
        )
        .def(
          "getWallNames",
          &WallFieldSaverPrmsPy::getWallNameList,
          "Returns the list of wall names associated with these saver parameters.\n"
          "@rtype: list of str\n"
          "@return: Wall name list."
        )
      ;

      boost::python::class_<
        WallVectorFieldSaverPrmsPy,
        boost::python::bases<WallFieldSaverPrmsPy>
      >(
        "WallVectorFieldSaverPrms",
        "Parameters for saving vector wall-data to file.",
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
            arg("wallName"),
            arg("fieldName"),
            arg("fileName"),
            arg("fileFormat"),
            arg("beginTimeStep"),
            arg("endTimeStep"),
            arg("timeStepIncr")
          )
        )
      )
        .def(
          boost::python::init<
            const boost::python::tuple &,
            const std::string &,
            const std::string &,
            const std::string &,
            int,
            int,
            int
          >(
            (
              arg("wallName"),
              arg("fieldName"),
              arg("fileName"),
              arg("fileFormat"),
              arg("beginTimeStep"),
              arg("endTimeStep"),
              arg("timeStepIncr")
            )
          )
        )
        .def(
          boost::python::init<
            const boost::python::list &,
            const std::string &,
            const std::string &,
            const std::string &,
            int,
            int,
            int
          >(
            (
              arg("wallName"),
              arg("fieldName"),
              arg("fileName"),
              arg("fileFormat"),
              arg("beginTimeStep"),
              arg("endTimeStep"),
              arg("timeStepIncr")
            ),
          "Defines parameters for a vector wall field saver\n"
          "@type wallName: str or tuple or list of str\n"
          "@kwarg wallName: Name of the wall(s) for"
          " which data are saved.\n"
          "@type fieldName: string\n"
          "@kwarg fieldName: Name of the data field, e.g. 'Position'\n"
          "@type fileName: string\n"
          "@kwarg fileName: Name of the file where data are saved\n"
          "@type fileFormat: string\n"
          "@kwarg fileFormat: Format of the data, e.g. 'RAW_SERIES'\n"
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
