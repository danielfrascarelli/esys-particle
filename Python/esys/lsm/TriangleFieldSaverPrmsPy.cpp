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

#include "Python/esys/lsm/TriangleFieldSaverPrmsPy.h"
#include <boost/version.hpp>
#include <boost/python.hpp>


namespace esys
{
  namespace lsm
  {
    TriangleScalarFieldSaverPrmsPy::TriangleScalarFieldSaverPrmsPy(
      const std::string& meshName,
      const std::string& fieldName,
      const std::string& fileName,
      const std::string& fileFormat,
      int beginTimeStep,
      int endTimeStep,
      int timeStepIncr
    ):
      FieldSaverPrmsPy(
        fieldName,
        fileName,
        fileFormat,
        beginTimeStep,
        endTimeStep,
        timeStepIncr
      ),
      m_MeshName(meshName)
    {}

    TriangleVectorFieldSaverPrmsPy::TriangleVectorFieldSaverPrmsPy(
      const std::string& meshName,
      const std::string& fieldName,
      const std::string& fileName,
      const std::string& fileFormat,
      int beginTimeStep,
      int endTimeStep,
      int timeStepIncr
    ):
      FieldSaverPrmsPy(
        fieldName,
        fileName,
        fileFormat,
        beginTimeStep,
        endTimeStep,
        timeStepIncr
      ),
      m_MeshName(meshName){}

    using boost::python::arg;
    void exportTriangleFieldSaverPrms()
    {
      // Disable autogeneration of C++ signatures (Boost 1.34.0 and higher)
      // for Epydoc which stumbles over indentation in the automatically generated strings.
      boost::python::docstring_options no_autogen(true,false);

      boost::python::class_<
        TriangleScalarFieldSaverPrmsPy,
        boost::python::bases<FieldSaverPrmsPy>
      >(
        "TriangleScalarFieldSaverPrms",
        "Class describing parameters for saving data located on a triangle mesh to"
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
            arg("meshName"),
            arg("fieldName"),
            arg("fileName"),
            arg("fileFormat"),
            arg("beginTimeStep"),
            arg("endTimeStep"),
            arg("timeStepIncr")
          ),
	  "Parameters defining a scalar field saver for triangle meshes\n"
          "@type meshName: string\n"
          "@kwarg meshName: Name of the mesh for which data are saved\n"
          "@type fieldName: string\n"
          "@kwarg fieldName: Name of the data field - 'pressure'\n"
          "@type fileName: string\n"
          "@kwarg fileName: Name of the file where data are saved\n"
          "@type fileFormat: string\n"
          "@kwarg fileFormat: Format of the data - 'RAW'\n"
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
        TriangleVectorFieldSaverPrmsPy,
        boost::python::bases<FieldSaverPrmsPy>
      >(
        "TriangleVectorFieldSaverPrms",
        "Class describing parameters for saving data located on a triangle mesh to"
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
            arg("meshName"),
            arg("fieldName"),
            arg("fileName"),
            arg("fileFormat"),
            arg("beginTimeStep"),
            arg("endTimeStep"),
            arg("timeStepIncr")
          ),
	  "Parameters defining a vector field saver for triangle meshes\n"
          "@type meshName: string\n"
          "@kwarg meshName: Name of the mesh for which data are saved\n"
          "@type fieldName: string\n"
          "@kwarg fieldName: Name of the data field - 'force'\n"
          "@type fileName: string\n"
          "@kwarg fileName: Name of the file where data are saved\n"
          "@type fileFormat: string\n"
          "@kwarg fileFormat: Format of the data - 'RAW'\n"
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
  } // namespace lsm
} // namespace esys
