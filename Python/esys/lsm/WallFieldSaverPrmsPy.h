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


#ifndef ESYS_LSMWALLFIELDSAVERPRMSPY_H
#define ESYS_LSMWALLFIELDSAVERPRMSPY_H

#include "Python/esys/lsm/FieldSaverPrmsPy.h"

#include <boost/python.hpp>
#include <string>
#include <vector>

namespace esys
{
  namespace lsm
  {
    class WallFieldSaverPrmsPy : public FieldSaverPrmsPy
    {
    public:
      typedef std::vector<std::string> StringVector;

      WallFieldSaverPrmsPy(
        const std::string &wallName,
        const std::string &fieldName,
        const std::string &fileName,
        const std::string &fileFormat,
        int beginTimeStep,
        int endTimeStep,
        int timeStepIncr
      );

      WallFieldSaverPrmsPy(
        const boost::python::list &wallNameList,
        const std::string &fieldName,
        const std::string &fileName,
        const std::string &fileFormat,
        int beginTimeStep,
        int endTimeStep,
        int timeStepIncr
      );

      WallFieldSaverPrmsPy(
        const boost::python::tuple &wallNameList,
        const std::string &fieldName,
        const std::string &fileName,
        const std::string &fileFormat,
        int beginTimeStep,
        int endTimeStep,
        int timeStepIncr
      );

      boost::python::list getWallNameList() const;

      const std::vector<std::string> &getWallNameVector() const;

    protected:
      void setWallNames(const StringVector &wallNameVec);

      void setWallNames(const boost::python::list &wallNameList);

      void setWallNames(const boost::python::tuple &wallNameTuple);

    private:
      StringVector m_wallNameVector;
    };

    class WallVectorFieldSaverPrmsPy : public WallFieldSaverPrmsPy
    {
    public:
      WallVectorFieldSaverPrmsPy(
        const std::string &wallName,
        const std::string &fieldName,
        const std::string &fileName,
        const std::string &fileFormat,
        int beginTimeStep,
        int endTimeStep,
        int timeStepIncr
      );

      WallVectorFieldSaverPrmsPy(
        const boost::python::list &wallNameList,
        const std::string &fieldName,
        const std::string &fileName,
        const std::string &fileFormat,
        int beginTimeStep,
        int endTimeStep,
        int timeStepIncr
      );

      WallVectorFieldSaverPrmsPy(
        const boost::python::tuple &wallNameTuple,
        const std::string &fieldName,
        const std::string &fileName,
        const std::string &fileFormat,
        int beginTimeStep,
        int endTimeStep,
        int timeStepIncr
      );
    };

    void exportWallFieldSaverPrms();
  } // namespace lsm
} // namespace esys

#endif
