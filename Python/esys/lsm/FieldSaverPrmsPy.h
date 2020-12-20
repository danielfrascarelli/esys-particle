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

#ifndef ESYS_LSMFIELDSAVERPRMSPY_H
#define ESYS_LSMFIELDSAVERPRMSPY_H

//--- STL includes ---
#include <string>

//--- Boost includes ---
#include <boost/python.hpp>

namespace esys
{
  namespace lsm
  {
    class FieldSaverPrmsPy
    {
    public:
      FieldSaverPrmsPy(
        const std::string &fieldName,
        const std::string &fileName,
        const std::string &fileFormat,
        int beginTimeStep,
        int endTimeStep,
        int timeStepIncr
      );

      const std::string &getFieldName() const
      {
        return m_fieldName;
      }

      const std::string &getFileName() const
      {
        return m_fileName;
      }

      const std::string &getFileFormat() const
      {
        return m_fileFormat;
      }

      int getBeginTimeStep() const
      {
        return m_beginTimeStep;
      }
      
      int getEndTimeStep() const
      {
        return m_endTimeStep;
      }
      
      int getTimeStepIncr() const
      {
        return m_timeStepIncr;
      }

    private:
      std::string m_fieldName;
      std::string m_fileName;
      std::string m_fileFormat;
      int         m_beginTimeStep;
      int         m_endTimeStep;
      int         m_timeStepIncr;
    }; // class
    
    void exportFieldSaverPrms();
  } // namespace lsm
} // namespace esys

#endif //ESYS_LSMCHECKPOINTPARAMSPY_H
