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

#ifndef ESYS_LSMCHECKPOINTPARAMSPY_H
#define ESYS_LSMCHECKPOINTPARAMSPY_H

//--- STL includes ---
#include <string>

//--- Boost includes ---
#include <boost/python.hpp>

namespace esys
{
  namespace lsm
  {
   
    /*!
      \class CheckPointPrmsPy
      \brief
      
      $Revision$
      $Date$
    */
    class CheckPointPrmsPy
    {
    private:
      std::string m_fileNamePrefix;
      int         m_beginTimeStep;
      int         m_endTimeStep;
      int         m_timeStepIncr;
      
    protected:
      std::string getFileName(int, int rank=0) const;
      
    public:
      CheckPointPrmsPy(const std::string&,int, int, int);
      std::string getFileNamePrefix() const {return m_fileNamePrefix;};
      int getBeginTimeStep() const  {return m_beginTimeStep;};
      int getEndTimeStep() const {return m_endTimeStep;};
      int getTimeStepIncr() const {return m_timeStepIncr;};
      boost::python::list getFileNameList() const; 
    }; // class
    
    /*!
      \class RestartCheckPointPrmsPy
      \brief Parameter class for restart checkpointers, differs from CheckPointPrmsPy by having an additional "binary" flag
    */
    class RestartCheckPointPrmsPy : public CheckPointPrmsPy
    {
    private:
      int m_Precision;

    public:
      RestartCheckPointPrmsPy(const std::string&,int, int, int, int);
      RestartCheckPointPrmsPy(const std::string&,int, int, int);

      int getPrecision() const {return  m_Precision;};
     
    };


    void exportCheckPointPrms();
  } // namespace lsm
} // namespace esys

#endif //ESYS_LSMCHECKPOINTPARAMSPY_H
