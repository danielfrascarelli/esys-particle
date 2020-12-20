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


#ifndef ESYS_LSMCHECKPOINTINFO_H
#define ESYS_LSMCHECKPOINTINFO_H

#include <vector>
#include <iostream>

namespace esys
{
  namespace lsm
  {
    class GeometryInfo;
    typedef std::vector<std::string> StringVector;
    /**
     * Summary of LSM check-point information.
     */
    class CheckPointInfo
    {
    public:
      CheckPointInfo();

      ~CheckPointInfo();
      
      bool operator==(const CheckPointInfo &cpInfo) const;

      const GeometryInfo &getGeometryInfo() const;
      void setGeometryInfo(const GeometryInfo &geoInfo);
      
      const StringVector &getLatticeDataFiles() const;
      void setLatticeDataFiles(const StringVector &fileNames);
      
      int getNumTimeSteps() const;
      void setNumTimeSteps(int numTimeSteps);

      int getTimeStep() const;
      void setTimeStep(int timeStep);

      double getTimeStepSize() const;
      void setTimeStepSize(double timeStepSize);

      void read(std::istream &iStream);
      void write(std::ostream &oStream) const;

    protected:
      CheckPointInfo(const CheckPointInfo &cpInfo);
      CheckPointInfo &operator=(const CheckPointInfo &cpInfo);

    private:
      class Impl;
      Impl  *m_pImpl;
    };
  }
}

#endif
