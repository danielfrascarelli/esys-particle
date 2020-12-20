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


#include "Foundation/StringUtil.h"
#include "Parallel/CheckPointInfo.h"
#include "Geometry/GeometryInfo.h"
#include "Foundation/version.h"
#include "bzrversion.h"

namespace esys
{
  namespace lsm
  {
    class CheckPointInfo::Impl
    {
    public:
      Impl()
        : m_numTimeSteps(0),
          m_timeStepSize(0.0),
          m_timeStep(0),
          m_geoInfo(),
          m_fileNames(),
	  m_bzrVersion(s_bzr_revision)
      {
      }

      ~Impl()
      {
      }

      bool operator==(const Impl &impl) const
      {
        return
          (
            (m_numTimeSteps == impl.m_numTimeSteps)
            &&
            (m_timeStepSize == impl.m_timeStepSize)
            &&
            (m_timeStep     == impl.m_timeStep)
            &&
            (m_geoInfo      == impl.m_geoInfo)
            &&
            (m_fileNames    == impl.m_fileNames)
          );
      }
      
      void read(std::istream &iStream)
      {
	string magic;
	int version;
        iStream	>> magic >> version; 
	iStream 
          >> m_numTimeSteps
          >> m_timeStepSize
          >> m_timeStep;
          m_geoInfo.read(iStream);

          std::string fileNames;
          /*
           * Skip over any blank lines.
           */
          while (fileNames.size() <= 0) {
            std::getline(iStream, fileNames);
            fileNames = StringUtil::trim(fileNames);
          }
          m_fileNames = StringUtil::splitStrings(fileNames, " ");
      }

      /*!
	write checkpoint info to output stream

	\param oStream the output stream
      */
      void write(std::ostream &oStream) const
      {
        const char delim = '\n';
        oStream
	  << "V " << lsm_version_info::CheckPointVersion << delim
          << m_numTimeSteps << delim
          << m_timeStepSize << delim
          << m_timeStep     << delim;
	// ESyS-Particle version info
	oStream << lsm_version_info::ESySParticleVersion << m_bzrVersion << delim ;
        m_geoInfo.writeWithoutVersion(oStream);
        oStream
          << delim
          << StringUtil::join(
              m_fileNames.begin(),
              m_fileNames.end(),
              std::string(" "),
              StringUtil::StdOStreamOp<StringVector::const_iterator>()
            );
      }

      int          m_numTimeSteps;
      double       m_timeStepSize;
      int          m_timeStep;
      GeometryInfo m_geoInfo;
      StringVector m_fileNames;
      int          m_bzrVersion;
    };

    CheckPointInfo::CheckPointInfo() : m_pImpl(new CheckPointInfo::Impl)
    {
    }

    CheckPointInfo::~CheckPointInfo()
    {
      delete m_pImpl;
    }
    
    bool CheckPointInfo::operator==(const CheckPointInfo &cpInfo) const
    {
      return (*m_pImpl == *(cpInfo.m_pImpl));
    }

    const GeometryInfo &CheckPointInfo::getGeometryInfo() const
    {
      return m_pImpl->m_geoInfo;
    }

    void CheckPointInfo::setGeometryInfo(const GeometryInfo &geoInfo)
    {
      m_pImpl->m_geoInfo = geoInfo;
    }
    
    const StringVector &CheckPointInfo::getLatticeDataFiles() const
    {
      return m_pImpl->m_fileNames;
    }

    void CheckPointInfo::setLatticeDataFiles(const StringVector &fileNames)
    {
      m_pImpl->m_fileNames = fileNames;
    }
    
    void CheckPointInfo::setNumTimeSteps(int numTimeSteps)
    {
      m_pImpl->m_numTimeSteps = numTimeSteps;
    }

    void CheckPointInfo::setTimeStep(int timeStep)
    {
      m_pImpl->m_timeStep = timeStep;
    }

    int CheckPointInfo::getTimeStep() const
    {
      return m_pImpl->m_timeStep;
    }
    
    void CheckPointInfo::setTimeStepSize(double timeStepSize)
    {
      m_pImpl->m_timeStepSize = timeStepSize;
    }

    void CheckPointInfo::read(std::istream &iStream)
    {
      m_pImpl->read(iStream);
    }

    void CheckPointInfo::write(std::ostream &oStream) const
    {
      m_pImpl->write(oStream);
    }
  }
}
