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

#ifndef ESYS_LSMCHECKPOINTLOADER_H
#define ESYS_LSMCHECKPOINTLOADER_H

#include "Model/BondedInteraction.h"
#include "Parallel/IterativeReader.h"
#include "Model/Particle.h"
#include "Model/BondedInteractionCpData.h"

#include <vector>
#include <string>
#include <fstream>

namespace esys
{
  namespace lsm
  {
    typedef std::vector<std::string> StringVector;

    /**
      * Objects of this class load particle and interaction check-point
      * data from file and initialise a specified object via the
      * CheckPointLoader::loadInto method.
      */
    class CheckPointLoader
    {
    public:
      class ParticleData : public CParticle
      {
      public:
        ParticleData() : CParticle()
        {
        }

        void read(std::istream &iStream)
        {
          loadCheckPointData(iStream);
        }
      };

      class ConnectionData : public BondedInteractionCpData
      {
      public:
        ConnectionData() : BondedInteractionCpData()
        {
        }
        
        void read(std::istream &iStream)
        {
          loadCheckPointData(iStream);
        }
      };
      
      class ParticleReader : public IterativeReader<IStreamIterator<ParticleData> >
      {
      public:
        typedef IterativeReader<IStreamIterator<ParticleData> >::Iterator Iterator;
        
        ParticleReader(std::istream &iStream) : IterativeReader<IStreamIterator<ParticleData> >(iStream)
        {
        }
        
        virtual void initialise()
        {
          int numParticles = 0;
          getIStream() >> numParticles;
          setNumElements(numParticles);
          IterativeReader<IStreamIterator<ParticleData> >::initialise();
        }
      };

      class ConnectionReader : public IterativeReader<IStreamIterator<ConnectionData> >
      {
      public:
        ConnectionReader(std::istream &iStream) : IterativeReader<IStreamIterator<ConnectionData> >(iStream)
        {
        }

        virtual void initialise()
        {
          int numConnections = 0;
          getIStream() >> numConnections;
          setNumElements(numConnections);
          IterativeReader<IStreamIterator<ConnectionData> >::initialise();
        }
      };

      CheckPointLoader(const StringVector &fileNames) : m_fileNames(fileNames)
      {
      }

      template<class TmplLsmData>
      void loadInto(TmplLsmData &lsmData)
      {
        for (
          StringVector::const_iterator it = m_fileNames.begin();
          it != m_fileNames.end();
          it++
        )
        {
          std::ifstream iStream(it->c_str());
          ParticleReader pReader(iStream);
          lsmData.template addParticles<ParticleReader::Iterator,CParticle>(pReader.getIterator());

          int numConnectionGroups = 0;
          iStream >> numConnectionGroups;
          for (int i = 0; i < numConnectionGroups; i++) {
            ConnectionReader cReader(iStream);
            lsmData.addConnections(cReader.getIterator());
          }
        }
      }

    private:
      StringVector m_fileNames;
    };
  }
}

#endif
