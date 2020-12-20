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


#include "Parallel/GeometryReader.h"
#include "Foundation/console.h"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace esys
{
  namespace lsm
  {
//===============================================================================

    ParticleReader::ParticleReader(std::istream &iStream, bool is2d)
      : IterativeReader<ParticleIterator>(iStream),
        m_is2d(is2d)
    {
    }

    void ParticleReader::initialise()
    {
      std::string token;
      while ((!getIStream().eof()) && (token != "BeginParticles")) {
        getIStream() >> token;
      }
      if (token != "BeginParticles") {
        throw std::runtime_error("Could not find 'BeginParticles' token in stream.");
      }

      getIStream() >> m_particleType;

      int numParticles = 0;
      getIStream() >> numParticles;

      setNumElements(numParticles);

      IterativeReader<ParticleIterator>::initialise();
    }

    ParticleIterator *ParticleReader::createNewIterator()
    {
      return new ParticleIterator(getIStream(), getNumElements(), m_is2d);
    }
            
    const std::string &ParticleReader::getParticleType()
    {
      if (!this->isInitialised()) {
        initialise();
      }
      return m_particleType;
    }
    

//===============================================================================
    SimpleConnectionData::SimpleConnectionData()
      : m_particle1Id(-1),
        m_particle2Id(-1),
        m_tag(-1)
    {
    }

    SimpleConnectionData::SimpleConnectionData(Id p1Id, Id p2Id, Tag tag)
      : m_particle1Id(p1Id),
        m_particle2Id(p2Id),
        m_tag(tag)
    {
    }

    bool SimpleConnectionData::operator==(
      const SimpleConnectionData &particleData
    ) const
    {
      return
        (
          (getP1Id()       == particleData.getP1Id())
          &&
          (getP2Id()       == particleData.getP2Id())
          &&
          (getTag()      == particleData.getTag())
        );
    }

    const SimpleConnectionData::Id &SimpleConnectionData::getP1Id() const
    {
      return m_particle1Id;
    }

    const SimpleConnectionData::Id &SimpleConnectionData::getP2Id() const
    {
      return m_particle2Id;
    }

    const SimpleConnectionData::Tag &SimpleConnectionData::getTag() const
    {
      return m_tag;
    }

    void SimpleConnectionData::read(std::istream &istream)
    {
      istream
        >> m_particle1Id
        >> m_particle2Id
        >> m_tag;
    }

    void SimpleConnectionData::write(std::ostream &oStream) const
    {
      oStream
        << getP1Id() << " "
        << getP2Id() << " "
        << getTag();
    }

    std::istream &operator>>(std::istream &iStream, SimpleConnectionData &connectionData)
    {
      connectionData.read(iStream);
      return iStream;
    }

    std::ostream &operator<<(std::ostream &oStream, const SimpleConnectionData &connectionData)
    {
      connectionData.write(oStream);
      return oStream;
    }

//===============================================================================

    ConnectionReader::ConnectionReader(std::istream &iStream)
      : IterativeReader<IStreamIterator<SimpleConnectionData> >(iStream)
    {
    }

    void ConnectionReader::initialise()
    {
      std::string token;
      while ((!getIStream().eof()) && (token != "BeginConnect")) {
        getIStream() >> token;
      }
      if (token != "BeginConnect") {
        throw std::runtime_error("Could not find 'BeginConnect' token in stream.");
      }
      
      int numConnections = 0;
      getIStream() >> numConnections;

      setNumElements(numConnections);
      
      IterativeReader<IStreamIterator<SimpleConnectionData> >::initialise();
    }

//===============================================================================    
    class GeometryReader::Impl
    {
    public:
      Impl(const std::string &fileName);

      Impl(std::istream &iStream);

      ~Impl();

      void initialise();
      
      void initialiseFile();
      
      void initialiseStream();

      typedef std::auto_ptr<ParticleReader>   ParticleReaderPtr;
      typedef std::auto_ptr<ConnectionReader> ConnectionReaderPtr;
      typedef std::auto_ptr<std::istream>     IStreamPtr;

      std::string         m_fileName;
      GeometryInfo        m_geoInfo;
      IStreamPtr          m_iStreamPtr;
      std::istream        *m_pIStream;
      ParticleReaderPtr   m_particleReaderPtr;
      ConnectionReaderPtr m_connectionReaderPtr;
    };

    GeometryReader::Impl::Impl(const std::string &fileName)
      : m_fileName(fileName),
        m_geoInfo(),
        m_iStreamPtr(),
        m_pIStream(NULL),
        m_particleReaderPtr(),
        m_connectionReaderPtr()
    {
    }

    GeometryReader::Impl::Impl(std::istream &iStream)
      : m_fileName(),
        m_geoInfo(),
        m_iStreamPtr(),
        m_pIStream(&iStream),
        m_particleReaderPtr(),
        m_connectionReaderPtr()
    {
    }

    void GeometryReader::Impl::initialiseStream()
    {
      m_geoInfo.read(*m_pIStream);

      m_particleReaderPtr   = ParticleReaderPtr(new ParticleReader(*m_pIStream, m_geoInfo.is2d()));
      m_connectionReaderPtr = ConnectionReaderPtr(new ConnectionReader(*m_pIStream));
    }

    void GeometryReader::Impl::initialiseFile()
    {
      m_iStreamPtr = IStreamPtr(new std::ifstream(m_fileName.c_str()));
      m_pIStream = m_iStreamPtr.get();
      
      // check if file is open
      if(!(*m_pIStream)){
        throw std::runtime_error("Can not open geometry description file " + m_fileName);
      }
      console.Debug() << "Reading geometry file " << m_fileName << "  open \n";
      initialiseStream();
    }

    void GeometryReader::Impl::initialise()
    {
      if (m_pIStream == NULL) {
        initialiseFile();
      }
      else
      {
        initialiseStream();
      }
    }

    GeometryReader::Impl::~Impl()
    {
    }

    GeometryReader::GeometryReader(const std::string &fileName)
      : m_pImpl(new Impl(fileName))
    {
      m_pImpl->initialise();
    }

    GeometryReader::GeometryReader(std::istream &iStream)
      : m_pImpl(new Impl(iStream))
    {
      m_pImpl->initialise();
    }

    GeometryReader::~GeometryReader()
    {
      delete m_pImpl;
    }

    const GeometryInfo &GeometryReader::getGeometryInfo() const
    {
      return m_pImpl->m_geoInfo;
    }

    const std::string &GeometryReader::getFileName() const
    {
      return m_pImpl->m_fileName;
    }

    const std::string &GeometryReader::getParticleType()
    {
      return m_pImpl->m_particleReaderPtr->getParticleType();
    }

    GeometryReader::ParticleIterator &GeometryReader::getParticleIterator()
    {
      return m_pImpl->m_particleReaderPtr->getIterator();
    }

    GeometryReader::ConnectionIterator &GeometryReader::getConnectionIterator()
    {
      return m_pImpl->m_connectionReaderPtr->getIterator();
    }
  } // namespace lsm
} // namespace esys
