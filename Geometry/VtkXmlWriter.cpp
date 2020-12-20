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


#include "Geometry/VtkXmlWriter.h"

#include <vector>
#include <map>
#include <sstream>
#include <fstream>
#include <stdexcept>

namespace esys
{
  namespace lsm
  {
    typedef std::vector<const SimpleParticle *>   ParticleVector;
    typedef std::vector<const BasicInteraction *> ConnectionVector;
    typedef std::map<int,int> IdIndexMap;
  //=============================================================================
    class ParticleDataVisitor::Impl
    {
    public:
      Impl()
      {
      }

      void append(const SimpleParticle &particle)
      {
        m_idIndexMap[particle.getID()] = m_particleVector.size();
        m_particleVector.push_back(&particle);
      }

      void append(const BasicInteraction &connection)
      {
        m_connectionVector.push_back(&connection);
      }

      ParticleVector    m_particleVector;
      ConnectionVector m_connectionVector;
      IdIndexMap        m_idIndexMap;
    };

  //=============================================================================
    ParticleDataVisitor::ParticleDataVisitor()
      : m_implPtr(new Impl())
    {
    }

    void ParticleDataVisitor::visitSimpleParticle(const SimpleParticle &particle)
    {
      m_implPtr->append(particle);
    }

    void ParticleDataVisitor::visitParticle(const Particle &particle)
    {
      visitSimpleParticle(particle);
    }

    void ParticleDataVisitor::visitBasicInteraction(const Connection &connection)
    {
      m_implPtr->append(connection);
    }

    void ParticleDataVisitor::visitConnection(const Connection &connection)
    {
      visitBasicInteraction(connection);
    }

    size_t ParticleDataVisitor::getNumParticles() const
    {
      return m_implPtr->m_particleVector.size();
    }

    size_t ParticleDataVisitor::getNumConnections() const
    {
      return m_implPtr->m_connectionVector.size();
    }

    template <class TmplContainer>
    class ConstContainerIterator
    {
    public:
      ConstContainerIterator(const TmplContainer &container)
        : m_it(container.begin()),
          m_end(container.end())
      {
      }

      bool hasNext() const
      {
        return (m_it != m_end);
      }

      typename TmplContainer::const_reference next()
      {
        typename TmplContainer::const_reference valueRef = (*m_it);
        ++m_it;
        return valueRef;
      }

    private:
      typename TmplContainer::const_iterator m_it;
      typename TmplContainer::const_iterator m_end;
    };

    class ParticleIterator
    {
    public:
      ParticleIterator(const ParticleVector &particleVector)
        : m_it(particleVector)
      {
      }

      bool hasNext() const
      {
        return m_it.hasNext();
      }

      const SimpleParticle &next()
      {
        return *m_it.next();
      }

    private:
      ConstContainerIterator<ParticleVector> m_it;
    };

    class ConnectionIterator
    {
    public:
      ConnectionIterator(const ConnectionVector &connectionVector)
        : m_it(connectionVector)
      {
      }

      bool hasNext() const
      {
        return m_it.hasNext();
      }

      const BasicInteraction &next()
      {
        return *m_it.next();
      }

    private:
      ConstContainerIterator<ConnectionVector> m_it;
    };

    void ParticleDataVisitor::writeCentrePoints(std::ostream &oStream) const
    {
      ParticleIterator it(m_implPtr->m_particleVector);
      while (it.hasNext())
      {
        oStream << it.next().getPos() << "\n";
      }
    }

    void ParticleDataVisitor::writeRadii(std::ostream &oStream) const
    {
      ParticleIterator it(m_implPtr->m_particleVector);
      while (it.hasNext())
      {
        oStream << it.next().getRad() << "\n";
      }
    }

    void ParticleDataVisitor::writeTags(std::ostream &oStream) const
    {
      ParticleIterator it(m_implPtr->m_particleVector);
      while (it.hasNext())
      {
        oStream << it.next().getTag() << "\n";
      }
    }

    void ParticleDataVisitor::writeIds(std::ostream &oStream) const
    {
      ParticleIterator it(m_implPtr->m_particleVector);
      while (it.hasNext())
      {
        oStream << it.next().getID() << "\n";
      }
    }

    int ParticleDataVisitor::getIndex(int particleId) const
    {
      IdIndexMap::const_iterator it = m_implPtr->m_idIndexMap.find(particleId);
      if (it == m_implPtr->m_idIndexMap.end())
      {
        std::stringstream msg;
        msg << "Could not find particle id " << particleId << " in index map.";
        throw std::runtime_error(msg.str());
      }
      return it->second;
    }

    void ParticleDataVisitor::writeParticleIndexConnections(std::ostream &oStream) const
    {
      ConnectionIterator it(m_implPtr->m_connectionVector);
      while (it.hasNext())
      {
        const BasicInteraction &connection = it.next();
        oStream 
          << getIndex(connection.first())
          << " "
          << getIndex(connection.second())
          << "\n";
      }
    }

    void ParticleDataVisitor::writeConnectionTags(std::ostream &oStream) const
    {
      ConnectionIterator it(m_implPtr->m_connectionVector);
      while (it.hasNext())
      {
        oStream 
          << it.next().getTag()
          << "\n";
      }
    }

  //=============================================================================
    class VtkXmlWriter::Impl
    {
    public:
      Impl() : m_pParticleData(NULL)
      {
      }
      
      const ParticleDataVisitor *m_pParticleData;
    };

  //=============================================================================
    VtkXmlWriter::VtkXmlWriter()
      : m_implPtr(new Impl())
    {
    }

    VtkXmlWriter::~VtkXmlWriter()
    {
    }
    
    void VtkXmlWriter::setData(const ParticleDataVisitor &particleData)
    {
      m_implPtr->m_pParticleData = &particleData;
    }
    
    size_t VtkXmlWriter::getNumParticles() const
    {
      return m_implPtr->m_pParticleData->getNumParticles();
    }

    size_t VtkXmlWriter::getNumConnections() const
    {
      return m_implPtr->m_pParticleData->getNumConnections();
    }

    void VtkXmlWriter::writePoints(std::ostream &oStream)
    {
      oStream << "<Points>\n";
      oStream << "<DataArray NumberOfComponents=\"3\" type=\"Float32\" format=\"ascii\">\n";
      m_implPtr->m_pParticleData->writeCentrePoints(oStream);
      oStream << "</DataArray>\n";
      oStream << "</Points>\n";
    }

    void VtkXmlWriter::writePointData(std::ostream &oStream)
    {
      oStream << "<PointData Scalars=\"radius\">\n";

      oStream << "<DataArray type=\"Float32\" Name=\"radius\" NumberOfComponents=\"1\" format=\"ascii\">\n";
      m_implPtr->m_pParticleData->writeRadii(oStream);
      oStream << "</DataArray>\n";
      oStream << "<DataArray type=\"Int32\" Name=\"particleTag\" NumberOfComponents=\"1\" format=\"ascii\">\n";
      m_implPtr->m_pParticleData->writeTags(oStream);
      oStream << "</DataArray>\n";
      oStream << "<DataArray type=\"Int32\" Name=\"Id\" NumberOfComponents=\"1\" format=\"ascii\">\n";
      m_implPtr->m_pParticleData->writeIds(oStream);
      oStream << "</DataArray>\n";

      oStream << "</PointData>\n";
    }

    void VtkXmlWriter::writeCells(std::ostream &oStream)
    {
      oStream << "<Cells>\n";
      oStream << "<DataArray type=\"Int32\" NumberOfComponents=\"1\" Name=\"connectivity\" format=\"ascii\">\n";
      m_implPtr->m_pParticleData->writeParticleIndexConnections(oStream);
      oStream << "</DataArray>";
      
      oStream << "<DataArray type=\"Int32\" NumberOfComponents=\"1\" Name=\"offsets\" format=\"ascii\">\n";
      for (size_t i = 1; i < getNumConnections()*2; i+=2) oStream << i+1 << "\n";
      oStream << "</DataArray>\n";

      oStream << "<DataArray type=\"UInt8\" NumberOfComponents=\"1\" Name=\"types\" format=\"ascii\">\n";
      const int CELL_LINE_TYPE = 3;
      for (size_t i = 0; i < getNumConnections(); i++) oStream << CELL_LINE_TYPE << "\n";
      oStream << "</DataArray>\n";

      oStream << "</Cells>\n";
    }

    void VtkXmlWriter::writeCellData(std::ostream &oStream)
    {
      oStream << "<CellData>\n";
      oStream << "<DataArray type=\"Int32\" NumberOfComponents=\"1\" Name=\"bondTag\" format=\"ascii\">\n";
      m_implPtr->m_pParticleData->writeConnectionTags(oStream);
      oStream << "</DataArray>\n";
      oStream << "</CellData>\n";
    }

    void VtkXmlWriter::write(std::ostream &oStream)
    {
      oStream << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
      oStream << "<UnstructuredGrid>\n";
      oStream 
        << "<Piece NumberOfPoints=\"" << getNumParticles()
        << "\" NumberOfCells=\"" << getNumConnections() << "\">\n";
      writePoints(oStream);
      writePointData(oStream);
      writeCells(oStream);
      writeCellData(oStream);
      oStream << "</Piece>\n";
      oStream << "</UnstructuredGrid>\n";
      oStream << "</VTKFile>\n";
    }

    void VtkXmlWriter::writeToFile(const std::string &fileName)
    {
      std::ofstream oStream(fileName.c_str());
      oStream << "<?xml version=\"1.0\"?>\n";
      write(oStream);
    }
  }
}
