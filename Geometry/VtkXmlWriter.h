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


#ifndef ESYS_LSMVTKXMLWRITER_H
#define ESYS_LSMVTKXMLWRITER_H

#include <Geometry/SimpleParticle.h>
#include <Geometry/BasicInteraction.h>

#include <iostream>
#include <boost/shared_ptr.hpp>

namespace esys
{
  namespace lsm
  {
    class ParticleDataVisitor
    {
      typedef SimpleParticle    Particle;
      typedef BasicInteraction  Connection;
    public:
      ParticleDataVisitor();

      void visitSimpleParticle(const Particle &particle);

      void visitParticle(const Particle &particle);

      void visitBasicInteraction(const Connection &connection);

      void visitConnection(const Connection &connection);

      size_t getNumParticles() const;

      size_t getNumConnections() const;

      int getIndex(int particleId) const;

      void writeCentrePoints(std::ostream &oStream) const;

      void writeRadii(std::ostream &oStream) const;

      void writeTags(std::ostream &oStream) const;

      void writeIds(std::ostream &oStream) const;

      void writeParticleIndexConnections(std::ostream &oStream) const;
      
      void writeConnectionTags(std::ostream &oStream) const;

    private:
      class Impl;
      typedef boost::shared_ptr<Impl> ImplPtr;
      ImplPtr m_implPtr;
    };

    /**
     *
     */
    class VtkXmlWriter
    {
    public:
      VtkXmlWriter();

      virtual ~VtkXmlWriter();

      void setData(const ParticleDataVisitor &particleData);

      size_t getNumParticles() const;

      size_t getNumConnections() const;
      
      virtual void writePoints(std::ostream &oStream);
      virtual void writePointData(std::ostream &oStream);
      virtual void writeCells(std::ostream &oStream);
      virtual void writeCellData(std::ostream &oStream);

      virtual void write(std::ostream &oStream);
      
      virtual void writeToFile(const std::string &fileName);

      private:

      class Impl;
      typedef boost::shared_ptr<Impl> ImplPtr;
      ImplPtr m_implPtr;
    };
  }
}

#endif
