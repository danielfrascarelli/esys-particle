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


#include "Geometry/BlockGenerator.h"
#include "Geometry/GridIterator.h"
#include "Geometry/SimpleParticle.h"

#include <float.h>

namespace esys
{
  namespace lsm
  {
    BlockGenerator::BlockGenerator(
      NTable            &nTable,
      ParticlePool      &particlePool,
      const BoundingBox &bBox,
      const BoolVector  &periodicDimensions,
      double            tolerance
    ) : ParticleGenerator(nTable, particlePool),
        m_bBox(bBox),
        m_periodicDimensions(periodicDimensions),
        m_tolerance(tolerance)
    {
    }

    BlockGenerator::~BlockGenerator()
    {
    }
    
    const BoundingBox &BlockGenerator::getBBox() const
    {
      return m_bBox;
    }

    size_t BlockGenerator::getNumParticles() const
    {
      return m_idSet.size();
    }

    int BlockGenerator::getNextId()
    {
      return static_cast<int>(getNTable().getNumParticles());
    }

    double BlockGenerator::getTolerance() const
    {
      return m_tolerance;
    }

    bool BlockGenerator::is2d() const
    {
      return ((getBBox().getMaxPt().Z() - getBBox().getMinPt().Z()) <= 0);
    }

    bool BlockGenerator::particleFitsInBBox(const SimpleParticle &particle) const
    {
      return
        (
          (
            m_periodicDimensions[0]
            ||
            (
              m_bBox.contains(particle.getPos() - Vec3(particle.getRad(), 0, 0), getTolerance())
              &&
              m_bBox.contains(particle.getPos() + Vec3(particle.getRad(), 0, 0), getTolerance())
            )
          )
          &&
          (
            m_periodicDimensions[1]
            ||
            (
              m_bBox.contains(particle.getPos() - Vec3(0, particle.getRad(), 0), getTolerance())
              &&
              m_bBox.contains(particle.getPos() + Vec3(0, particle.getRad(), 0), getTolerance())
            )
          )
          &&
          (
            is2d() || m_periodicDimensions[2]
            ||
            (
              m_bBox.contains(particle.getPos() - Vec3(0, 0, particle.getRad()), getTolerance())
              &&
              m_bBox.contains(particle.getPos() + Vec3(0, 0, particle.getRad()), getTolerance())
            )
          )
        );
    }

    bool BlockGenerator::particleFitsWithNeighbours(const SimpleParticle &particle) const
    {
      const ParticleVector neighbours = 
        getNTable().getNeighbourVector(
          particle.getPos(),
          particle.getRad() + getTolerance()
        );
      ParticleVector::const_iterator iter = neighbours.begin();
      for (; iter != neighbours.end(); iter++) {
        const double dist = (particle.getPos() - (*iter)->getPos()).norm();
        if (dist < ((particle.getRad() + (*iter)->getRad()) - getTolerance())) {
          return false;
        }
      }
      return true;
    }

    bool BlockGenerator::particleFits(const SimpleParticle &particle) const
    {
      return (particleFitsInBBox(particle) && particleFitsWithNeighbours(particle));
    }

    void BlockGenerator::insertParticle(const SimpleParticle &particle)
    {
      SimpleParticle *pParticle = getParticlePool().construct(particle);
      m_particleVector.push_back(pParticle);
      m_idSet.insert(pParticle->getID());
      getNTable().insert(pParticle);
    }

    bool BlockGenerator::contains(const SimpleParticle &particle) const
    {
      return (m_idSet.find(particle.getID()) != m_idSet.end());
    }

    SimpleParticle BlockGenerator::generateParticle(const Vec3 &point)
    {
      return SimpleParticle(point, getRadius(), getNextId());
    }

    void BlockGenerator::generateSeedParticles()
    {
      GridIterator pointIt = GridIterator(getBBox(), getGridRadius());
      while (pointIt.hasNext()) {
        SimpleParticle particle = generateParticle(pointIt.next());
        if (particleFits(particle)) {
          insertParticle(particle);
        }
      }
    }

  }
}
