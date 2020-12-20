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


#ifndef ESYS_LSMBLOCKGENERATOR_H
#define ESYS_LSMBLOCKGENERATOR_H

#include <Geometry/ParticleGenerator.h>
#include <Geometry/SimpleParticle.h>
#include <Foundation/BoundingBox.h>

#include <vector>
#include <set>

namespace esys
{
  namespace lsm
  {
    typedef std::vector<bool> BoolVector;
    /**
     *
     */
    class BlockGenerator : public ParticleGenerator
    {
    public:
      BlockGenerator(
        NTable            &nTable,
        ParticlePool      &particlePool,
        const BoundingBox &bBox,
        const BoolVector  &periodicDimensions,
        double            tolerance
      );

      virtual ~BlockGenerator();

      virtual void generate() = 0;

      virtual void generateSeedParticles();

      virtual SimpleParticle generateParticle(const Vec3 &point);

      virtual double getRadius() const = 0;
      
      virtual double getGridRadius() const = 0;

      size_t getNumParticles() const;

      int getNextId();

      virtual bool particleFits(const SimpleParticle &particle) const;

      bool is2d() const;

      bool particleFitsInBBox(const SimpleParticle &particle) const;

      bool particleFitsWithNeighbours(const SimpleParticle &particle) const;

      void insertParticle(const SimpleParticle &particle);

      double getTolerance() const;

      const BoundingBox &getBBox() const;

      bool contains(const SimpleParticle &particle) const;

      typedef NTable::ParticleVector   ParticleVector;
      typedef NTable::ParticleIterator ParticleIterator;
      
      ParticleIterator getParticleIterator()
      {
        return ParticleIterator(m_particleVector);
      }

      typedef std::set<int> IdSet;
    private:
      BoundingBox    m_bBox;
      BoolVector     m_periodicDimensions;
      ParticleVector m_particleVector;
      double         m_tolerance;
      IdSet          m_idSet;
    };
  }
}

#endif
