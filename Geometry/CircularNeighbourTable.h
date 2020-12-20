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


#ifndef ESYS_LSMCIRCULARNEIGHBOURTABLE_H
#define ESYS_LSMCIRCULARNEIGHBOURTABLE_H

#include "Geometry/NeighbourTable.h"
#include <boost/pool/object_pool.hpp>
#include <boost/shared_ptr.hpp>

#include <sstream>
#include <stdexcept>
#include <set>

namespace esys
{
  namespace lsm
  {
    /**
     *
     */
    template <class TmplParticle>
    class CircularNeighbourTable : public NeighbourTable<TmplParticle>
    {
    public:
      typedef NeighbourTable<TmplParticle>                  Inherited;
      typedef typename Inherited::Particle                  Particle;
      typedef typename Inherited::ParticleVector            ParticleVector;
      typedef std::set<typename ParticleVector::value_type> ParticleSet;
      typedef boost::object_pool<Particle>                  ParticlePool;
      typedef boost::shared_ptr<ParticlePool>               ParticlePoolPtr;
      typedef std::vector<bool>                             BoolVector;

    public:
      CircularNeighbourTable(
        const BoundingBox  &bBox,
        double             gridSpacing,
        const BoolVector   &periodicDimensions = BoolVector(3, false),
        double             circBorderWidth = 0.0
      );

      CircularNeighbourTable(
        const BoundingBox  &bBox,
        double             gridSpacing,
        ParticlePoolPtr    particlePoolPtr,
        const BoolVector   &periodicDimensions = BoolVector(3, false),
        double             circBorderWidth = 0.0
      );

    public:
      void checkPeriodicDimensions();

      virtual ~CircularNeighbourTable();

      void setCircularBorderWidth(
        double circBorderWidth,
        double gridSpacing
      );
      
      void setCircularBorderWidth(double circBorderWidth);

      void resize(
        const BoundingBox &bBox,
        double gridSpacing,
        double circBorderWidth
      );

      void resize(
        const BoundingBox &bBox,
        double gridSpacing
      );
      
      void insertClone(Particle *pParticle, const Vec3 &newPosition);

      bool havePeriodicDimensions() const;

      Vec3 getModdedPosn(const Vec3 &posn) const;
      
      void insert(Particle *pParticle);

      void insert(Particle &particle);

      size_t getNumClonedParticles() const;

      size_t getNumParticles() const;

      const BoolVector &getPeriodicDimensions() const;

    protected:
      bool isClone(Particle *p) const;
      
      ParticleVector getNonClonedParticles();
      
      void clearClonedParticles();

    private:
      BoolVector      m_periodicDimensions;
      ParticlePoolPtr m_particlePoolPtr;
      ParticleSet     m_clonedParticleSet;
      int             m_circGridWidth;
      int             m_periodicDimIndex;
    };
  }
}

#include "Geometry/CircularNeighbourTable.hpp"

#endif
