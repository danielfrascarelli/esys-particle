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

#include "Geometry/ParticleCollection.h"
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <stdexcept>
#include <boost/limits.hpp>

namespace esys
{
  namespace lsm
  {
    template <typename TmplParticle>
    ParticleCollection<TmplParticle>::ParticleCollection()
      : m_particlePoolPtr(),
        m_particleVector()
    {
      m_particlePoolPtr = ParticlePoolPtr(new ParticlePool(2048));
    }

    template <typename TmplParticle>
    ParticleCollection<TmplParticle>::ParticleCollection(
      ParticlePoolPtr particlePoolPtr
    )
      : m_particlePoolPtr(particlePoolPtr),
        m_particleVector()
    {
    }

    template <typename TmplParticle>
    ParticleCollection<TmplParticle>::ParticleCollection(const ParticleCollection &p)
      : m_particlePoolPtr(p.m_particlePoolPtr),
        m_particleVector(p.m_particleVector)
    {
    }

    template <typename TmplParticle>
    ParticleCollection<TmplParticle> &
    ParticleCollection<TmplParticle>::operator=(const ParticleCollection &p)
    {
      m_particlePoolPtr = p.m_particlePoolPtr;
      m_particleVector  = p.m_particleVector;
    }

    template <typename TmplParticle>
    void
    ParticleCollection<TmplParticle>::noCheckInsertRef(Particle &p)
    {
      m_particleVector.push_back(&p);
    }

    template <typename TmplParticle>
    void
    ParticleCollection<TmplParticle>::insertRef(Particle &p)
    {
      if (m_particlePoolPtr->is_from(&p))
      {
        noCheckInsertRef(p);
      }
      else
      {
        throw
          std::runtime_error(
            "ParticleCollection<TmplParticle>::insertRef: Tried to insert"
            " reference to non-created particle."
          );
      }
    }

    template <typename TmplParticle>
    typename ParticleCollection<TmplParticle>::Particle &
    ParticleCollection<TmplParticle>::createParticle(const Particle &p)
    {
      Particle *newP = m_particlePoolPtr->construct(p);
      noCheckInsertRef(*newP);
      return *newP;
    }

    template <typename TmplParticle>
    ParticleCollection<TmplParticle>::~ParticleCollection()
    {
    }

    template <typename TmplParticle>
    int ParticleCollection<TmplParticle>::getNumParticles() const
    {
      return m_particleVector.size();
    }

    template <typename TmplParticle>
    void ParticleCollection<TmplParticle>::translateBy(const Vec3 &vec)
    {
      ParticleIterator it = getParticleIterator();
      while (it.hasNext())
      {
        it.next().translateBy(vec);
      }
    }

    template <typename TmplParticle>
    void ParticleCollection<TmplParticle>::rotate(
        const Vec3 &rotation,
        const Vec3 &posn
    )
    {
      ParticleIterator it = getParticleIterator();
      while (it.hasNext())
      {
        it.next().rotate(rotation, posn);
      }
    }

    template <typename TmplParticle>
    void ParticleCollection<TmplParticle>::incrementIdBy(
      typename Particle::Id idIncr
    )
    {
      ParticleIterator it = getParticleIterator();
      while (it.hasNext())
      {
        Particle &p = it.next();
        p.setId(p.getId() + idIncr);
      }
    }

    template <typename TmplParticle>
    BoundingBox ParticleCollection<TmplParticle>::getParticleBBox() const
    {
      Vec3 minPt = Vec3(std::numeric_limits<double>::max());
      Vec3 maxPt = -minPt;
      ParticleConstIterator it = getParticleIterator();
      while (it.hasNext())
      {
        const Particle &next = it.next();
        minPt = cmin(minPt, next.getPos() - next.getRad());
        maxPt = cmax(maxPt, next.getPos() + next.getRad());
      }
      return BoundingBox(minPt, maxPt);
    }
  }
}
