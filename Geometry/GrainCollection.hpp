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

#include "Geometry/GrainCollection.h"
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <boost/limits.hpp>

namespace esys
{
  namespace lsm
  {
    template <typename TmplGrain>
    GrainCollection<TmplGrain>::GrainCollection()
      : m_particlePoolPtr(new ParticlePool(4096)),
        m_grainPoolPtr(new GrainPool(4096)),
        m_grainVector()
    {
    }

    template <typename TmplGrain>
    GrainCollection<TmplGrain>::GrainCollection(
      ParticlePoolPtr particlePoolPtr
    )
      : m_particlePoolPtr(particlePoolPtr),
        m_grainPoolPtr(new GrainPool(4096)),
        m_grainVector()
    {
    }

    template <typename TmplGrain>
    GrainCollection<TmplGrain>::GrainCollection(
      ParticlePoolPtr particlePoolPtr,
      GrainPoolPtr grainPoolPtr
    )
      : m_particlePoolPtr(particlePoolPtr),
        m_grainPoolPtr(grainPoolPtr),
        m_grainVector()
    {
    }

    template <typename TmplGrain>
    typename GrainCollection<TmplGrain>::ParticlePoolPtr
    GrainCollection<TmplGrain>::getParticlePoolPtr()
    {
      return m_particlePoolPtr;
    }

    template <typename TmplGrain>
    typename GrainCollection<TmplGrain>::GrainPoolPtr
    GrainCollection<TmplGrain>::getGrainPoolPtr()
    {
      return m_grainPoolPtr;
    }

    template <typename TmplGrain>
    GrainCollection<TmplGrain>::~GrainCollection()
    {
    }

    template <typename TmplGrain>
    int GrainCollection<TmplGrain>::getNumGrains() const
    {
      return m_grainVector.size();
    }

    template <typename TmplGrain>
    int GrainCollection<TmplGrain>::getNumParticles() const
    {
      int i = 0;
      for (
        GrainConstIterator it = getGrainIterator();
        it.hasNext();
        i += it.next().getNumParticles()
      )
      {
      }
      return i;
    }

    template <typename TmplGrain>
    void
    GrainCollection<TmplGrain>::insertRef(Grain &g)
    {
      if (m_grainPoolPtr->is_from(&g))
      {
        m_grainVector.push_back(&g);
      }
      else
      {
        throw
          std::runtime_error(
            "GrainCollection<TmplGrain>::insertRef: Tried to insert"
            " reference to non-created grain."
          );
      }
    }

    template <typename TmplGrain>
    typename GrainCollection<TmplGrain>::Grain &
    GrainCollection<TmplGrain>::createGrain()
    {
      Grain *pGrain = m_grainPoolPtr->construct(getParticlePoolPtr());
      insertRef(*pGrain);
      return *pGrain;
    }

    template <typename TmplGrain>
    typename GrainCollection<TmplGrain>::Grain &
    GrainCollection<TmplGrain>::createGrain(typename Grain::Id id)
    {
      Grain *pGrain = m_grainPoolPtr->construct(id, getParticlePoolPtr());
      insertRef(*pGrain);
      return *pGrain;
    }

    template <typename TmplGrain>
    typename GrainCollection<TmplGrain>::Grain &
    GrainCollection<TmplGrain>::createGrain(const Grain &g)
    {
      Grain *pGrain = m_grainPoolPtr->construct(g);
      insertRef(*pGrain);
      return *pGrain;
    }

    template <typename TmplGrain>
    typename GrainCollection<TmplGrain>::GrainIterator
    GrainCollection<TmplGrain>::getGrainIterator()
    {
      return GrainIterator(VectorIterator(m_grainVector));
    }

    template <typename TmplGrain>
    typename GrainCollection<TmplGrain>::GrainConstIterator
    GrainCollection<TmplGrain>::getGrainIterator() const
    {
      return GrainConstIterator(VectorConstIterator(m_grainVector));
    }

  }
}
