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


#include "Geometry/Grain.h"
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <boost/limits.hpp>

namespace esys
{
  namespace lsm
  {
    template <typename TmplParticleCollection>
    Grain<TmplParticleCollection>::Grain() : Inherited(), m_id(-1)
    {
    }

    template <typename TmplParticleCollection>
    Grain<TmplParticleCollection>::Grain(ParticlePoolPtr particlePoolPtr)
      : Inherited(particlePoolPtr),
        m_id(-1)
    {
    }

    template <typename TmplParticleCollection>
    Grain<TmplParticleCollection>::Grain(Id id) : Inherited(), m_id(id)
    {
    }

    template <typename TmplParticleCollection>
    Grain<TmplParticleCollection>::Grain(Id id, ParticlePoolPtr particlePoolPtr)
      : Inherited(particlePoolPtr),
        m_id(id)
    {
    }

    template <typename TmplParticleCollection>
    Grain<TmplParticleCollection>::Grain(const Grain &g)
      : Inherited(g), m_id(g.getId())
    {
    }

    template <typename TmplParticleCollection>
    Grain<TmplParticleCollection> &
    Grain<TmplParticleCollection>::operator=(const Grain &g)
    {
      Inherited::operator=(g);
      setId(g.getId());
      return *this;
    }

    template <typename TmplParticleCollection>
    typename Grain<TmplParticleCollection>::Id
    Grain<TmplParticleCollection>::getId() const
    {
      return m_id;
    }
    
    template <typename TmplParticleCollection>
    void Grain<TmplParticleCollection>::setId(Id id)
    {
      m_id = id;
    }
    
    template <typename TmplParticleCollection>
    void Grain<TmplParticleCollection>::setParticleIds(
      typename Particle::Id minId
    )
    {
      ParticleIterator it = this->getParticleIterator();
      while (it.hasNext())
      {
        it.next().setId(minId++);
      }
    }
  }
}
