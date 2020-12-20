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

#ifndef ESYS_LSMGRAIN_H
#define ESYS_LSMGRAIN_H

#include "Geometry/ParticleCollection.h"

namespace esys
{
  namespace lsm
  {
    /**
     *
     */
    template <typename TmplParticleCollection>
    class Grain : public TmplParticleCollection
    {
    public:
      typedef TmplParticleCollection                    Inherited;
      typedef Inherited                                 ParticleCollection;
      typedef typename Inherited::Particle              Particle;
      typedef typename Inherited::ParticlePool          ParticlePool;
      typedef typename Inherited::ParticlePoolPtr       ParticlePoolPtr;
      typedef typename Inherited::ParticleIterator      ParticleIterator;
      typedef typename Inherited::ParticleConstIterator ParticleConstIterator;
      typedef long Id;

      Grain();

      Grain(Id id);

      Grain(ParticlePoolPtr particlePoolPtr);

      Grain(Id id, ParticlePoolPtr particlePoolPtr);

      Grain(const Grain &g);

      Grain &operator=(const Grain &g);

      Id getId() const;

      void setId(Id id);

      void setParticleIds(typename Particle::Id minId);

    private:
      Id m_id;
    };
  }
}

#include "Geometry/Grain.hpp"

#endif
