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


#ifndef ESYS_LSMPACKER_H
#define ESYS_LSMPACKER_H

#include "Geometry/CircularNeighbourTable.h"
#include "Geometry/ParticleCollection.h"
#include <boost/pool/object_pool.hpp>
#include <boost/shared_ptr.hpp>

#include <set>

namespace esys
{
  namespace lsm
  {
    /**
     *
     */
    template <typename TmplParticleCollection>
    class Packer
    {
    public:
      typedef TmplParticleCollection                     ParticleCollection;
      typedef Packer<ParticleCollection>                 PackerBase;
      typedef boost::shared_ptr<ParticleCollection>      ParticleCollectionPtr;
      typedef typename ParticleCollection::Particle      Particle;
      typedef boost::object_pool<Particle>               ParticlePool;
      typedef boost::shared_ptr<ParticlePool>            ParticlePoolPtr;
      typedef CircularNeighbourTable<Particle>           NTable;
      typedef boost::shared_ptr<NTable>                  NTablePtr;
      typedef
        typename ParticleCollection::ParticleIterator
        ParticleIterator;
      typedef
        typename ParticleCollection::ParticleConstIterator
        ParticleConstIterator;

      Packer(NTablePtr nTablePtr);

      Packer(ParticlePoolPtr particlePoolPtr, NTablePtr nTablePtr);

      virtual ~Packer();

      virtual void generate() = 0;

      int getNumParticles() const;

      int getNextParticleId();
      
      void setNTablePtr(NTablePtr nTablePtr);
      NTable &getNTable();
      const NTable &getNTable() const;

      ParticlePoolPtr getParticlePoolPtr();
      ParticlePool &getParticlePool();
      const ParticlePool &getParticlePool() const;

      ParticleCollection &getParticleCollection();
      const ParticleCollection &getParticleCollection() const;

      Particle &constructParticle(const Particle &particle);

      ParticleIterator getParticleIterator();

      ParticleConstIterator getParticleIterator() const;

      bool contains(const Particle &particle) const;

      Particle &createAndInsertParticle(const Particle &particle);

    protected:

      typedef std::set<int> IdSet;

    private:
      NTablePtr             m_nTablePtr;
      ParticlePoolPtr       m_particlePoolPtr;
      ParticleCollectionPtr m_particleCollectionPtr;
      IdSet                 m_idSet;
    };
  }
}

#include "Geometry/Packer.hpp"

#endif
