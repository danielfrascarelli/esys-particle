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

#ifndef ESYS_LSMGRAINCOLLECTION_H
#define ESYS_LSMGRAINCOLLECTION_H

#include "Foundation/StlIterator.h"
#include <boost/shared_ptr.hpp>
#include <boost/pool/object_pool.hpp>

#include <vector>

namespace esys
{
  namespace lsm
  {
    /**
     *
     */
    template <typename TmplGrain>
    class GrainCollection
    {
    public:
      typedef TmplGrain                                    Grain;
      typedef typename Grain::Particle                     Particle;
      typedef typename Grain::ParticleCollection           ParticleCollection;
      typedef typename ParticleCollection::ParticlePool    ParticlePool;
      typedef typename ParticleCollection::ParticlePoolPtr ParticlePoolPtr;
      typedef typename ParticleCollection::ParticleIterator      ParticleIterator;
      typedef typename ParticleCollection::ParticleConstIterator ParticleConstIterator;
      typedef boost::object_pool<Grain>                    GrainPool;
      typedef boost::shared_ptr<GrainPool>                 GrainPoolPtr;
      typedef std::vector<Grain *>                         GrainVector;

    protected:
      typedef ForwardIterator<GrainVector>                 VectorIterator;
      typedef ForwardConstIterator<GrainVector>            VectorConstIterator;

    public:

      class GrainIterator : public VectorIterator
      {
      public:
        typedef Grain& value_type;
        GrainIterator(const VectorIterator &it)
         : VectorIterator(it)
        {
        }

        value_type next()
        {
          return *(VectorIterator::next());
        }

        value_type current() const
        {
          return *(VectorIterator::current());
        }
      };

      class GrainConstIterator : public VectorConstIterator
      {
      public:
        typedef const Grain& value_type;
        GrainConstIterator (const VectorConstIterator &it)
         : VectorConstIterator(it)
        {
        }

        GrainConstIterator (const VectorIterator &it)
         : VectorConstIterator(it)
        {
        }

        value_type next()
        {
          return *(VectorConstIterator::next());
        }

        value_type current() const
        {
          return *(VectorConstIterator::current());
        }
      };

      GrainCollection();

      GrainCollection(ParticlePoolPtr particlePoolPtr);

      GrainCollection(ParticlePoolPtr particlePoolPtr, GrainPoolPtr grainPoolPtr);

      virtual ~GrainCollection();

      /**
       * Returns the number of grains in this collection.
       */
      int getNumGrains() const;

      /**
       * Returns the number of particles contained in all
       * grains of this collection.
       */
      int getNumParticles() const;

      /**
       * Stores reference to specified grain.
       *
       * @param g Inserts reference to grain g in this collection.
       * @throws std::runtime_error if g was not created by this
       *   collection's GrainPool.
       */
      void insertRef(Grain &g);
      
      /**
       * Creates an empty grain.
       * @return reference to new grain.
       */
      Grain &createGrain();

      /**
       * Creates an empty (no particles) grain.
       * @param id Create a grain with this id.
       * @return reference to new grain.
       */
      Grain &createGrain(typename Grain::Id id);

      /**
       * Returns a copy-constructed grain.
       * @param g Copy created from this argument.
       * @return reference to new grain.
       */
      Grain &createGrain(const Grain &g);

      GrainIterator getGrainIterator();

      GrainConstIterator getGrainIterator() const;

    protected:
      ParticlePoolPtr getParticlePoolPtr();

      GrainPoolPtr getGrainPoolPtr();

    private:
      ParticlePoolPtr m_particlePoolPtr;
      GrainPoolPtr    m_grainPoolPtr;
      GrainVector     m_grainVector;
    };
  }
}

#include "Geometry/GrainCollection.hpp"

#endif
