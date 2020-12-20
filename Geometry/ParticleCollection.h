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

#ifndef ESYS_LSMPARTICLECOLLECTION_H
#define ESYS_LSMPARTICLECOLLECTION_H

#include "Foundation/BoundingBox.h"
#include "Foundation/StlIterator.h"
#include "Geometry/Vec3L.h"

#include <boost/shared_ptr.hpp>
#include <boost/pool/object_pool.hpp>

#include <vector>
#include <float.h>

namespace esys
{
  namespace lsm
  {
    /**
     *
     */
    template <typename TmplParticle>
    class ParticleCollection
    {
    public:
      typedef TmplParticle Particle;
      typedef boost::object_pool<Particle>         ParticlePool;
      typedef boost::shared_ptr<ParticlePool>      ParticlePoolPtr;

    private:
      typedef std::vector<Particle *>              ParticleVector;
      typedef ForwardIterator<ParticleVector>      VectorIterator;
      typedef ForwardConstIterator<ParticleVector> VectorConstIterator;

    public:

      class ParticleIterator : public VectorIterator
      {
      public:
        typedef Particle& value_type;
        ParticleIterator(const VectorIterator &it)
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

      class ParticleConstIterator : public VectorConstIterator
      {
      public:
        typedef const Particle& value_type;
        ParticleConstIterator(const VectorConstIterator &it)
         : VectorConstIterator(it)
        {
        }

        ParticleConstIterator(const VectorIterator &it)
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

      ParticleCollection();

      ParticleCollection(ParticlePoolPtr particlePoolPtr);

      ParticleCollection(const ParticleCollection &p);

      ParticleCollection &operator=(const ParticleCollection &p);
      
      virtual ~ParticleCollection();

      int getNumParticles() const;

      BoundingBox getParticleBBox() const;

      ParticleIterator getParticleIterator()
      {
        return ParticleIterator(VectorIterator(m_particleVector));
      }

      ParticleConstIterator getParticleIterator() const
      {
        return ParticleConstIterator(VectorConstIterator(m_particleVector));
      }

      /**
       * Translates all particle positions by the specified mount.
       *
       * @param vec Translation increment.
       */
      void translateBy(const Vec3 &vec);

      /**
       * Rotates all particles according to the specified rotation.
       *
       * @param rotation Specifies rotation-axis and magnitude,
       *                 ie angle = rotation.norm() radians.
       * @param posn Specifies position of the rotation vector.
       */
      void rotate(const Vec3 &rotation, const Vec3 &posn);

      /**
       * Increments all particle Id's by the specifed amount.
       *
       * @param idIncr Increment.
       */
      void incrementIdBy(typename Particle::Id idIncr);

      /**
       * Adds the specifed particle reference to this collection.
       * @param p Reference to p is inserted.
       */
      void insertRef(Particle &p);

      /**
       * Creates a new particle constructed particle from p.
       * @param p Create particle copy constructed from p.
       * @return returns reference to newly constructed particle.
       */
      Particle &createParticle(const Particle &p);

    protected:
      /**
       * Adds the specifed particle reference to this collection.
       * @param p Reference to p is inserted.
       */
      void noCheckInsertRef(Particle &p);

    private:
      ParticlePoolPtr m_particlePoolPtr;
      ParticleVector  m_particleVector;
    };
  }
}

#include "Geometry/ParticleCollection.hpp"

#endif
