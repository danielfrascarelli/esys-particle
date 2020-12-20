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


#ifndef ESYS_LSMSIMPLEBLOCK_H
#define ESYS_LSMSIMPLEBLOCK_H

#include "Foundation/BoundingBox.h"
#include "Geometry/ParticleCollection.h"
#include "Geometry/SimpleParticle.h"
#include "Geometry/NeighbourTable.h"
#include "Geometry/BasicInteraction.h"
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
    template <typename TmplParticle = SimpleParticle>
    class SimpleBlockGenerator
    {
    public:
      SimpleBlockGenerator(
          unsigned int numX,
          unsigned int numY,
          unsigned int numZ,
          double radius = 0.5
      );

      virtual ~SimpleBlockGenerator();

      double getRadius() const;
    
      template <typename TmplParticleCollection>
      void createParticles(TmplParticleCollection &particleCollection);

    protected:
      Vec3 getPos(const Vec3L &idx);

      int getId(const Vec3L &idx);

    private:
      double       m_radius;
      Vec3L        m_dimCounts;
    };

    /**
     *
     */
    template <typename TmplParticle = SimpleParticle>
    class SimpleBlock : public ParticleCollection<TmplParticle>
    {
    public:
      typedef typename ParticleCollection<TmplParticle>::Particle Particle;
      SimpleBlock(
          unsigned int numX,
          unsigned int numY,
          unsigned int numZ,
          double radius = 0.5
      );

      virtual ~SimpleBlock();

      double getRadius() const;

    protected:
      void createParticles();

    private:
      SimpleBlockGenerator<Particle> m_generator;
    };
  }
}

#include "Geometry/SimpleBlock.hpp"

#endif
