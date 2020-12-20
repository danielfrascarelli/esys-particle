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


#ifndef ESYS_LSMCLOSEPACKBLOCK_H
#define ESYS_LSMCLOSEPACKBLOCK_H

#include "Geometry/ClosePackIterator.h"
#include "Geometry/ParticleCollection.h"

namespace esys
{
  namespace lsm
  {
    /**
     *
     */
    template <typename TmplClosePackIterator,typename TmplParticle>
    class ClosePackBlockGenerator
    {
    public:
      typedef TmplClosePackIterator CentrePointIterator;
      
      ClosePackBlockGenerator(
          unsigned int numX,
          unsigned int numY,
          unsigned int numZ,
          double radius = 0.5,
          ClosePackOrientation orientation = DEFAULT_ORIENT
      );

      virtual ~ClosePackBlockGenerator();

      double getRadius() const;

      template <typename TmplParticleCollection>
      void createParticles(TmplParticleCollection &particleCollection);

    protected:

    private:
      double               m_radius;
      Vec3L                m_dimCounts;
      ClosePackOrientation m_orientation;
    };

    /**
     *
     */
    template <typename TmplClosePackIterator, typename TmplParticle>
    class ClosePackBlock : public ParticleCollection<TmplParticle>
    {
    public:
      typedef typename ParticleCollection<TmplParticle>::Particle Particle;
      typedef TmplClosePackIterator ClosePackIterator;
      typedef ClosePackBlockGenerator<TmplClosePackIterator, Particle> BlockGenerator;
      ClosePackBlock(
          unsigned int numX,
          unsigned int numY,
          unsigned int numZ,
          double radius = 0.5,
          ClosePackOrientation orientation = DEFAULT_ORIENT
      );

      virtual ~ClosePackBlock();

      double getRadius() const;

    protected:
      void createParticles();

    private:
      BlockGenerator m_generator;
    };
  }
}

#include "Geometry/ClosePackBlock.hpp"

#endif
