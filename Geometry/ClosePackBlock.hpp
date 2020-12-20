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


#ifndef ESYS_LSMCLOSEPACKBLOCK_HPP
#define ESYS_LSMCLOSEPACKBLOCK_HPP

#include "Geometry/ClosePackBlock.h"
#include "Geometry/ClosePackIterator.h"

namespace esys
{
  namespace lsm
  {
    template <typename TmplClosePackIterator,typename TmplParticle>
    ClosePackBlockGenerator<TmplClosePackIterator,TmplParticle>::ClosePackBlockGenerator(
      unsigned int numX,
      unsigned int numY,
      unsigned int numZ,
      double radius,
      ClosePackOrientation orientation
    )
      : m_radius(radius),
        m_dimCounts(numX, numY, numZ),
        m_orientation(orientation)
    {
    }

    template <typename TmplClosePackIterator,typename TmplParticle>
    template <typename TmplParticleCollection>
    void
    ClosePackBlockGenerator<TmplClosePackIterator,TmplParticle>::createParticles(
      TmplParticleCollection &particleCollection
    )
    {
      int id = 0;
      CentrePointIterator it =
        CentrePointIterator(
          m_dimCounts[0],
          m_dimCounts[1],
          m_dimCounts[2],
          getRadius(),
          ((m_orientation == DEFAULT_ORIENT) && (m_dimCounts[2] <= 1)) ? XYZ : m_orientation
        );

      while (it.hasNext())
      {
        particleCollection.createParticle(
          TmplParticle(
            it.next(),
            getRadius(),
            id,
            0
          )
        );
        id++;
      }
    }

    template <typename TmplClosePackIterator,typename TmplParticle>
    ClosePackBlockGenerator<TmplClosePackIterator,TmplParticle>::~ClosePackBlockGenerator()
    {
    }

    template <typename TmplClosePackIterator,typename TmplParticle>
    double ClosePackBlockGenerator<TmplClosePackIterator,TmplParticle>::getRadius() const
    {
      return m_radius;
    }







    template <typename TmplClosePackIterator,typename TmplParticle>
    ClosePackBlock<TmplClosePackIterator,TmplParticle>::ClosePackBlock(
      unsigned int numX,
      unsigned int numY,
      unsigned int numZ,
      double radius,
      ClosePackOrientation orientation
    )
      : ParticleCollection<TmplParticle>(),
        m_generator(numX, numY, numZ, radius)
    {
      createParticles();
    }

    template <typename TmplClosePackIterator,typename TmplParticle>
    ClosePackBlock<TmplClosePackIterator,TmplParticle>::~ClosePackBlock()
    {
    }
    
    template <typename TmplClosePackIterator,typename TmplParticle>
    void ClosePackBlock<TmplClosePackIterator,TmplParticle>::createParticles()
    {
      m_generator.createParticles(*this);
    }

    template <typename TmplClosePackIterator,typename TmplParticle>
    double ClosePackBlock<TmplClosePackIterator,TmplParticle>::getRadius() const
    {
      return m_generator.getRadius();
    }
  }
}

#endif
