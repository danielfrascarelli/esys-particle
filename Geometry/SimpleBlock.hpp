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


#include "Geometry/SimpleBlock.h"
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <boost/limits.hpp>

namespace esys
{
  namespace lsm
  {
    template <typename TmplParticle>
    SimpleBlockGenerator<TmplParticle>::SimpleBlockGenerator(
      unsigned int numX,
      unsigned int numY,
      unsigned int numZ,
      double radius
    )
      : m_radius(radius),
        m_dimCounts(numX, numY, numZ)
    {
    }

    template <typename TmplParticle>
    Vec3 SimpleBlockGenerator<TmplParticle>::getPos(const Vec3L &idx)
    {
      return
        Vec3(
          idx[0]*2.0*getRadius() + getRadius(),
          idx[1]*2.0*getRadius() + getRadius(),
          idx[2]*2.0*getRadius() + ((m_dimCounts[2] > 1) ? getRadius() : 0.0)
        );
    }

    template <typename TmplParticle>
    int SimpleBlockGenerator<TmplParticle>::getId(const Vec3L &idx)
    {
      return
        idx[0]
        +
        idx[1]*m_dimCounts[0]
        +
        idx[2]*m_dimCounts[0]*m_dimCounts[1];
    }

    template <typename TmplParticle>
    template <typename TmplParticleCollection>
    void
    SimpleBlockGenerator<TmplParticle>::createParticles(
      TmplParticleCollection &particleCollection
    )
    {
      Vec3L idx(0, 0, 0);
      for (idx[2]=0; idx[2] < m_dimCounts[2]; (idx[2])++)
      {
        for (idx[1]=0; idx[1] < m_dimCounts[1]; (idx[1])++)
        {
          for (idx[0]=0; idx[0] < m_dimCounts[0]; (idx[0])++)
          {
            particleCollection.createParticle(
              TmplParticle(
                getPos(idx),
                getRadius(),
                getId(idx),
                0
              )
            );
          }
        }
      }
    }

    template <typename TmplParticle>
    SimpleBlockGenerator<TmplParticle>::~SimpleBlockGenerator()
    {
    }

    template <typename TmplParticle>
    double SimpleBlockGenerator<TmplParticle>::getRadius() const
    {
      return m_radius;
    }







    template <typename TmplParticle>
    SimpleBlock<TmplParticle>::SimpleBlock(
      unsigned int numX,
      unsigned int numY,
      unsigned int numZ,
      double radius
    )
      : ParticleCollection<TmplParticle>(),
        m_generator(numX, numY, numZ, radius)
    {
      createParticles();
    }

    template <typename TmplParticle>
    SimpleBlock<TmplParticle>::~SimpleBlock()
    {
    }
    
    template <typename TmplParticle>
    void SimpleBlock<TmplParticle>::createParticles()
    {
      m_generator.createParticles(*this);
    }

    template <typename TmplParticle>
    double SimpleBlock<TmplParticle>::getRadius() const
    {
      return m_generator.getRadius();
    }
  }
}
