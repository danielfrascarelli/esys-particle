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


#ifndef ESYS_LSMCIRCULARNEIGHBOURTABLE_HPP
#define ESYS_LSMCIRCULARNEIGHBOURTABLE_HPP

#include "Geometry/NeighbourTable.h"
#include <boost/pool/object_pool.hpp>
#include <boost/shared_ptr.hpp>

#include <sstream>
#include <stdexcept>
#include <set>

namespace esys
{
  namespace lsm
  {
    /**
     *
     */
    template <class TmplParticle>
    CircularNeighbourTable<TmplParticle>::CircularNeighbourTable(
      const BoundingBox  &bBox,
      double             gridSpacing,
      const BoolVector   &periodicDimensions,
      double             circBorderWidth
    )
      : Inherited(bBox, gridSpacing),
        m_periodicDimensions(periodicDimensions),
        m_particlePoolPtr(new ParticlePool(4096)),
        m_clonedParticleSet(),
        m_circGridWidth(1),
        m_periodicDimIndex(-1)
    {
      checkPeriodicDimensions();
      if (circBorderWidth <= 0.0) {
        circBorderWidth = this->getGridSpacing();
      }
      setCircularBorderWidth(circBorderWidth, this->getGridSpacing());
    }

    template <class TmplParticle>
    CircularNeighbourTable<TmplParticle>::CircularNeighbourTable(
      const BoundingBox  &bBox,
      double             gridSpacing,
      ParticlePoolPtr    particlePoolPtr,
      const BoolVector   &periodicDimensions,
      double             circBorderWidth
    )
      : Inherited(bBox, gridSpacing),
        m_periodicDimensions(periodicDimensions),
        m_particlePoolPtr(particlePoolPtr),
        m_clonedParticleSet(),
        m_circGridWidth(1),
        m_periodicDimIndex(-1)
    {
      checkPeriodicDimensions();
      if (circBorderWidth <= 0.0) {
        circBorderWidth = this->getGridSpacing();
      }
      setCircularBorderWidth(circBorderWidth, this->getGridSpacing());
    }

    template <class TmplParticle>
    void CircularNeighbourTable<TmplParticle>::checkPeriodicDimensions()
    {
      if (m_periodicDimensions.size() != 3) {
        std::stringstream msg;
        msg
          << "CircularNeighbourTable::CircularNeighbourTable -"
          << " size of periodic dimensions argument ("
          << m_periodicDimensions.size()
          << ") is not equal to 3";
        throw std::runtime_error(msg.str());
      }
      int numPeriodic = 0;
      for (int i = 0; i < 3; i++) {
        if (m_periodicDimensions[i])
        {
          m_periodicDimIndex = i;
          numPeriodic++;
        }
      }

      if (numPeriodic > 1) {
        std::stringstream msg;
        msg
          << "CircularNeighbourTable::CircularNeighbourTable -"
          << " only a single dimension may be periodic.";
        throw std::runtime_error(msg.str());
      }
    }

    template <class TmplParticle>
    CircularNeighbourTable<TmplParticle>::~CircularNeighbourTable()
    {
    }

    template <class TmplParticle>
    void CircularNeighbourTable<TmplParticle>::setCircularBorderWidth(
      double circBorderWidth,
      double gridSpacing
    )
    {
      m_circGridWidth = int(ceil(circBorderWidth/gridSpacing));
    }

    template <class TmplParticle>
    void CircularNeighbourTable<TmplParticle>::setCircularBorderWidth(
      double circBorderWidth
    )
    {
      setCircularBorderWidth(circBorderWidth, this->getGridSpacing());
    }

    template <class TmplParticle>
    void CircularNeighbourTable<TmplParticle>::resize(
      const BoundingBox &bBox,
      double gridSpacing,
      double circBorderWidth
    )
    {
      if (havePeriodicDimensions())
      {
        circBorderWidth =
          min(
            (bBox.getSizes()[m_periodicDimIndex])/2.0,
            circBorderWidth
          );
      }
      setCircularBorderWidth(circBorderWidth, gridSpacing);
      ParticleVector particles = getNonClonedParticles();
      clearClonedParticles();
      this->clearAndRecomputeGrid(bBox, gridSpacing);
      for (
        typename ParticleVector::iterator it = particles.begin();
        it != particles.end();
        it++
      )
      {
        this->insert(*it);
      }
    }

    template <class TmplParticle>
    void CircularNeighbourTable<TmplParticle>::resize(
      const BoundingBox &bBox,
      double gridSpacing
    )
    {
      this->resize(bBox, gridSpacing, gridSpacing);
    }
      
    template <class TmplParticle>
    void CircularNeighbourTable<TmplParticle>::insertClone(
      Particle *pParticle,
      const Vec3 &newPosition
    )
    {
      Particle *pNewParticle = m_particlePoolPtr->construct(*pParticle);
      pNewParticle->moveTo(newPosition);
      Inherited::insert(pNewParticle);
      m_clonedParticleSet.insert(pNewParticle);
    }

    template <class TmplParticle>
    bool CircularNeighbourTable<TmplParticle>::havePeriodicDimensions() const
    {
      return (m_periodicDimIndex >= 0);
    }

    template <class TmplParticle>
    Vec3 CircularNeighbourTable<TmplParticle>::getModdedPosn(
      const Vec3 &posn
    ) const
    {
      if (
        havePeriodicDimensions()
      )
      {
        const int &dimIdx = m_periodicDimIndex;
        if (
          (posn[dimIdx] < this->getBBox().getMinPt()[dimIdx])
          ||
          (posn[dimIdx] > this->getBBox().getMaxPt()[dimIdx])
        )
        {
          Vec3 moddedPosn = posn;
          const double dimSize = this->getBBox().getSizes()[dimIdx];
          moddedPosn[dimIdx] -= this->getBBox().getMinPt()[dimIdx];
          moddedPosn[dimIdx] =
            (moddedPosn[dimIdx] > 0.0)
            ?
            (
              this->getBBox().getMinPt()[dimIdx]
              +
              moddedPosn[dimIdx] - (floor(moddedPosn[dimIdx]/dimSize)*dimSize)
            )
            :
            (
              this->getBBox().getMaxPt()[dimIdx]
              -
              (fabs(moddedPosn[dimIdx]) - (floor(fabs(moddedPosn[dimIdx])/dimSize)*dimSize))
            );

          return moddedPosn;
        }
      }
      return posn;
    }

    template <class TmplParticle>
    void CircularNeighbourTable<TmplParticle>::insert(Particle *pParticle)
    {
      pParticle->moveTo(getModdedPosn(pParticle->getPos()));
      const Vec3L minIdx = this->getVecIndex(pParticle->getPos() - pParticle->getRad());
      const Vec3L maxIdx = this->getVecIndex(pParticle->getPos() + pParticle->getRad());

      this->insertInTable(pParticle, minIdx, maxIdx);
      this->addInserted(pParticle);

      if (havePeriodicDimensions())
      {
        for (int dimIdx = 0; dimIdx < 3; dimIdx++) {
          if (m_periodicDimensions[dimIdx]) {
            if (minIdx[dimIdx] < (this->getMinVecIndex()[dimIdx] + m_circGridWidth)) {
              Vec3 shift = Vec3::ZERO;
              shift[dimIdx] = this->getBBox().getSizes()[dimIdx];
              this->insertClone(pParticle, pParticle->getPos() + shift);
            }
            if (maxIdx[dimIdx] > (this->getMaxVecIndex()[dimIdx] - m_circGridWidth)) {
              Vec3 shift = Vec3::ZERO;
              shift[dimIdx] = this->getBBox().getSizes()[dimIdx];
              this->insertClone(pParticle, pParticle->getPos() - shift);
            }
          }
        }
      }
    }

    template <class TmplParticle>
    void CircularNeighbourTable<TmplParticle>::insert(Particle &particle)
    {
      this->insert(&particle);
    }

    template <class TmplParticle>
    size_t CircularNeighbourTable<TmplParticle>::getNumClonedParticles() const
    {
      return m_clonedParticleSet.size();
    }

    template <class TmplParticle>
    size_t CircularNeighbourTable<TmplParticle>::getNumParticles() const
    {
      return this->size() - getNumClonedParticles();
    }

    template <class TmplParticle>
    const typename CircularNeighbourTable<TmplParticle>::BoolVector &
    CircularNeighbourTable<TmplParticle>::getPeriodicDimensions() const
    {
      return m_periodicDimensions;
    }

    template <class TmplParticle>
    bool CircularNeighbourTable<TmplParticle>::isClone(
      Particle *p
    ) const
    {
      return (m_clonedParticleSet.find(p) != m_clonedParticleSet.end());
    }
      
    template <class TmplParticle>
    typename CircularNeighbourTable<TmplParticle>::ParticleVector
    CircularNeighbourTable<TmplParticle>::getNonClonedParticles()
    {
      ParticleVector all = this->getInsertedParticles();
      ParticleVector nonCloned;
      nonCloned.reserve(all.size()/2);
      for (
        typename ParticleVector::iterator it = all.begin();
        it != all.end();
        it++
      )
      {
        if (!isClone(*it))
        {
          nonCloned.push_back(*it);
        }
      }
      return nonCloned;
    }

    template <class TmplParticle>
    void CircularNeighbourTable<TmplParticle>::clearClonedParticles()
    {
      for (
        typename ParticleSet::iterator it = m_clonedParticleSet.begin();
        it != m_clonedParticleSet.end();
        it++
      )
      {
        m_particlePoolPtr->destroy(*it);
      }
      m_clonedParticleSet.clear();
    }
  }
}

#endif
