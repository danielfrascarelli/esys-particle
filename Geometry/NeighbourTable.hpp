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


#ifndef ESYS_LSMNEIGHBOURTABLE_HPP
#define ESYS_LSMNEIGHBOURTABLE_HPP

namespace esys
{
  namespace lsm
  {
    template <class TmplParticle>
    NeighbourTable<TmplParticle>::NeighbourTable(
      const BoundingBox &bBox,
      double gridSpacing
    )
      : m_dimensions(),
        m_minIndex(),
        m_maxIndex(Vec3L(-1, -1, -1)),
        m_gridSpacing(gridSpacing),
        m_bBox(bBox),
        m_insertedParticles(),
        m_tablePtr()
    {
      resize(bBox, gridSpacing);
    }

    template <class TmplParticle>
    NeighbourTable<TmplParticle>::NeighbourTable(
      const NeighbourTable &nTable
    )
      : m_dimensions(nTable.m_dimensions),
        m_minIndex(nTable.m_minIndex),
        m_maxIndex(nTable.m_maxIndex),
        m_gridSpacing(nTable.m_gridSpacing),
        m_bBox(nTable.m_bBox),
        m_insertedParticles(nTable.m_insertedParticles),
        m_tablePtr()
    {
      m_tablePtr =
        ParticleVectorArrayPtr(
          new ParticleVector[nTable.getNumCells()]
        );
      for (int i = 0; i < nTable.getNumCells(); i++)
      {
        m_tablePtr[i] = nTable.m_tablePtr[i];
      }
    }

    template <class TmplParticle>
    NeighbourTable<TmplParticle>::~NeighbourTable()
    {
    }
      
    template <class TmplParticle>
    void NeighbourTable<TmplParticle>::clear()
    {
      for (int i = getMinVecIndex().X(); i <= getMaxVecIndex().X(); i++) {
        for (int j = getMinVecIndex().Y(); j <= getMaxVecIndex().Y(); j++) {
          for (int k = getMinVecIndex().Z(); k <= getMaxVecIndex().Z(); k++) {
            m_tablePtr[getScalarIndex(i, j, k)].clear();
          }
        }
      }
      m_insertedParticles.clear();
    }

    template <class TmplParticle>
    double NeighbourTable<TmplParticle>::getGridSpacing() const
    {
      return m_gridSpacing;
    }

    template <class TmplParticle>
    void NeighbourTable<TmplParticle>::resize(
      const BoundingBox &bBox,
      double gridSpacing
    )
    {
      ParticleVector particles = getInsertedParticles();
      clearAndRecomputeGrid(bBox, gridSpacing);
      for (
        typename ParticleVector::iterator it = particles.begin();
        it != particles.end();
        it++
      )
      {
        insert(*it);
      }
    }

    template <class TmplParticle>
    const Vec3L &
    NeighbourTable<TmplParticle>::getDimensions() const
    {
      return m_dimensions;
    }

    template <class TmplParticle>
    const BoundingBox &
    NeighbourTable<TmplParticle>::getBBox() const
    {
      return m_bBox;
    }

    template <class TmplParticle>
    const Vec3 &
    NeighbourTable<TmplParticle>::getMinPt() const
    {
      return getBBox().getMinPt();
    }

    template <class TmplParticle>
    size_t NeighbourTable<TmplParticle>::size() const
    {
      return m_insertedParticles.size();
    }

    template <class TmplParticle>
    int NeighbourTable<TmplParticle>::getScalarIndex(
      int xIdx,
      int yIdx,
      int zIdx
    ) const
    {
      return
        xIdx*m_dimensions.Z()*m_dimensions.Y()
        +
        yIdx*m_dimensions.Z()
        +
        zIdx;
    }

    template <class TmplParticle>
    int NeighbourTable<TmplParticle>::getScalarIndex(const Vec3L &index) const
    {
      return getScalarIndex(index.X(), index.Y(), index.Z());
    }

    template <class TmplParticle>
    int
    NeighbourTable<TmplParticle>::getScalarIndex(const Vec3 &pt) const
    {
      return getScalarIndex(getVecIndex(pt));
    }

    template <class TmplParticle>
    const Vec3L &
    NeighbourTable<TmplParticle>::getMinVecIndex() const
    {
      return m_minIndex;
    }

    template <class TmplParticle>
    const Vec3L &
    NeighbourTable<TmplParticle>::getMaxVecIndex() const
    {
      return m_maxIndex;
    }

    template <class TmplParticle>
    Vec3L
    NeighbourTable<TmplParticle>::getVecIndex(const Vec3 &pt) const
    {
      const Vec3 relPos = Vec3((pt - getMinPt())/m_gridSpacing);
      const Vec3L index = Vec3L(int(floor(relPos.X())), int(floor(relPos.Y())), int(floor(relPos.Z())));
      return getMinVecIndex().max(getMaxVecIndex().min(index));
    }

    template <class TmplParticle>
    typename NeighbourTable<TmplParticle>::ParticleVector
    NeighbourTable<TmplParticle>::getNeighbourVector(
      const Vec3 &pt,
      double radius
    ) const
    {
      ParticleVector neighbours;
      neighbours.reserve(128);
      const Vec3L min = getVecIndex(pt - radius);
      const Vec3L max = getVecIndex(pt + radius);
      for (int i = min.X(); i <= max.X(); i++) {
        for (int j = min.Y(); j <= max.Y(); j++) {
          for (int k = min.Z(); k <= max.Z(); k++) {
            neighbours.insert(
              neighbours.end(),
              m_tablePtr[getScalarIndex(i, j, k)].begin(),
              m_tablePtr[getScalarIndex(i, j, k)].end()
            );
          }
        }
      }
      return neighbours;
    }

    template <class TmplParticle>
    typename NeighbourTable<TmplParticle>::ParticleVector
    NeighbourTable<TmplParticle>::getUniqueNeighbourVector(
      const Vec3 &pt,
      double radius
    ) const
    {
      ParticleVector neighbours = getNeighbourVector(pt, radius);
      std::sort(neighbours.begin(), neighbours.end());
      typename ParticleVector::iterator uniqueEnd =
        std::unique(neighbours.begin(), neighbours.end());
      neighbours.erase(
        uniqueEnd,
        neighbours.end()
      );

      return neighbours;
    }

    template <class TmplParticle>
    typename NeighbourTable<TmplParticle>::ParticleVector
    NeighbourTable<TmplParticle>::getNeighbourVector(
      const Vec3 &pt
    ) const
    {
      return m_tablePtr[getScalarIndex(pt)];
    }

    template <class TmplParticle>
    void NeighbourTable<TmplParticle>::insert(Particle *pParticle)
    {
      const Vec3L minIdx = getVecIndex(pParticle->getPos() - pParticle->getRad());
      const Vec3L maxIdx = getVecIndex(pParticle->getPos() + pParticle->getRad());
      insertInTable(pParticle, minIdx, maxIdx);
      addInserted(pParticle);
    }

    template <class TmplParticle>
    void NeighbourTable<TmplParticle>::insert(Particle &particle)
    {
      insert(&particle);
    }

    template <class TmplParticle>
    typename NeighbourTable<TmplParticle>::ParticleIterator
    NeighbourTable<TmplParticle>::getParticleIterator()
    {
      return ParticleIterator(m_insertedParticles);
    }

    template <class TmplParticle>
    typename NeighbourTable<TmplParticle>::ParticleConstIterator
    NeighbourTable<TmplParticle>::getParticleIterator() const
    {
      return ParticleConstIterator(m_insertedParticles);
    }

    template <class TmplParticle>
    void NeighbourTable<TmplParticle>::insertInTable(
      Particle *pParticle,
      const Vec3L &minIdx,
      const Vec3L &maxIdx
    )
    {
      for (int i = minIdx.X(); i <= maxIdx.X(); i++) {
        for (int j = minIdx.Y(); j <= maxIdx.Y(); j++) {
          for (int k = minIdx.Z(); k <= maxIdx.Z(); k++) {
            m_tablePtr[getScalarIndex(i, j, k)].push_back(pParticle);
          }
        }
      }
    }

    template <class TmplParticle>
    void NeighbourTable<TmplParticle>::addInserted(Particle *pParticle)
    {
      m_insertedParticles.push_back(pParticle);
    }

    template <class TmplParticle>
    int NeighbourTable<TmplParticle>::getNumCells() const
    {
      return getDimensions()[0]*getDimensions()[1]*getDimensions()[2];
    }

    template <class TmplParticle>
    typename NeighbourTable<TmplParticle>::ParticleVector
    NeighbourTable<TmplParticle>::getInsertedParticles() const
    {
      return m_insertedParticles;
    }

    template <class TmplParticle>
    void NeighbourTable<TmplParticle>::clearAndRecomputeGrid(
      const BoundingBox &bBox,
      double gridSpacing
    )
    {
      clear();
      m_bBox = bBox;
      m_gridSpacing = gridSpacing;

      const Vec3 dims  = m_bBox.getSizes()/gridSpacing;
      m_dimensions =
        Vec3L(
          int(floor(dims[0])),
          int(floor(dims[1])),
          int(floor(dims[2]))
        );
      m_dimensions = m_dimensions.max(Vec3L(1, 1, 1));

      m_tablePtr =
        ParticleVectorArrayPtr(
          new ParticleVector[
            m_dimensions.X()*m_dimensions.Y()*m_dimensions.Z()
          ]
        );
      m_minIndex = Vec3L(0, 0, 0);
      m_maxIndex = (m_dimensions - 1);
    }
  }
}

#endif
