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


#include "Foundation/StringUtil.h"
#include <limits>

namespace esys
{
  namespace lsm
  {
    template <typename TmplParticle, typename TmplConnection>
    DistConnections<TmplParticle,TmplConnection>::DistConnections(
      double maxDist,
      Tag defaultTag,
      const BoundingBox &bBox,
      const BoolVector &circDimensions
    )
      : m_connectionPoolPtr(new ConnectionPool(4096)),
        m_connectionSet(),
        m_nTablePtr(),
        m_minRadius(std::numeric_limits<double>::max()),
        m_maxRadius(-std::numeric_limits<double>::max()),
        m_maxDist(maxDist),
        m_minPt(bBox.getMinPt()),
        m_maxPt(bBox.getMaxPt()),
        m_defaultTag(defaultTag)
    {
      const double gridSize =
        (
          max(
            bBox.getSizes()[0],
            max(
              bBox.getSizes()[1],
              bBox.getSizes()[2]
            )
          )
        )/5.0;
      m_nTablePtr =
        NTablePtr(
          new NTable(
            bBox,
            gridSize,
            circDimensions,
            2*gridSize
          )
      );
    }

    template <typename TmplParticle, typename TmplConnection>
    DistConnections<TmplParticle,TmplConnection>::~DistConnections()
    {
    }

    template <typename TmplParticle, typename TmplConnection>
    int DistConnections<TmplParticle,TmplConnection>::getNumParticles() const
    {
      return m_nTablePtr->size();
    }

    template <typename TmplParticle, typename TmplConnection>
    int DistConnections<TmplParticle,TmplConnection>::getNumConnections() const
    {
      return m_connectionSet.size();
    }

    template <typename TmplParticle, typename TmplConnection>
    double DistConnections<TmplParticle,TmplConnection>::getMinRadius() const
    {
      return m_minRadius;
    }

    template <typename TmplParticle, typename TmplConnection>
    double DistConnections<TmplParticle,TmplConnection>::getMaxRadius() const
    {
      return m_maxRadius;
    }

    template <typename TmplParticle, typename TmplConnection>
    typename DistConnections<TmplParticle,TmplConnection>::ParticleConstIterator
    DistConnections<TmplParticle,TmplConnection>::getParticleIterator() const
    {
      return m_nTablePtr->getIterator();
    }

    template <typename TmplParticle, typename TmplConnection>
    void
    DistConnections<TmplParticle,TmplConnection>::createConnection(
      const Particle &p1,
      const Particle &p2,
      Tag tag
    )
    {
      m_connectionSet.insert(
        m_connectionPoolPtr->construct(p1.getId(), p2.getId(), tag)
      );
    }
    template <typename TmplParticle>
    class CmpParticleId
    {
    public:
      bool operator()(const TmplParticle &p1, const TmplParticle &p2) const
      {
        return (p1.getId() < p2.getId());
      }

      bool operator()(const TmplParticle *p1, const TmplParticle *p2) const
      {
        return (p1->getId() < p2->getId());
      }
    };

    template <typename TmplParticle, typename TmplConnection>
    template <typename TmplParticleIterator>
    void
    DistConnections<TmplParticle,TmplConnection>::create(
      TmplParticleIterator it,
      Tag tag
    )
    {
      typedef std::set<Particle *, CmpParticleId<Particle> > ParticleSet;
      ParticleSet pSet;
      while (it.hasNext())
      {
        Particle &p = it.next();
        insert(p);
        pSet.insert(&p);
      }
      m_nTablePtr->resize(getParticleBBox(), 4.1*getMinRadius(), 2.1*getMaxRadius());

      for (
        typename ParticleSet::const_iterator pIt = pSet.begin();
        pIt != pSet.end();
        pIt++
      )
      {
        typename NTable::ParticleVector nVector =
          m_nTablePtr->getNeighbourVector(
            (*pIt)->getPos(),
            (*pIt)->getRad() + m_maxDist
          );
        for (
          typename NTable::ParticleVector::const_iterator nIt = nVector.begin();
          nIt != nVector.end();
          nIt++
        )
        {
          Particle *p1 = (*pIt);
          Particle *p2 = (*nIt);

          if (
            (
              (pSet.find(p1) != pSet.end())
              &&
              (pSet.find(p2) != pSet.end())
              &&
              (p1->getId() < p2->getId())
            )
            ||
            (
              ((pSet.find(p1)==pSet.end()) && (pSet.find(p2)!= pSet.end()))
              ||
              ((pSet.find(p1)!=pSet.end()) && (pSet.find(p2)== pSet.end()))
            )
          )
          {
            p1 =
              ((*pIt)->getId() < (*nIt)->getId())
              ?
              (*pIt)
              :
              (*nIt);
            p2 =
              ((*pIt)->getId() < (*nIt)->getId())
              ?
              (*nIt)
              :
              (*pIt);
            const double radiusSum =
              m_maxDist + p1->getRad() + p2->getRad();
            const double radiusSumSqrd = radiusSum*radiusSum;
  
            if (
              (p1->getPos() - p2->getPos()).norm2()
              <=
              (radiusSumSqrd)
            )
            {
#if 0
              console.Debug()
                << "creating connection: \n"
                << StringUtil::toString(*p1)
                << "->"
                << StringUtil::toString(*p2) << "\n";
#endif
              createConnection(*p1, *p2, tag);
            }
          }
        }
      }
    }

    template <typename TmplParticle, typename TmplConnection>
    template <typename TmplParticleIterator>
    void
    DistConnections<TmplParticle,TmplConnection>::create(
      TmplParticleIterator it
    )
    {
      create(it, getDefaultTag());
    }
    
    template <typename TmplParticle, typename TmplConnection>
    void
    DistConnections<TmplParticle,TmplConnection>::insert(Particle &p)
    {
      if (p.getRad() < m_minRadius)
      {
        m_minRadius = p.getRad();
      }
      if (p.getRad() > m_maxRadius)
      {
        m_maxRadius = p.getRad();
      }

      m_nTablePtr->insert(p);

      for (int i = 0; i < 3; i++)
      {
        if (!(m_nTablePtr->getPeriodicDimensions()[i]))
        {
          if (p.getPos()[i]-p.getRad() < m_minPt[i])
          {
            m_minPt[i] = p.getPos()[i]-p.getRad();
          }
          if (p.getPos()[i]+p.getRad() > m_maxPt[i])
          {
            m_maxPt[i] = p.getPos()[i]+p.getRad();
          }
        }
      }
    }

    template <typename TmplParticle, typename TmplConnection>
    typename DistConnections<TmplParticle,TmplConnection>::Tag
    DistConnections<TmplParticle,TmplConnection>::getDefaultTag() const
    {
      return m_defaultTag;
    }

    template <typename TmplParticle, typename TmplConnection>
    void
    DistConnections<TmplParticle,TmplConnection>::setDefaultTag(Tag defaultTag)
    {
      m_defaultTag = defaultTag;
    }

    template <typename TmplParticle, typename TmplConnection>
    BoundingBox
    DistConnections<TmplParticle,TmplConnection>::getParticleBBox() const
    {
      return BoundingBox(m_minPt, m_maxPt);
    }

  }
}
