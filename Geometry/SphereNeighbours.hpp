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
#include <functional>
#include <limits>

namespace esys
{
  namespace lsm
  {
    template <typename TmplSphere, typename TmplIdPair>
    SphereNeighbours<TmplSphere,TmplIdPair>::SphereNeighbours(
      double maxDist,
      const BoundingBox &bBox,
      const BoolVector &circDimensions
    )
      : m_connectionPoolPtr(new IdPairPool(4096)),
        m_connectionSet(),
        m_nTablePtr(),
        m_minRadius(std::numeric_limits<double>::max()),
        m_maxRadius(-std::numeric_limits<double>::max()),
        m_maxDist(maxDist),
        m_minPt(bBox.getMinPt()),
        m_maxPt(bBox.getMaxPt())
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

    template <typename TmplSphere, typename TmplIdPair>
    SphereNeighbours<TmplSphere,TmplIdPair>::~SphereNeighbours()
    {
    }

    template <typename TmplSphere, typename TmplIdPair>
    int SphereNeighbours<TmplSphere,TmplIdPair>::getNumSpheres() const
    {
      return m_nTablePtr->size();
    }

    template <typename TmplSphere, typename TmplIdPair>
    int SphereNeighbours<TmplSphere,TmplIdPair>::getNumIdPairs() const
    {
      return m_connectionSet.size();
    }

    template <typename TmplSphere, typename TmplIdPair>
    double SphereNeighbours<TmplSphere,TmplIdPair>::getMinRadius() const
    {
      return m_minRadius;
    }

    template <typename TmplSphere, typename TmplIdPair>
    double SphereNeighbours<TmplSphere,TmplIdPair>::getMaxRadius() const
    {
      return m_maxRadius;
    }

    template <typename TmplSphere, typename TmplIdPair>
    typename SphereNeighbours<TmplSphere,TmplIdPair>::SphereConstIterator
    SphereNeighbours<TmplSphere,TmplIdPair>::getSphereIterator() const
    {
      return m_nTablePtr->getIterator();
    }

    template <typename TmplSphere, typename TmplIdPair>
    const typename SphereNeighbours<TmplSphere,TmplIdPair>::IdPair &
    SphereNeighbours<TmplSphere,TmplIdPair>::createIdPair(
      const Sphere &p1,
      const Sphere &p2
    )
    {
      return 
        **(m_connectionSet.insert(
          m_connectionPoolPtr->construct(p1.getId(), p2.getId())
        ).first);
    }

    template <typename TmplSphere>
    class CmpSphereId
    {
    public:
      bool operator()(const TmplSphere &p1, const TmplSphere &p2) const
      {
        return (p1.getId() < p2.getId());
      }

      bool operator()(const TmplSphere *p1, const TmplSphere *p2) const
      {
        return (p1->getId() < p2->getId());
      }
    };

    template <typename TmplType>
    class Deref
    {
    public:
      typedef const TmplType& result_type;
      typedef const TmplType* argument_type;
  
      result_type
      operator()(argument_type a) const
      { return *a; }
    };
    
    template <typename TmplSphere, typename TmplIdPair>
    template <typename TmplSphereIterator>
    typename SphereNeighbours<TmplSphere,TmplIdPair>::IdPairVector
    SphereNeighbours<TmplSphere,TmplIdPair>::getNeighbours(
      TmplSphereIterator it
    )
    {
      // Put the spheres in the neighbour table
      typedef std::set<Sphere *, CmpSphereId<Sphere> > SphereSet;
      SphereSet pSet;
      while (it.hasNext())
      {
        Sphere &p = it.next();
        insert(p);
        pSet.insert(&p);
      }
      ConstIdPairSet idPairSet;
      // Resize ntable according to min and max sphere radii.
      m_nTablePtr->resize(
        getSphereBBox(),
        4.1*getMinRadius(),
        2.1*getMaxRadius()
      );

      // For each sphere in the iterator it, determine it's
      // neighours within m_maxDist distance.
      for (
        typename SphereSet::const_iterator pIt = pSet.begin();
        pIt != pSet.end();
        pIt++
      )
      {
        // get the vector of spheres which are contained in nearby cells.
        typename NTable::ParticleVector nVector =
          m_nTablePtr->getNeighbourVector(
            (*pIt)->getPos(),
            (*pIt)->getRad() + m_maxDist
          );
        // for each of these cell-spheres, determine if they
        // are closer than m_maxDist.
        for (
          typename NTable::ParticleVector::const_iterator nIt = nVector.begin();
          nIt != nVector.end();
          nIt++
        )
        {
          Sphere *p1 = (*pIt);
          Sphere *p2 = (*nIt);

          if (
            (
              (p1->getId() < p2->getId())
              &&
              (pSet.find(p1) != pSet.end())
              &&
              (pSet.find(p2) != pSet.end())
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
            const double radiusSumPlusTol =
              m_maxDist + p1->getRad() + p2->getRad();
            const double radiusSumPlusTolSqrd =
              radiusSumPlusTol*radiusSumPlusTol;
  
            if (
              (p1->getPos() - p2->getPos()).norm2()
              <=
              (radiusSumPlusTolSqrd)
            )
            {
              idPairSet.insert(&createIdPair(*p1, *p2));
            }
          }
        }
      }
      IdPairVector idPairVector;
      idPairVector.reserve(idPairSet.size());
      std::transform(
        idPairSet.begin(),
        idPairSet.end(),
        std::back_insert_iterator<IdPairVector>(idPairVector),
        Deref<IdPair>()
      );
      return idPairVector;
    }

    template <typename TmplSphere, typename TmplIdPair>
    void
    SphereNeighbours<TmplSphere,TmplIdPair>::insert(Sphere &p)
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

    template <typename TmplSphere, typename TmplIdPair>
    BoundingBox
    SphereNeighbours<TmplSphere,TmplIdPair>::getSphereBBox() const
    {
      return BoundingBox(m_minPt, m_maxPt);
    }
  }
}
