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


#ifndef ESYS_LSMSPHERENEIGHBOURS_H
#define ESYS_LSMSPHERENEIGHBOURS_H

#include "Geometry/CircularNeighbourTable.h"
#include "Geometry/BasicInteraction.h"

#include <boost/shared_ptr.hpp>
#include <boost/pool/object_pool.hpp>

#include <set>
#include <vector>
#include <float.h>

namespace esys
{
  namespace lsm
  {
    /**
     *
     */
    template <typename TmplSphere, typename TmplIdPairVector>
    class SphereNeighbours
    {
    public:
      typedef int                               Id;
      typedef TmplSphere                        Sphere;
      typedef TmplIdPairVector                  IdPairVector;
      typedef typename IdPairVector::value_type IdPair;

      class Cmp
      {
      public:
        bool operator()(const IdPair &c1, const IdPair &c2) const
        {
          return
            (
              (c1.first < c2.first)
              ||
              (
                (c1.first == c2.first)
                &&
                (
                  (c1.second < c2.second)
                )
              )
            );
        }
        bool operator()(const IdPair *c1, const IdPair *c2) const
        {
          return (*this)(*c1, *c2);
        }
      };
    public:
      typedef std::set<IdPair *,Cmp>                 IdPairSet;
      typedef std::set<const IdPair *,Cmp>           ConstIdPairSet;
      typedef std::vector<Sphere *>                  SphereVector;
      typedef CircularNeighbourTable<Sphere>         NTable;
      typedef typename NTable::ParticleIterator      SphereIterator;
      typedef typename NTable::ParticleConstIterator SphereConstIterator;
      
    public:
      typedef typename NTable::BoolVector BoolVector;

      SphereNeighbours(
        double maxDist,
        const BoundingBox &bBox = BoundingBox(Vec3(-10,-10,-10), Vec3(10,10,10)),
        const BoolVector  &circDimensions = BoolVector(3, false)
      );

      ~SphereNeighbours();

      int getNumSpheres() const;

      int getNumIdPairs() const;

      double getMinRadius() const;

      double getMaxRadius() const;

      SphereConstIterator getSphereIterator() const;

      BoundingBox getSphereBBox() const;

      template<typename TmplSphereIterator>
      IdPairVector getNeighbours(TmplSphereIterator it);

      typedef ForwardConstIterator<IdPairSet> IdPairConstIterator;

      class ConstIterator : public IdPairConstIterator
      {
      public:
        typedef const IdPair& value_type;
        typedef const IdPair& reference;
        ConstIterator(const IdPairSet &set)
         : IdPairConstIterator(set)
        {
        }

        value_type next()
        {
          return *(IdPairConstIterator::next());
        }

        value_type current() const
        {
          return *(IdPairConstIterator::current());
        }
      };
      typedef ConstIterator Iterator;

      Iterator getIterator() const
      {
        return Iterator(m_connectionSet);
      }

    protected:
      void insert(Sphere &p);

      const IdPair &createIdPair(const Sphere &p1, const Sphere &p2);

    private:
      typedef boost::shared_ptr<NTable>     NTablePtr;
      typedef boost::object_pool<IdPair>    IdPairPool;
      typedef boost::shared_ptr<IdPairPool> IdPairPoolPtr;

      IdPairPoolPtr m_connectionPoolPtr;
      IdPairSet     m_connectionSet;
      NTablePtr     m_nTablePtr;
      double        m_minRadius;
      double        m_maxRadius;
      double        m_maxDist;
      Vec3          m_minPt;
      Vec3          m_maxPt;
    };
  }
}

#include "Geometry/SphereNeighbours.hpp"

#endif
