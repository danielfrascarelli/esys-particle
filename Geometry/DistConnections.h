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


#ifndef ESYS_LSMDISTCONNECTIONS_H
#define ESYS_LSMDISTCONNECTIONS_H

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
    template <typename TmplParticle, typename TmplConnection>
    class DistConnections
    {
    public:
      typedef TmplParticle   Particle;
      typedef TmplConnection Connection;
      typedef int            Tag;

      class Cmp
      {
      public:
        bool operator()(const Connection &c1, const Connection &c2) const
        {
          return
            (
              (c1.first() < c2.first())
              ||
              (
                (c1.first() == c2.first())
                &&
                (
                  (c1.second() < c2.second())
                  ||
                  (
                    (c1.second() == c2.second())
                    &&
                    (c1.getTag() < c2.getTag())
                  )
                )
              )
            );
        }
        bool operator()(const Connection *c1, const Connection *c2) const
        {
          return (*this)(*c1, *c2);
        }
      };
    public:
      typedef std::set<Connection *,Cmp>             ConnectionSet;
      typedef std::vector<Particle *>                ParticleVector;
      typedef CircularNeighbourTable<Particle>       NTable;
      typedef typename NTable::ParticleIterator      ParticleIterator;
      typedef typename NTable::ParticleConstIterator ParticleConstIterator;
      
    public:
      typedef typename NTable::BoolVector BoolVector;

      DistConnections(
        double maxDist,
        Tag defaultTag = 0,
        const BoundingBox &bBox = BoundingBox(Vec3(-10,-10,-10), Vec3(10,10,10)),
        const BoolVector  &circDimensions = BoolVector(3, false)
      );

      ~DistConnections();

      int getNumParticles() const;

      int getNumConnections() const;

      double getMinRadius() const;

      double getMaxRadius() const;

      ParticleConstIterator getParticleIterator() const;

      BoundingBox getParticleBBox() const;

      template<typename TmplParticleIterator>
      void create(TmplParticleIterator it);

      template<typename TmplParticleIterator>
      void create(TmplParticleIterator it, Tag tag);

      Tag getDefaultTag() const;
      void setDefaultTag(Tag defaultTag);

      typedef ForwardConstIterator<ConnectionSet> ConnectionConstIterator;

      class ConstIterator : public ConnectionConstIterator
      {
      public:
        typedef const Connection& value_type;
        typedef const Connection& reference;
        ConstIterator(const ConnectionSet &set)
         : ConnectionConstIterator(set)
        {
        }

        value_type next()
        {
          return *(ConnectionConstIterator::next());
        }

        value_type current() const
        {
          return *(ConnectionConstIterator::current());
        }
      };
      typedef ConstIterator Iterator;

      Iterator getIterator() const
      {
        return Iterator(m_connectionSet);
      }

    protected:
      void insert(Particle &p);

      void createConnection(const Particle &p1, const Particle &p2, Tag tag);

    private:
      typedef boost::shared_ptr<NTable>         NTablePtr;
      typedef boost::object_pool<Connection>    ConnectionPool;
      typedef boost::shared_ptr<ConnectionPool> ConnectionPoolPtr;

      ConnectionPoolPtr m_connectionPoolPtr;
      ConnectionSet     m_connectionSet;
      NTablePtr         m_nTablePtr;
      double            m_minRadius;
      double            m_maxRadius;
      double            m_maxDist;
      Vec3              m_minPt;
      Vec3              m_maxPt;
      Tag               m_defaultTag;
    };
  }
}

#include "Geometry/DistConnections.hpp"

#endif
