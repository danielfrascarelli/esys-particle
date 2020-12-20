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


#ifndef ESYS_LSMNEIGHBOURTABLE_H
#define ESYS_LSMNEIGHBOURTABLE_H

#include <Foundation/BoundingBox.h>
#include <Foundation/StlIterator.h>
#include <Geometry/Vec3L.h>
#include <vector>
#include <algorithm>
#include <boost/shared_array.hpp>

namespace esys
{
  namespace lsm
  {
    /**
     *
     */
    template <class TmplParticle>
    class NeighbourTable
    {
    public:
      typedef TmplParticle            Particle; 
      typedef std::vector<Particle *> ParticleVector;

      NeighbourTable(const BoundingBox &bBox, double gridSpacing);

      NeighbourTable(const NeighbourTable &nTable);
      
      virtual ~NeighbourTable();
      
      void clear();
      
      double getGridSpacing() const;

      void resize(const BoundingBox &bBox, double gridSpacing);

      const Vec3L &getDimensions() const;

      const BoundingBox &getBBox() const;

      const Vec3 &getMinPt() const;

      /**
       * Return the number of particles inserted into this table.
       */
      size_t size() const;

      int getScalarIndex(int xIdx, int yIdx, int zIdx) const;

      int getScalarIndex(const Vec3L &index) const;

      int getScalarIndex(const Vec3 &pt) const;

      const Vec3L &getMinVecIndex() const;

      const Vec3L &getMaxVecIndex() const;

      Vec3L getVecIndex(const Vec3 &pt) const;

      ParticleVector getNeighbourVector(const Vec3 &pt, double radius) const;

      ParticleVector getUniqueNeighbourVector(const Vec3 &pt, double radius) const;

      ParticleVector getNeighbourVector(const Vec3 &pt) const;

      void insert(Particle *pParticle);

      void insert(Particle &particle);

      typedef ForwardIterator<ParticleVector> ParticleIterator;
      typedef ForwardConstIterator<ParticleVector> ParticleConstIterator;

      ParticleIterator getParticleIterator();

      ParticleConstIterator getParticleIterator() const;

    protected:

      void insertInTable(Particle *pParticle, const Vec3L &minIdx, const Vec3L &maxIdx);

      void addInserted(Particle *pParticle);

      int getNumCells() const;

      ParticleVector getInsertedParticles() const;

      void clearAndRecomputeGrid(const BoundingBox &bBox, double gridSpacing);
      
    private:
      typedef boost::shared_array<ParticleVector> ParticleVectorArrayPtr;

      Vec3L                  m_dimensions;
      Vec3L                  m_minIndex;
      Vec3L                  m_maxIndex;
      double                 m_gridSpacing;
      BoundingBox            m_bBox;
      ParticleVector         m_insertedParticles;
      ParticleVectorArrayPtr m_tablePtr;
    };
  }
}

#include "Geometry/NeighbourTable.hpp"

#endif
