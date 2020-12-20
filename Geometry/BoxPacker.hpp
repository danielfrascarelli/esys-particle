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


#include "Geometry/GridIterator.h"
#include <float.h>

namespace esys
{
  namespace lsm
  {
    template <typename TmplPackerBase>
    BoxPacker<TmplPackerBase>::BoxPacker(
      ParticlePoolPtr   particlePoolPtr,
      NTablePtr         nTablePtr,
      const BoundingBox &bBox,
      const BoolVector  &periodicDimensions,
      double            tolerance
    ) : Inherited(particlePoolPtr, nTablePtr),
        m_bBox(bBox),
        m_periodicDimensions(periodicDimensions),
        m_tolerance(tolerance)
    {
    }

    template <typename TmplPackerBase>
    BoxPacker<TmplPackerBase>::~BoxPacker()
    {
    }

    template <typename TmplPackerBase>
    const BoundingBox &BoxPacker<TmplPackerBase>::getBBox() const
    {
      return m_bBox;
    }

    template <typename TmplPackerBase>
    double BoxPacker<TmplPackerBase>::getTolerance() const
    {
      return m_tolerance;
    }

    template <typename TmplPackerBase>
    bool BoxPacker<TmplPackerBase>::is2d() const
    {
      return ((getBBox().getMaxPt().Z() - getBBox().getMinPt().Z()) <= 0);
    }

    template <typename TmplPackerBase>
    const BoolVector &
    BoxPacker<TmplPackerBase>::getPeriodicDimensions() const
    {
      return m_periodicDimensions;
    }

    template <typename TmplPackerBase>
    bool BoxPacker<TmplPackerBase>::particleFitsInBBox(
      const Particle &particle
    ) const
    {
      return
        (
          (
            m_periodicDimensions[0]
            ||
            (
              m_bBox.contains(
                particle.getPos() - Vec3(particle.getRad(), 0, 0),
                getTolerance()
              )
              &&
              m_bBox.contains(
                particle.getPos() + Vec3(particle.getRad(), 0, 0),
                getTolerance()
              )
            )
          )
          &&
          (
            m_periodicDimensions[1]
            ||
            (
              m_bBox.contains(
                particle.getPos() - Vec3(0, particle.getRad(), 0),
                getTolerance()
              )
              &&
              m_bBox.contains(
                particle.getPos() + Vec3(0, particle.getRad(), 0),
                getTolerance()
              )
            )
          )
          &&
          (
            is2d() || m_periodicDimensions[2]
            ||
            (
              m_bBox.contains(
                particle.getPos() - Vec3(0, 0, particle.getRad()),
                getTolerance()
              )
              &&
              m_bBox.contains(
                particle.getPos() + Vec3(0, 0, particle.getRad()),
                getTolerance()
              )
            )
          )
        );
    }

    template <typename TmplPackerBase>
    bool BoxPacker<TmplPackerBase>::particleFitsWithNeighbours(
      const Particle &particle
    ) const
    {
      const typename NTable::ParticleVector neighbours = 
        this->getNTable().getNeighbourVector(
          particle.getPos(),
          particle.getRad() + getTolerance()
        );
      typename NTable::ParticleVector::const_iterator iter = neighbours.begin();
      for (; iter != neighbours.end(); iter++) {
        const double interCentreDistSqrd =
          (particle.getPos() - (*iter)->getPos()).norm2();
        const double radiusSum =
          ((particle.getRad() + (*iter)->getRad()) - getTolerance());
        if (interCentreDistSqrd < (radiusSum*radiusSum)) {
          return false;
        }
      }
      return true;
    }

    template <typename TmplPackerBase>
    bool BoxPacker<TmplPackerBase>::particleFitsInBBoxWithNeighbours(
      const Particle &particle
    ) const
    {
      return
        (
          particleFitsInBBox(particle)
          &&
          particleFitsWithNeighbours(particle)
        );
    }
  }
}
