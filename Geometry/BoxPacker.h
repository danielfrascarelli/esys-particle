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


#ifndef ESYS_LSMBOXPACKER_H
#define ESYS_LSMBOXPACKER_H

#include <Geometry/Packer.h>
#include <Foundation/vec3.h>
#include <Foundation/BoundingBox.h>

#include <vector>

namespace esys
{
  namespace lsm
  {
    typedef std::vector<bool> BoolVector;
    /**
     *
     */
    template <typename TmplPackerBase>
    class BoxPacker : public TmplPackerBase
    {
    public:
      typedef TmplPackerBase                      Inherited;
      typedef typename Inherited::Particle        Particle;
      typedef typename Inherited::NTable          NTable;
      typedef typename Inherited::NTablePtr       NTablePtr;
      typedef typename Inherited::ParticlePool    ParticlePool;
      typedef typename Inherited::ParticlePoolPtr ParticlePoolPtr;

      BoxPacker(
        ParticlePoolPtr   particlePoolPtr,
        NTablePtr         nTablePtr,
        const BoundingBox &bBox,
        const BoolVector  &periodicDimensions,
        double            tolerance
      );

      virtual ~BoxPacker();

      virtual void generate() = 0;

      bool particleFitsInBBox(const Particle &particle) const;

      bool is2d() const;

      bool particleFitsInBBoxWithNeighbours(const Particle &particle) const;

      bool particleFitsWithNeighbours(const Particle &particle) const;

      double getTolerance() const;

      const BoundingBox &getBBox() const;


      const BoolVector &getPeriodicDimensions() const;

    private:
      BoundingBox    m_bBox;
      BoolVector     m_periodicDimensions;
      double         m_tolerance;
    };
  }
}

#include "Geometry/BoxPacker.hpp"

#endif
