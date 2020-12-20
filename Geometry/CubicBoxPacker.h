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


#ifndef ESYS_LSMCUBICBOXPACKER_H
#define ESYS_LSMCUBICBOXPACKER_H

#include <Geometry/BoxPacker.h>

namespace esys
{
  namespace lsm
  {
    typedef std::vector<bool> BoolVector;
    /**
     *
     */
    template <typename TmplParticleGenerator, typename TmplBoxPackerBase>
    class CubicBoxPacker : public TmplBoxPackerBase
    {
    public:
      typedef TmplParticleGenerator                ParticleGenerator;
      typedef boost::shared_ptr<ParticleGenerator> ParticleGeneratorPtr;
      typedef TmplBoxPackerBase                    Inherited;
      typedef Inherited                            BoxPackerBase;
      typedef typename Inherited::Particle         Particle;
      typedef typename Inherited::NTable           NTable;
      typedef typename Inherited::NTablePtr        NTablePtr;
      typedef typename Inherited::ParticlePool     ParticlePool;
      typedef typename Inherited::ParticlePoolPtr  ParticlePoolPtr;

      CubicBoxPacker(
        ParticleGeneratorPtr particleGeneratorPtr,
        ParticlePoolPtr      particlePoolPtr,
        NTablePtr            nTablePtr,
        const BoundingBox    &bBox,
        const BoolVector     &periodicDimensions,
        double               tolerance,
        double               cubicPackRadius
      );

      virtual ~CubicBoxPacker();

      Particle getCandidateParticle(const Vec3 &point);

      Particle getCandidateParticle(const Vec3 &point, double radius);

      double getCubicPackingRadius() const;

      const ParticleGenerator &getParticleGenerator() const;

      ParticleGenerator &getParticleGenerator();

      void setParticleGenerator(ParticleGenerator &particleGenerator);

      void setParticleGenerator(ParticleGeneratorPtr particleGenerator);

      void generateCubicPacking();

      virtual void generate();

    private:
      double               m_cubicPackRadius;
      ParticleGeneratorPtr m_particleGeneratorPtr;
      ParticleGenerator    *m_pParticleGenerator;
    };
  }
}

#include "Geometry/CubicBoxPacker.hpp"

#endif
