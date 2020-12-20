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


namespace esys
{
  namespace lsm
  {
    template <typename TmplParticleGenerator, typename TmplBoxPackerBase>
    CubicBoxPacker<TmplParticleGenerator,TmplBoxPackerBase>::CubicBoxPacker(
      ParticleGeneratorPtr particleGeneratorPtr,
      ParticlePoolPtr      particlePoolPtr,
      NTablePtr            nTablePtr,
      const BoundingBox    &bBox,
      const BoolVector     &periodicDimensions,
      double               tolerance,
      double               cubicPackRadius
    ) : Inherited(
          particlePoolPtr,
          nTablePtr,
          bBox,
          periodicDimensions,
          tolerance
        ),
        m_cubicPackRadius(cubicPackRadius),
        m_particleGeneratorPtr(particleGeneratorPtr),
        m_pParticleGenerator(particleGeneratorPtr.get())
    {
    }

    template <typename TmplParticleGenerator, typename TmplBoxPackerBase>
    CubicBoxPacker<TmplParticleGenerator,TmplBoxPackerBase>::~CubicBoxPacker()
    {
    }

    template <typename TmplParticleGenerator, typename TmplBoxPackerBase>
    const typename CubicBoxPacker<TmplParticleGenerator,TmplBoxPackerBase>::ParticleGenerator &
    CubicBoxPacker<TmplParticleGenerator,TmplBoxPackerBase>::getParticleGenerator() const
    {
      return *m_pParticleGenerator;
    }

    template <typename TmplParticleGenerator, typename TmplBoxPackerBase>
    typename CubicBoxPacker<TmplParticleGenerator,TmplBoxPackerBase>::ParticleGenerator &
    CubicBoxPacker<TmplParticleGenerator,TmplBoxPackerBase>::getParticleGenerator()
    {
      return *m_pParticleGenerator;
    }

    template <typename TmplParticleGenerator, typename TmplBoxPackerBase>
    void
    CubicBoxPacker<TmplParticleGenerator,TmplBoxPackerBase>::setParticleGenerator(
      ParticleGenerator &particleGenerator
    )
    {
      m_pParticleGenerator = &particleGenerator;
      m_particleGeneratorPtr = ParticleGeneratorPtr();
    }

    template <typename TmplParticleGenerator, typename TmplBoxPackerBase>
    void
    CubicBoxPacker<TmplParticleGenerator,TmplBoxPackerBase>::setParticleGenerator(
      ParticleGeneratorPtr particleGeneratorPtr
    )
    {
      m_particleGeneratorPtr = particleGeneratorPtr;
      m_pParticleGenerator = m_particleGeneratorPtr.get();
    }

    template <typename TmplParticleGenerator, typename TmplBoxPackerBase>
    double
    CubicBoxPacker<TmplParticleGenerator,TmplBoxPackerBase>::getCubicPackingRadius() const
    {
      return m_cubicPackRadius;
    }

    template <typename TmplParticleGenerator, typename TmplBoxPackerBase>
    typename CubicBoxPacker<TmplParticleGenerator,TmplBoxPackerBase>::Particle
    CubicBoxPacker<TmplParticleGenerator,TmplBoxPackerBase>::getCandidateParticle(
      const Vec3 &point,
      double radius
    )
    {
      return getParticleGenerator().getParticle(point, radius);
    }

    template <typename TmplParticleGenerator, typename TmplBoxPackerBase>
    typename CubicBoxPacker<TmplParticleGenerator,TmplBoxPackerBase>::Particle
    CubicBoxPacker<TmplParticleGenerator,TmplBoxPackerBase>::getCandidateParticle(
      const Vec3 &point
    )
    {
      return getParticleGenerator().getParticle(point);
    }

    template <typename TmplParticleGenerator, typename TmplBoxPackerBase>
    void
    CubicBoxPacker<TmplParticleGenerator,TmplBoxPackerBase>::generateCubicPacking()
    {
      GridIterator pointIt = GridIterator(this->getBBox(), getCubicPackingRadius());
      while (pointIt.hasNext()) {
        const Particle candidate =
          getCandidateParticle(pointIt.next(), getCubicPackingRadius());
        if (this->particleFitsInBBoxWithNeighbours(candidate)) {
          this->createAndInsertParticle(candidate);
        }
      }
    }

    template <typename TmplParticleGenerator, typename TmplBoxPackerBase>
    void
    CubicBoxPacker<TmplParticleGenerator,TmplBoxPackerBase>::generate()
    {
      generateCubicPacking();
    }
  }
}
