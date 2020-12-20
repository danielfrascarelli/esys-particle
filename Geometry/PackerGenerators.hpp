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


#include "Foundation/Rng.h"

namespace esys
{
  namespace lsm
  {
    template <typename TmplParticle>
    ConstRadiusGen<TmplParticle>::ConstRadiusGen(double radius)
      : m_radius(radius)
    {
    }

    template <typename TmplParticle>
    const double &ConstRadiusGen<TmplParticle>::getParticleRadius() const
    {
      return m_radius;
    }

    template <typename TmplParticle>
    const double &ConstRadiusGen<TmplParticle>::getMinFitRadius() const
    {
      return this->getParticleRadius();
    }

    template <typename TmplParticle>
    const double &ConstRadiusGen<TmplParticle>::getMaxFitRadius() const
    {
      return this->getParticleRadius();
    }

    template <typename TmplParticle>
    typename ConstRadiusGen<TmplParticle>::Particle
    ConstRadiusGen<TmplParticle>::getParticle(const Vec3 &posn) const
    {
      return Particle(posn, getParticleRadius());
    }

    template <typename TmplParticle>
    typename ConstRadiusGen<TmplParticle>::Particle
    ConstRadiusGen<TmplParticle>::getParticle(
      const Vec3 &posn,
      const double &maxRadius
    ) const
    {
      return getParticle(posn);
    }

    template <typename TmplParticle>
    bool ConstRadiusGen<TmplParticle>::isValidFitRadius(
      const double &fitRadius
    ) const
    {
      return (fitRadius == getParticleRadius());
    }
    //========================================================================
    //========================================================================
    //========================================================================
    template <typename TmplParticle>
    RangeRadiusGen<TmplParticle>::RangeRadiusGen(
      double minFitRadius,
      double maxFitRadius
    )
      : m_minFitRadius(minFitRadius),
        m_maxFitRadius(maxFitRadius)
    {
    }

    template <typename TmplParticle>
    RangeRadiusGen<TmplParticle>::~RangeRadiusGen()
    {
    }

    template <typename TmplParticle>
    const double &RangeRadiusGen<TmplParticle>::getMinFitRadius() const
    {
      return m_minFitRadius;
    }

    template <typename TmplParticle>
    const double &RangeRadiusGen<TmplParticle>::getMaxFitRadius() const
    {
      return m_maxFitRadius;
    }

    template <typename TmplParticle>
    bool RangeRadiusGen<TmplParticle>::isValidFitRadius(
      const double &fitRadius
    ) const
    {
      return
        (
          (fitRadius >= this->getMinFitRadius())
          &&
          (fitRadius <= this->getMaxFitRadius())
        );
    }
    //========================================================================
    //========================================================================
    //========================================================================
    template <typename TmplParticle>
    RndRadiusGen<TmplParticle>::RndRadiusGen(
      double minFitRadius,
      double maxFitRadius
    ) : Inherited(minFitRadius, maxFitRadius)
    {
    }

    template <typename TmplParticle>
    double RndRadiusGen<TmplParticle>::getRandomRadius() const
    {
      return
        this->getMinFitRadius()
        +
        (this->getMaxFitRadius()-this->getMinFitRadius())
        *
        rng::s_zeroOneUniform();
    }

    template <typename TmplParticle>
    typename RndRadiusGen<TmplParticle>::Particle
    RndRadiusGen<TmplParticle>::getParticle(const Vec3 &posn) const
    {
      return Particle(posn, getRandomRadius());
    }

    template <typename TmplParticle>
    typename RndRadiusGen<TmplParticle>::Particle
    RndRadiusGen<TmplParticle>::getParticle(
      const Vec3 &posn,
      double suggestedRadius
    ) const
    {
      return getParticle(posn);
    }
    //========================================================================
    //========================================================================
    //========================================================================
    template <typename TmplGrain>
    GrainRndRadiusGen<TmplGrain>::GrainRndRadiusGen(
      double minGrainRadius,
      double maxGrainRadius
    ) : Inherited(minGrainRadius, maxGrainRadius)
    {
    }

    template <typename TmplGrain>
    GrainRndRadiusGen<TmplGrain>::~GrainRndRadiusGen()
    {
    }

    template <typename TmplGrain>
    const double &GrainRndRadiusGen<TmplGrain>::getMinGrainRadius() const
    {
      return this->getMinFitRadius();
    }

    template <typename TmplGrain>
    const double &GrainRndRadiusGen<TmplGrain>::getMaxGrainRadius() const
    {
      return this->getMaxFitRadius();
    }

    //========================================================================
    //========================================================================
    //========================================================================

    template <typename TmplGrain>
    SingleParticleGrainGen<TmplGrain>::SingleParticleGrainGen(
      double minGrainRadius,
      double maxGrainRadius
    ) : Inherited(minGrainRadius, maxGrainRadius)
    {
    }

    template <typename TmplGrain>
    const double &SingleParticleGrainGen<TmplGrain>::getMinParticleRadius() const
    {
      return this->getMinGrainRadius();
    }

    template <typename TmplGrain>
    const double &SingleParticleGrainGen<TmplGrain>::getMaxParticleRadius() const
    {
      return this->getMaxGrainRadius();
    }

    template <typename TmplGrain>
    typename SingleParticleGrainGen<TmplGrain>::Grain
    SingleParticleGrainGen<TmplGrain>::getGrain(const Particle &p)
    {
      Grain g;
      g.createParticle(p);
      return g;
    }
  }
}
