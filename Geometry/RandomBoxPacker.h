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


#ifndef ESYS_LSMRANDOMBOXPACKER_H
#define ESYS_LSMRANDOMBOXPACKER_H

#include "Foundation/console.h"
#include "Foundation/Rng.h"
#include "Geometry/CubicBoxPacker.h"
#include "Geometry/SphereFitter.h"
#include "Geometry/Plane3D.h"

#include <vector>
#include <boost/shared_ptr.hpp>

namespace esys
{
  namespace lsm
  {
    template <typename TmplFitterTraits>
    class FittedParticleIterator
    {
    public:
      typedef TmplFitterTraits                   FitterTraits;
      typedef typename FitterTraits::Plane3D       Plane3D;
      typedef typename FitterTraits::PlaneVector PlaneVector;
      typedef typename FitterTraits::Packer      Packer;
      typedef typename Packer::Particle          Particle;
      typedef typename Packer::ParticleVector    ParticleVector;
      
      class FitTraits
      {
      public:
        typedef Packer                           Validator;
        typedef typename Packer::Particle        Particle;
        typedef typename Packer::ParticleVector  ParticleVector;
        typedef typename FitterTraits::Plane3D     Plane3D;
      };

      typedef SphereFitter<FitTraits>            Fitter;
      typedef boost::shared_ptr<Fitter>          FitterPtr;
      typedef std::vector<FitterPtr>             FitterPtrVector;
      typedef MoveToSurfaceFitter<FitTraits>     Move2SurfaceFitter;
      typedef ThreeDSphereFitter<FitTraits>      ThreeDFitter;
      typedef TwoDSphereFitter<FitTraits>        TwoDFitter;
      typedef TwoDPlaneSphereFitter<FitTraits>   TwoDPlaneFitter;
      typedef ThreeDPlaneSphereFitter<FitTraits> ThreeDPlaneFitter;

      FittedParticleIterator(
        Packer            &packer,
        int               maxInsertionFailures,
        const PlaneVector &fitPlaneVector
      );

      void initialiseFitterPtrVector();

      int getMaxInsertionFailures() const;

      const FitterPtrVector &getFitterPtrVector() const;

      FitterPtrVector &getFitterPtrVector();

      const PlaneVector &getFitPlaneVector() const;

      const Packer &getPacker() const;

      Packer &getPacker();

      Plane3D getClosestFitPlane(const Particle &particle) const;

      double getRandom(double min, double max) const;

      Vec3 getRandomPoint() const;

      Particle getCandidateParticle(const Vec3 &point);

      ParticleVector getClosestNeighbours(const Particle& particle, int numClosest);

      Particle &generateNext();

      bool hasNext();

      Particle next();

      void logInfo();

    private:
      Packer          *m_pPacker;
      PlaneVector     m_fitPlaneVector;
      int             m_maxInsertionFailures;      
      int             m_lastFailCount;
      int             m_successCount;
      Particle        m_next;
      FitterPtrVector m_fitterPtrVector;
    };

    /**
     *
     */
    template <typename TmplParticleGenerator,template <typename TmplPartGen> class TmplCubicBoxPackerWrap>
    class RandomBoxPacker : public TmplCubicBoxPackerWrap<TmplParticleGenerator>::CubicBoxPackerBase
    {
    public:
      typedef
        typename TmplCubicBoxPackerWrap<TmplParticleGenerator>::CubicBoxPackerBase
        Inherited;
      typedef typename Inherited::ParticleGenerator     ParticleGenerator;
      typedef typename Inherited::ParticleGeneratorPtr  ParticleGeneratorPtr;
      typedef typename Inherited::Particle              Particle;
      typedef typename Inherited::NTable                NTable;
      typedef typename Inherited::NTablePtr             NTablePtr;
      typedef typename NTable::ParticleVector           ParticleVector;
      typedef typename Inherited::ParticlePool          ParticlePool;
      typedef typename Inherited::ParticlePoolPtr       ParticlePoolPtr;
      typedef std::vector<Plane3D>                        PlaneVector;

      class StufferTraits
      {
      public:
        typedef RandomBoxPacker              Packer;
        typedef esys::lsm::Plane3D             Plane3D;
        typedef typename Packer::PlaneVector PlaneVector;
      };
      typedef FittedParticleIterator<StufferTraits> StuffedParticleIterator;

      RandomBoxPacker(
        ParticleGeneratorPtr particleGeneratorPtr,
        ParticlePoolPtr      particlePoolPtr,
        NTablePtr            nTablePtr,
        const BoundingBox    &bBox,
        const BoolVector     &periodicDimensions,
        double               tolerance,
        double               cubicPackRadius,
        int                  maxInsertionFailures
      );

      RandomBoxPacker(
        ParticleGeneratorPtr particleGeneratorPtr,
        ParticlePoolPtr      particlePoolPtr,
        NTablePtr            nTablePtr,
        const BoundingBox    &bBox,
        const BoolVector     &periodicDimensions,
        double               tolerance,
        double               cubicPackRadius,
        int                  maxInsertionFailures,
        const PlaneVector    &fitPlaneVector
      );

      virtual ~RandomBoxPacker();

      PlaneVector getDefaultFitPlaneVector() const;

      bool particleIsValid(const Particle &particle) const;

      virtual void generate();

      double getRandom(double min, double max) const;

      Vec3 getRandomPoint() const;

      ParticleVector getClosestNeighbours(const Particle& particle, int numClosest);

      int getMaxInsertionFailures() const;

      void generateRandomFill();

      const PlaneVector &getFitPlaneVector() const;

      Plane3D getClosestFitPlane(const Particle &particle) const;

    private:
      PlaneVector  m_fitPlaneVector;
      int          m_maxInsertionFailures;
    };
  }
}

#include "Geometry/RandomBoxPacker.hpp"

#endif
