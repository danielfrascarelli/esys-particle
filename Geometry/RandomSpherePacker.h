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


#ifndef ESYS_LSMRANDOMSPHEREPACKER_H
#define ESYS_LSMRANDOMSPHEREPACKER_H

#include "Foundation/console.h"
#include "Foundation/Rng.h"
#include "Foundation/BoundingSphere.h"
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
    class SphereFittedPIterator
    {
    public:
      typedef TmplFitterTraits                   FitterTraits;
      typedef typename FitterTraits::Plane3D       Plane;
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

      typedef SphereFitter<FitTraits>             Fitter;
      typedef boost::shared_ptr<Fitter>           FitterPtr;
      typedef std::vector<FitterPtr>              FitterPtrVector;
      typedef MoveToSurfaceFitter<FitTraits>      Move2SurfaceFitter;
      typedef ThreeDSphereFitter<FitTraits>       ThreeDFitter;
      typedef TwoDSphereFitter<FitTraits>         TwoDSFitter;
      typedef TwoDSphereSphereFitter<FitTraits>   TwoDSSphereFitter;
      typedef ThreeDSphereSphereFitter<FitTraits> ThreeDSSphereFitter;

      SphereFittedPIterator(
        Packer               &packer,
        int                  maxInsertionFailures,
        const BoundingSphere &bSphere
      );

      void initialiseFitterPtrVector();

      const BoundingSphere &getBSphere() const;

      int getMaxInsertionFailures() const;

      const FitterPtrVector &getFitterPtrVector() const;

      FitterPtrVector &getFitterPtrVector();

      const Packer &getPacker() const;

      Packer &getPacker();

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
      FitterPtrVector m_fitterPtrVector;
      int             m_maxInsertionFailures;
      int             m_lastFailCount;
      int             m_successCount;
      Particle        m_next;
      BoundingSphere  m_bSphere;
    };

    /**
     *
     */
    template <typename TmplParticleGenerator,template <typename TmplPartGen> class TmplCubicBoxPackerWrap>
    class RandomSpherePacker : public TmplCubicBoxPackerWrap<TmplParticleGenerator>::CubicBoxPackerBase
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

      class StufferTraits
      {
      public:
        typedef RandomSpherePacker Packer;
        typedef esys::lsm::Plane3D   Plane3D;
        typedef std::vector<Plane3D> PlaneVector;
      };
      typedef SphereFittedPIterator<StufferTraits> StuffedParticleIterator;

      RandomSpherePacker(
        ParticleGeneratorPtr particleGeneratorPtr,
        ParticlePoolPtr      particlePoolPtr,
        NTablePtr            nTablePtr,
        const BoundingSphere &bSphere,
        double               tolerance,
        double               cubicPackRadius,
        int                  maxInsertionFailures,
        bool                 do2d
      );

      virtual ~RandomSpherePacker();

      const BoundingSphere &getBSphere() const;

      bool particleIsValid(const Particle &particle) const;

      double getRandom(double min, double max) const;

      Vec3 getRandomPoint() const;

      ParticleVector getClosestNeighbours(const Particle& particle, int numClosest);

      int getMaxInsertionFailures() const;

      bool particleFitsInBSphere(const Particle &particle) const;

      bool particleFitsInBSphereWithNeighbours(const Particle &particle) const;

      void generateCubicPackingInSphere();

      void generateRandomFill();

      virtual void generate();

    private:
      BoundingSphere m_bSphere;
      int            m_maxInsertionFailures;
    };
  }
}

#include "Geometry/RandomSpherePacker.hpp"

#endif
