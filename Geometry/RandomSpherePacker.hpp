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


#include "Foundation/console.h"
#include "Foundation/StringUtil.h"
#include "Geometry/RandomSpherePacker.h"
#include "Geometry/SphereFitter.h"
#include "Geometry/ParticleComparer.h"

#include <algorithm>
#include <stdexcept>
#include <float.h>

namespace esys
{
  namespace lsm
  {
    template<typename TmplTraits>
    SphereFittedPIterator<TmplTraits>::SphereFittedPIterator(
      Packer               &packer,
      int                  maxInsertionFailures,
      const BoundingSphere &bSphere
    ) :
      m_pPacker(&packer),
      m_fitterPtrVector(),
      m_maxInsertionFailures(maxInsertionFailures),
      m_lastFailCount(0),
      m_successCount(0),
      m_next(Fitter::getInvalidParticle()),
      m_bSphere(bSphere)
    {
      if (m_next.isValid())
      {
        throw
          std::runtime_error(
            std::string("isValid method returns true for INVALID particle:")
            +
            StringUtil::toString(m_next)
          );
      }
      initialiseFitterPtrVector();
    }

    template<typename TmplTraits>
    int SphereFittedPIterator<TmplTraits>::getMaxInsertionFailures() const
    {
      return m_maxInsertionFailures;
    }

    template<typename TmplTraits>
    const typename SphereFittedPIterator<TmplTraits>::Packer &
    SphereFittedPIterator<TmplTraits>::getPacker() const
    {
      return *m_pPacker;
    }

    template<typename TmplTraits>
    typename SphereFittedPIterator<TmplTraits>::Packer &
    SphereFittedPIterator<TmplTraits>::getPacker()
    {
      return *m_pPacker;
    }

    template<typename TmplTraits>
    const BoundingSphere &
    SphereFittedPIterator<TmplTraits>::getBSphere() const
    {
      return m_bSphere;
    }

    template<typename TmplTraits>
    void
    SphereFittedPIterator<TmplTraits>::initialiseFitterPtrVector()
    {
      FitterPtrVector fitters;
      fitters.push_back(FitterPtr(new Move2SurfaceFitter(getPacker())));
      if (getPacker().is2d())
      {
        fitters.push_back(FitterPtr(new TwoDSFitter(getPacker())));
        fitters.push_back(FitterPtr(new TwoDSSphereFitter(getPacker(), getBSphere())));
      }
      else
      {
        fitters.push_back(FitterPtr(new ThreeDFitter(getPacker())));
        fitters.push_back(FitterPtr(new ThreeDSSphereFitter(getPacker(), getBSphere())));
      }
      m_fitterPtrVector = fitters;
    }

    template<typename TmplTraits>
    typename SphereFittedPIterator<TmplTraits>::FitterPtrVector &
    SphereFittedPIterator<TmplTraits>::getFitterPtrVector()
    {
      return m_fitterPtrVector;
    }

    template<typename TmplTraits>
    const typename SphereFittedPIterator<TmplTraits>::FitterPtrVector &
    SphereFittedPIterator<TmplTraits>::getFitterPtrVector() const
    {
      return m_fitterPtrVector;
    }

    template<typename TmplTraits>
    Vec3 SphereFittedPIterator<TmplTraits>::getRandomPoint() const
    {
      return getPacker().getRandomPoint();
    }

    template<typename TmplTraits>
    typename SphereFittedPIterator<TmplTraits>::Particle
    SphereFittedPIterator<TmplTraits>::getCandidateParticle(
      const Vec3 &point
    )
    {
      return getPacker().getCandidateParticle(point);
    }

    template<typename TmplTraits>
    typename SphereFittedPIterator<TmplTraits>::ParticleVector
    SphereFittedPIterator<TmplTraits>::getClosestNeighbours(
      const Particle& particle,
      int numClosest
    )
    {
      return getPacker().getClosestNeighbours(particle, numClosest);
    }

    template<typename TmplTraits>
    typename SphereFittedPIterator<TmplTraits>::Particle &
    SphereFittedPIterator<TmplTraits>::generateNext()
    {
      m_next = Fitter::getInvalidParticle();
      if (m_lastFailCount < getMaxInsertionFailures())
      {
        int numFails    = 0;
        //int insertCount = 0;
        FitterPtrVector fitters = getFitterPtrVector();
        Particle particle       = Fitter::getInvalidParticle();
        Particle fitParticle    = particle;
        Plane plane=Plane(Vec3(-1.0, 0, 0), Vec3(DBL_MAX/2, DBL_MAX/2, DBL_MAX/2));
        ParticleVector neighbours;
        while ((!fitParticle.isValid()) && (numFails < getMaxInsertionFailures()))
        {
          particle = getCandidateParticle(getRandomPoint());
          neighbours = getClosestNeighbours(particle, 4);
          
          for (
            typename FitterPtrVector::iterator it = getFitterPtrVector().begin();
            (
              (it != getFitterPtrVector().end())
              &&
              (!fitParticle.isValid())
            );
            it++
          )
          {
            fitParticle = (*it)->getFitParticle(particle, neighbours, plane);
          }
          numFails++;
        }
        m_lastFailCount = numFails;
        m_successCount += ((fitParticle.isValid()) ? 1 : 0);
        m_next = fitParticle;
      }
      return m_next;
    }

    template<typename TmplTraits>
    bool SphereFittedPIterator<TmplTraits>::hasNext()
    {
      return (m_next.isValid() || generateNext().isValid());
    }

    template<typename TmplTraits>
    typename SphereFittedPIterator<TmplTraits>::Particle
    SphereFittedPIterator<TmplTraits>::next()
    {
      if (!(m_next.isValid()))
      {
        generateNext();
      }
      const Particle next = m_next;
      m_next = Fitter::getInvalidParticle();
      return next;
    }

    template<typename TmplTraits>
    void SphereFittedPIterator<TmplTraits>::logInfo()
    {
      console.Info()
        << "BoundingSphere: minPt = " << getPacker().getBBox().getMinPt()
        << ", maxPt = " << getPacker().getBBox().getMaxPt() << "\n"
        << "BoundingSphere: sizes = " << getPacker().getBBox().getSizes() << "\n";

      for (
          typename FitterPtrVector::iterator it = getFitterPtrVector().begin();
          (it != getFitterPtrVector().end());
          it++      
      ) {
        console.Info() << (*(*it)).toString() << "\n";
      }
      console.Info() << "Total successful fits = " << m_successCount << "\n";
    }


























  
    template <typename TPartGen,template <typename TTPartGen> class TCubicPackWrap>
    RandomSpherePacker<TPartGen,TCubicPackWrap>::RandomSpherePacker(
      ParticleGeneratorPtr particleGeneratorPtr,
      ParticlePoolPtr      particlePoolPtr,
      NTablePtr            nTablePtr,
      const BoundingSphere &bSphere,
      double               tolerance,
      double               cubicPackRadius,
      int                  maxInsertionFailures,
      bool                 do2d
    )
     : Inherited(
         particleGeneratorPtr,
         particlePoolPtr,
         nTablePtr,
         do2d ? bSphere.get2dBBox() : bSphere.getBBox(),
         BoolVector(3,false),
         tolerance,
         cubicPackRadius
       ),
       m_bSphere(bSphere),
       m_maxInsertionFailures(maxInsertionFailures)
    {
    }

    template <typename TPartGen,template <typename TTPartGen> class TCubicPackWrap>
    RandomSpherePacker<TPartGen,TCubicPackWrap>::~RandomSpherePacker()
    {
    }

    template <typename TPartGen,template <typename TTPartGen> class TCubicPackWrap>
    const BoundingSphere &
    RandomSpherePacker<TPartGen,TCubicPackWrap>::getBSphere(
    ) const
    {
      return m_bSphere;
    }

    template <typename TPartGen,template <typename TTPartGen> class TCubicPackWrap>
    double
    RandomSpherePacker<TPartGen,TCubicPackWrap>::getRandom(
      double minVal,
      double maxVal
    ) const
    {
      return minVal + (maxVal-minVal)*rng::s_zeroOneUniform();
    }

    template <typename TPartGen,template <typename TTPartGen> class TCubicPackWrap>
    Vec3 RandomSpherePacker<TPartGen,TCubicPackWrap>::getRandomPoint() const
    {
      return 
        Vec3(
          getRandom(this->getBBox().getMinPt().X(), this->getBBox().getMaxPt().X()),
          getRandom(this->getBBox().getMinPt().Y(), this->getBBox().getMaxPt().Y()),
          getRandom(this->getBBox().getMinPt().Z(), this->getBBox().getMaxPt().Z())
        );
    }

    template <typename TPartGen,template <typename TTPartGen> class TCubicPackWrap>
    typename RandomSpherePacker<TPartGen,TCubicPackWrap>::ParticleVector
    RandomSpherePacker<TPartGen,TCubicPackWrap>::getClosestNeighbours(
      const Particle& particle,
      int numClosest
    )
    {
      ParticleVector neighbourVector =
        this->getNTable().getUniqueNeighbourVector(
          particle.getPos(),
          this->getParticleGenerator().getMaxFitRadius() + this->getTolerance()
        );
      std::sort(
        neighbourVector.begin(),
        neighbourVector.end(),
        ParticleComparer<Particle>(particle)
      );

      if (static_cast<int>(neighbourVector.size()) > numClosest) {
        neighbourVector.erase(
          neighbourVector.begin() + numClosest,
          neighbourVector.end()
        );
      }

      return neighbourVector;
    }

    template <typename TPartGen,template <typename TTPartGen> class TCubicPackWrap>
    bool RandomSpherePacker<TPartGen,TCubicPackWrap >::particleIsValid(
      const Particle &particle
    ) const
    {
      return 
        (
          this->getParticleGenerator().isValidFitRadius(particle.getRad())
          &&
          particleFitsInBSphereWithNeighbours(particle)
        );
    }

    template <typename TPartGen,template <typename TTPartGen> class TCubicPackWrap>
    int RandomSpherePacker<TPartGen,TCubicPackWrap>::getMaxInsertionFailures() const
    {
      return m_maxInsertionFailures;
    }

    template <typename TPartGen,template <typename TTPartGen> class TCubicPackWrap>
    void RandomSpherePacker<TPartGen,TCubicPackWrap>::generateRandomFill()
    {
      StuffedParticleIterator it =
        StuffedParticleIterator(
          *this,
          getMaxInsertionFailures(),
          getBSphere()
        );
      while (it.hasNext())
      {
        const Particle p = it.next();
        this->createAndInsertParticle(p);
      }
      it.logInfo();
    }

    template <typename TPartGen,template <typename TTPartGen> class TCubicPackWrap>
    bool
    RandomSpherePacker<TPartGen,TCubicPackWrap>::particleFitsInBSphere(
      const Particle &particle
    ) const
    {
      return
        getBSphere().contains(
          BoundingSphere(particle.getPos(), particle.getRadius()),
          this->getTolerance()
        );
    }

    template <typename TPartGen,template <typename TTPartGen> class TCubicPackWrap>
    bool
    RandomSpherePacker<TPartGen,TCubicPackWrap>::particleFitsInBSphereWithNeighbours(
      const Particle &particle
    ) const
    {
      return
        (
          particleFitsInBSphere(particle)
          &&
          this->particleFitsWithNeighbours(particle)
        );
    }

    template <typename TPartGen,template <typename TTPartGen> class TCubicPackWrap>
    void RandomSpherePacker<TPartGen,TCubicPackWrap>::generateCubicPackingInSphere()
    {
      GridIterator pointIt = GridIterator(this->getBBox(), this->getCubicPackingRadius());
      while (pointIt.hasNext()) {
        const Particle candidate =
          this->getCandidateParticle(pointIt.next(), this->getCubicPackingRadius());
        if (particleFitsInBSphereWithNeighbours(candidate)) {
          this->createAndInsertParticle(candidate);
        }
      }
    }

    template <typename TPartGen,template <typename TTPartGen> class TCubicPackWrap>
    void RandomSpherePacker<TPartGen,TCubicPackWrap>::generate()
    {
      generateCubicPackingInSphere();
      generateRandomFill();
    }
  }
}
