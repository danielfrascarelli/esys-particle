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
#include "Geometry/RandomBoxPacker.h"
#include "Geometry/ParticleComparer.h"
#include "Geometry/SphereFitter.h"

#include <algorithm>
#include <stdexcept>
#include <float.h>

namespace esys
{
  namespace lsm
  {
    template <typename TmplParticle>
    class PlaneComparer
    {
    public:
      typedef TmplParticle Particle;
      PlaneComparer(const Particle &particle)
        : m_pParticle(&particle)
      {
      }
      
      inline bool operator()(const Plane3D &plane1, const Plane3D &plane2) const
      {
        return 
          (
            plane1.sep(m_pParticle->getPos())
            <
            plane2.sep(m_pParticle->getPos())
          );
      }
    private:
      const Particle *m_pParticle;
    };

    template<typename TmplTraits>
    FittedParticleIterator<TmplTraits>::FittedParticleIterator(
      Packer            &packer,
      int               maxInsertionFailures,
      const PlaneVector &fitPlaneVector
    ) :
      m_pPacker(&packer),
      m_fitPlaneVector(fitPlaneVector),
      m_maxInsertionFailures(maxInsertionFailures),
      m_lastFailCount(0),
      m_successCount(0),
      m_next(Fitter::getInvalidParticle()),
      m_fitterPtrVector()
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
    const typename FittedParticleIterator<TmplTraits>::PlaneVector &
    FittedParticleIterator<TmplTraits>::getFitPlaneVector() const
    {
      return m_fitPlaneVector;
    }

    template<typename TmplTraits>
    int FittedParticleIterator<TmplTraits>::getMaxInsertionFailures() const
    {
      return m_maxInsertionFailures;
    }

    template<typename TmplTraits>
    const typename FittedParticleIterator<TmplTraits>::Packer &
    FittedParticleIterator<TmplTraits>::getPacker() const
    {
      return *m_pPacker;
    }

    template<typename TmplTraits>
    typename FittedParticleIterator<TmplTraits>::Packer &
    FittedParticleIterator<TmplTraits>::getPacker()
    {
      return *m_pPacker;
    }

    template<typename TmplTraits>
    void
    FittedParticleIterator<TmplTraits>::initialiseFitterPtrVector()
    {
      FitterPtrVector fitters;
      fitters.push_back(FitterPtr(new Move2SurfaceFitter(getPacker())));
      if (getPacker().is2d())
      {
        fitters.push_back(FitterPtr(new TwoDFitter(getPacker())));
        if (getFitPlaneVector().size() > 0) {
          fitters.push_back(FitterPtr(new TwoDPlaneFitter(getPacker())));
        }
      }
      else
      {
        fitters.push_back(FitterPtr(new ThreeDFitter(getPacker())));
        if (getFitPlaneVector().size() > 0) {
          fitters.push_back(FitterPtr(new ThreeDPlaneFitter(getPacker())));
        }
      }

      m_fitterPtrVector = fitters;
    }

    template<typename TmplTraits>
    typename FittedParticleIterator<TmplTraits>::FitterPtrVector &
    FittedParticleIterator<TmplTraits>::getFitterPtrVector()
    {
      return m_fitterPtrVector;
    }

    template<typename TmplTraits>
    const typename FittedParticleIterator<TmplTraits>::FitterPtrVector &
    FittedParticleIterator<TmplTraits>::getFitterPtrVector() const
    {
      return m_fitterPtrVector;
    }

    template<typename TmplTraits>
    typename FittedParticleIterator<TmplTraits>::Plane3D
    FittedParticleIterator<TmplTraits>::getClosestFitPlane(
      const Particle &particle
    ) const
    {
      PlaneVector fitPlanes = getFitPlaneVector();
      if (fitPlanes.size() > 0) {
        std::sort(fitPlanes.begin(), fitPlanes.end(), PlaneComparer<Particle>(particle));

        return fitPlanes[0];
      }
      return Plane3D(Vec3(-1.0, 0, 0), Vec3(DBL_MAX/2, DBL_MAX/2, DBL_MAX/2));
    }

    template<typename TmplTraits>
    Vec3 FittedParticleIterator<TmplTraits>::getRandomPoint() const
    {
      return getPacker().getRandomPoint();
    }

    template<typename TmplTraits>
    typename FittedParticleIterator<TmplTraits>::Particle
    FittedParticleIterator<TmplTraits>::getCandidateParticle(
      const Vec3 &point
    )
    {
      return getPacker().getCandidateParticle(point);
    }

    template<typename TmplTraits>
    typename FittedParticleIterator<TmplTraits>::ParticleVector
    FittedParticleIterator<TmplTraits>::getClosestNeighbours(
      const Particle& particle,
      int numClosest
    )
    {
      return getPacker().getClosestNeighbours(particle, numClosest);
    }

    template<typename TmplTraits>
    typename FittedParticleIterator<TmplTraits>::Particle &
    FittedParticleIterator<TmplTraits>::generateNext()
    {
      m_next = Fitter::getInvalidParticle();
      if (m_lastFailCount < getMaxInsertionFailures())
      {
        int numFails    = 0;
        //int insertCount = 0;
        FitterPtrVector fitters = getFitterPtrVector();
        Particle particle       = Fitter::getInvalidParticle();
        Particle fitParticle    = particle;
        Plane3D plane;
        ParticleVector neighbours;
        while ((!fitParticle.isValid()) && (numFails < getMaxInsertionFailures()))
        {
          particle = getCandidateParticle(getRandomPoint());
          neighbours = getClosestNeighbours(particle, 4);
          plane = getClosestFitPlane(particle);
          
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
        console.Info()
          << "FittedParticleIterator<T>::generateNext: numFails="
          << numFails << "\n";
        m_successCount += ((fitParticle.isValid()) ? 1 : 0);
        m_next = fitParticle;
      }
      return m_next;
    }

    template<typename TmplTraits>
    bool FittedParticleIterator<TmplTraits>::hasNext()
    {
      return (m_next.isValid() || generateNext().isValid());
    }

    template<typename TmplTraits>
    typename FittedParticleIterator<TmplTraits>::Particle
    FittedParticleIterator<TmplTraits>::next()
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
    void FittedParticleIterator<TmplTraits>::logInfo()
    {
      console.Info()
        << "BoundingBox: minPt = " << getPacker().getBBox().getMinPt()
        << ", maxPt = " << getPacker().getBBox().getMaxPt() << "\n"
        << "BoundingBox: sizes = " << getPacker().getBBox().getSizes() << "\n";

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
    RandomBoxPacker<TPartGen,TCubicPackWrap>::RandomBoxPacker(
      ParticleGeneratorPtr particleGeneratorPtr,
      ParticlePoolPtr      particlePoolPtr,
      NTablePtr            nTablePtr,
      const BoundingBox    &bBox,
      const BoolVector     &periodicDimensions,
      double               tolerance,
      double               cubicPackRadius,
      int                  maxInsertionFailures,
      const PlaneVector    &fitPlaneVector
    )
     : Inherited(
         particleGeneratorPtr,
         particlePoolPtr,
         nTablePtr,
         bBox,
         periodicDimensions,
         tolerance,
         cubicPackRadius
       ),
       m_fitPlaneVector(fitPlaneVector),
       m_maxInsertionFailures(maxInsertionFailures)
    {
    }

    template <typename TPartGen,template <typename TTPartGen> class TCubicPackWrap>
    RandomBoxPacker<TPartGen,TCubicPackWrap>::RandomBoxPacker(
      ParticleGeneratorPtr particleGeneratorPtr,
      ParticlePoolPtr      particlePoolPtr,
      NTablePtr            nTablePtr,
      const BoundingBox    &bBox,
      const BoolVector     &periodicDimensions,
      double               tolerance,
      double               cubicPackRadius,
      int                  maxInsertionFailures
    )
     : Inherited(
         particleGeneratorPtr,
         particlePoolPtr,
         nTablePtr,
         bBox,
         periodicDimensions,
         tolerance,
         cubicPackRadius
       ),
       m_fitPlaneVector(),
       m_maxInsertionFailures(maxInsertionFailures)
    {
      m_fitPlaneVector = getDefaultFitPlaneVector();
    }

    template <typename TPartGen,template <typename TTPartGen> class TCubicPackWrap>
    typename RandomBoxPacker<TPartGen,TCubicPackWrap>::PlaneVector
    RandomBoxPacker<TPartGen,TCubicPackWrap>::getDefaultFitPlaneVector() const
    {
      PlaneVector fitPlaneVector;
      if ((!this->getPeriodicDimensions()[1])) {
        fitPlaneVector.push_back(
          Plane3D(Vec3(0,  1, 0), this->getBBox().getMinPt())
        );
        fitPlaneVector.push_back(
          Plane3D(Vec3(0, -1, 0), this->getBBox().getMaxPt())
        );
      }
      if ((!this->getPeriodicDimensions()[0])) {
        fitPlaneVector.push_back(
          Plane3D(Vec3( 1, 0, 0), this->getBBox().getMinPt())
        );
        fitPlaneVector.push_back(
          Plane3D(Vec3(-1, 0, 0), this->getBBox().getMaxPt())
        );
      }
      if (
          !this->is2d()
          &&
          (!this->getPeriodicDimensions()[2])
      ) {
        fitPlaneVector.push_back(
          Plane3D(Vec3(0, 0,  1), this->getBBox().getMinPt())
        );
        fitPlaneVector.push_back(
          Plane3D(Vec3(0, 0, -1), this->getBBox().getMaxPt())
        );
      }
      return fitPlaneVector;
    }

    template <typename TPartGen,template <typename TTPartGen> class TCubicPackWrap>
    RandomBoxPacker<TPartGen,TCubicPackWrap>::~RandomBoxPacker()
    {
    }

    template <typename TPartGen,template <typename TTPartGen> class TCubicPackWrap>
    double
    RandomBoxPacker<TPartGen,TCubicPackWrap>::getRandom(
      double minVal,
      double maxVal
    ) const
    {
      return minVal + (maxVal-minVal)*rng::s_zeroOneUniform();
    }

    template <typename TPartGen,template <typename TTPartGen> class TCubicPackWrap>
    const typename RandomBoxPacker<TPartGen,TCubicPackWrap>::PlaneVector &
    RandomBoxPacker<TPartGen,TCubicPackWrap>::getFitPlaneVector() const
    {
      return m_fitPlaneVector;
    }


    template <typename TPartGen,template <typename TTPartGen> class TCubicPackWrap>
    Vec3 RandomBoxPacker<TPartGen,TCubicPackWrap>::getRandomPoint() const
    {
      return 
        Vec3(
          getRandom(this->getBBox().getMinPt().X(), this->getBBox().getMaxPt().X()),
          getRandom(this->getBBox().getMinPt().Y(), this->getBBox().getMaxPt().Y()),
          getRandom(this->getBBox().getMinPt().Z(), this->getBBox().getMaxPt().Z())
        );
    }

    template <typename TPartGen,template <typename TTPartGen> class TCubicPackWrap>
    typename RandomBoxPacker<TPartGen,TCubicPackWrap>::ParticleVector
    RandomBoxPacker<TPartGen,TCubicPackWrap>::getClosestNeighbours(
      const Particle& particle,
      int numClosest
    )
    {
      ParticleVector neighbourVector =
        this->getNTable().getUniqueNeighbourVector(
          particle.getPos(),
          this->getTolerance()
        );

      if (static_cast<int>(neighbourVector.size()) < numClosest)
      {
        neighbourVector =
          this->getNTable().getUniqueNeighbourVector(
            particle.getPos(),
            particle.getRad() + this->getTolerance()
          );
        if (static_cast<int>(neighbourVector.size()) < numClosest)
        {
          neighbourVector =
            this->getNTable().getUniqueNeighbourVector(
              particle.getPos(),
              this->getParticleGenerator().getMaxFitRadius() + this->getTolerance()
            );
        }
      }

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
    bool RandomBoxPacker<TPartGen,TCubicPackWrap >::particleIsValid(
      const Particle &particle
    ) const
    {
      return 
        (
          this->getParticleGenerator().isValidFitRadius(particle.getRad())
          &&
          Inherited::particleFitsInBBoxWithNeighbours(particle)
        );
    }

    template <typename TPartGen,template <typename TTPartGen> class TCubicPackWrap>
    int RandomBoxPacker<TPartGen,TCubicPackWrap>::getMaxInsertionFailures() const
    {
      return m_maxInsertionFailures;
    }

    template <typename TPartGen,template <typename TTPartGen> class TCubicPackWrap>
    void RandomBoxPacker<TPartGen,TCubicPackWrap>::generateRandomFill()
    {
      StuffedParticleIterator it =
        StuffedParticleIterator(
          *this,
          getMaxInsertionFailures(),
          getFitPlaneVector()
        );
      while (it.hasNext())
      {
        const Particle p = it.next();
        this->createAndInsertParticle(p);
      }
      it.logInfo();
    }

    template <typename TPartGen,template <typename TTPartGen> class TCubicPackWrap>
    void RandomBoxPacker<TPartGen,TCubicPackWrap>::generate()
    {
      this->generateCubicPacking();
      this->generateRandomFill();
    }

  }
}
