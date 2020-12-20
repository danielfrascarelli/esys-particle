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
#include "Geometry/RandomBlockGenerator.h"
#include "Geometry/SimpleParticle.h"
#include "Geometry/GridIterator.h"
#include "Geometry/SimpleNTable.h"
#include "Geometry/Sphere3d.h"
#include "Geometry/Plane3D.h"
#include "Geometry/ParticleFitter.h"

#include <algorithm>
#include <stdexcept>
#include <float.h>

namespace esys
{
  namespace lsm
  {
    RandomBlockGenerator::RandomBlockGenerator(
      NTable            &nTable,
      ParticlePool      &particlePool,
      const BoundingBox &bBox,
      const BoolVector  &periodicDimensions,
      double            tolerance,
      double            minSphereRadius,
      double            maxSphereRadius,
      const PlaneVector &fitPlaneVector,
      int               maxInsertionFailures
    )
      : BlockGenerator(nTable, particlePool, bBox, periodicDimensions, tolerance),
        m_minRadius(minSphereRadius),
        m_maxRadius(maxSphereRadius),
        m_fitPlaneVector(fitPlaneVector),
        m_maxInsertionFailures(maxInsertionFailures)
    {
    }

    RandomBlockGenerator::~RandomBlockGenerator()
    {
    }

    double RandomBlockGenerator::getRandom(double minVal, double maxVal) const
    {
      return 
        minVal
        +
        ((maxVal-minVal)*(static_cast<double>(rand()))/(static_cast<double>(RAND_MAX)));
    }

    double RandomBlockGenerator::getGridRadius() const
    {
      return m_maxRadius;
    }

    double RandomBlockGenerator::getRadius() const
    {
      return getRandom(m_minRadius, m_maxRadius);
    }

    class PlaneComparer
    {
    public:
      PlaneComparer(const SimpleParticle &particle)
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
      const SimpleParticle *m_pParticle;
    };

    const PlaneVector &RandomBlockGenerator::getFitPlaneVector() const
    {
      return m_fitPlaneVector;
    }

    Plane3D RandomBlockGenerator::getClosestFitPlane(const SimpleParticle &particle) const
    {
      PlaneVector fitPlanes = getFitPlaneVector();
      if (fitPlanes.size() > 0) {
        std::sort(fitPlanes.begin(), fitPlanes.end(), PlaneComparer(particle));

        return fitPlanes[0];
      }
      return Plane3D(Vec3(-1.0, 0, 0), Vec3(DBL_MAX/2, DBL_MAX/2, DBL_MAX/2));
    }

    Vec3 RandomBlockGenerator::getRandomPoint() const
    {
      return 
        Vec3(
          getRandom(getBBox().getMinPt().X(), getBBox().getMaxPt().X()),
          getRandom(getBBox().getMinPt().Y(), getBBox().getMaxPt().Y()),
          getRandom(getBBox().getMinPt().Z(), getBBox().getMaxPt().Z())
        );
    }

    RandomBlockGenerator::ParticleVector RandomBlockGenerator::getClosestNeighbors(
      const SimpleParticle& particle,
      int numClosest
    )
    {
      ParticleVector neighbourVector =
        getNTable().getUniqueNeighbourVector(particle.getPos(), m_maxRadius + getTolerance());
      std::sort(neighbourVector.begin(), neighbourVector.end(), ParticleComparer(particle));

      if (neighbourVector.size() > static_cast<size_t>(numClosest)) {
        neighbourVector.erase(neighbourVector.begin() + numClosest, neighbourVector.end());
      }

      return neighbourVector;
    }

    bool RandomBlockGenerator::particleFits(const SimpleParticle &particle) const
    {
      return 
        (
          (particle.getRad() >= m_minRadius)
          &&
          (particle.getRad() <= m_maxRadius)
          &&
          BlockGenerator::particleFits(particle)
        );
    }

    int RandomBlockGenerator::getMaxInsertionFailures() const
    {
      return m_maxInsertionFailures;
    }

    FitterPtrVector RandomBlockGenerator::getFitterPtrVector()
    {
      FitterPtrVector fitters;
      fitters.push_back(FitterPtr(new MoveToSurfaceFitter(*this)));
      if (is2d())
      {
        fitters.push_back(FitterPtr(new TwoDParticleFitter(*this)));
        if (getFitPlaneVector().size() > 0) {
          fitters.push_back(FitterPtr(new TwoDPlaneParticleFitter(*this)));
        }
      }
      else
      {
        fitters.push_back(FitterPtr(new ThreeDParticleFitter(*this)));
        if (getFitPlaneVector().size() > 0) {
          fitters.push_back(FitterPtr(new ThreeDPlaneParticleFitter(*this)));
        }
      }

      return fitters;
    }

    void RandomBlockGenerator::generateFillParticles()
    {
      int numFails    = 0;
      int insertCount = 0;
      FitterPtrVector fitters = getFitterPtrVector();
      SimpleParticle particle = ParticleFitter::INVALID;
      SimpleParticle fitParticle = ParticleFitter::INVALID;
      Plane3D plane;
      ParticleVector neighbours;
      while (numFails < getMaxInsertionFailures()) {
        fitParticle = ParticleFitter::INVALID;
        particle = generateParticle(getRandomPoint());
        neighbours = getClosestNeighbors(particle, 4);
        plane = getClosestFitPlane(particle);
        
        for (
          FitterPtrVector::iterator it = fitters.begin();
          (
            (it != fitters.end())
            &&
            (!fitParticle.isValid())
          );
          it++
        )
        {
          fitParticle = (*it)->getFitParticle(particle, neighbours, plane);
        }
        if (fitParticle.isValid()) {
          numFails = 0;
          insertCount++;
          insertParticle(fitParticle);
        } else {
          numFails++;
        }
      }
      console.Info()
        << "BoundingBox: minPt = " << getBBox().getMinPt()
        << ", maxPt = " << getBBox().getMaxPt() << "\n"
        << "BoundingBox: sizes = " << getBBox().getSizes() << "\n";
      for (
          FitterPtrVector::iterator it = fitters.begin();
          (it != fitters.end());
          it++      
      ) {
        console.Info() << (*(*it)).toString() << "\n";
      }
      console.Info() << "Total inserted = " << insertCount << "\n";
    }

    void RandomBlockGenerator::generate()
    {
      generateSeedParticles();
      generateFillParticles();
    }
  }
}
