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


#ifndef ESYS_LSMRANDOMBLOCKGENERATOR_H
#define ESYS_LSMRANDOMBLOCKGENERATOR_H

#include <Geometry/BlockGenerator.h>
#include <Geometry/Plane3D.h>

#include <vector>
#include <boost/shared_ptr.hpp>

namespace esys
{
  namespace lsm
  {
    class ParticleFitter;

    typedef std::vector<Plane3D>                PlaneVector;
    typedef boost::shared_ptr<ParticleFitter> FitterPtr;
    typedef std::vector<FitterPtr>            FitterPtrVector;

    /**
     *
     */
    class RandomBlockGenerator : public BlockGenerator
    {
    public:
      RandomBlockGenerator(
        NTable            &nTable,
        ParticlePool      &particlePool,
        const BoundingBox &bBox,
        const BoolVector  &periodicDimensions,
        double            tolerance,
        double            minSphereRadius,
        double            maxSphereRadius,
        const PlaneVector &fitPlaneVector,
        int               maxInsertionFailures
      );

      virtual ~RandomBlockGenerator();

      virtual bool particleFits(const SimpleParticle &particle) const;

      virtual void generate();

      double getRandom(double min, double max) const;

      virtual double getRadius() const;

      virtual double getGridRadius() const;

      Vec3 getRandomPoint() const;

      ParticleVector getClosestNeighbors(const SimpleParticle& particle, int numClosest);

      int getMaxInsertionFailures() const;

      FitterPtrVector getFitterPtrVector();

      void generateFillParticles();

      const PlaneVector &getFitPlaneVector() const;

      Plane3D getClosestFitPlane(const SimpleParticle &particle) const;

    private:
      double      m_minRadius;
      double      m_maxRadius;
      PlaneVector m_fitPlaneVector;
      int         m_maxInsertionFailures;
    };
  }
}

#endif
