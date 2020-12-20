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

#ifndef ESYS_LSMPARTICLEFITTER_H
#define ESYS_LSMPARTICLEFITTER_H

#include "Geometry/RandomBlockGenerator.h"
#include "Geometry/Sphere3d.h"
#include "Geometry/Sphere2d.h"

namespace esys
{
  namespace lsm
  {
    class ParticleFitter
    {
    public:
      typedef RandomBlockGenerator::ParticleVector ParticleVector;

      ParticleFitter(RandomBlockGenerator &blockGenerator)
        : m_pGenerator(&blockGenerator),
          m_successfulFitCount(0),
          m_getFitCount(0),
          m_failedFitCount(0)
      {
      }

      virtual ~ParticleFitter() {}

      virtual SimpleParticle getFitParticle(
        const SimpleParticle &particle,
        const ParticleVector &neighbours,
        const Plane3D &plane
      ) = 0;
      
      static const SimpleParticle INVALID;
      
      void incrGetFit()
      {
        m_getFitCount++;
      }

      void incrFailedFit()
      {
        m_failedFitCount++;
      }

      void incrSuccessfulFit()
      {
        m_successfulFitCount++;
      }

      virtual std::string getName() const = 0;
      
      void write(std::ostream &oStream) const
      {
        oStream
          << (getName() + ": ")
          << "fit count = " << m_getFitCount
          << ", success = " << m_successfulFitCount
          << ", fail    = " << m_failedFitCount;
      }
      
      std::string toString() const
      {
        std::stringstream sStream;
        write(sStream);
        return sStream.str();
      }

      virtual bool particleFits(const SimpleParticle &particle) const
      {
        return getGenerator().particleFits(particle);
      }

    protected:
      RandomBlockGenerator &getGenerator()
      {
        return *m_pGenerator;
      }

      const RandomBlockGenerator &getGenerator() const
      {
        return *m_pGenerator;
      }
    private:
      RandomBlockGenerator *m_pGenerator;
      int m_successfulFitCount;
      int m_getFitCount;
      int m_failedFitCount;
    };
    
    inline std::ostream &operator<<(std::ostream &oStream, const ParticleFitter &fitter)
    {
      fitter.write(oStream);
      return oStream;
    }

    class MoveToSurfaceFitter : public ParticleFitter
    {
    public:
      MoveToSurfaceFitter(RandomBlockGenerator &blockGenerator)
        : ParticleFitter(blockGenerator)
      {
      }
      
      virtual std::string getName() const
      {
        return "Mv to Surface";
      }
      
      SimpleParticle moveToSurface(
        const SimpleParticle &stationary,
        const SimpleParticle &movable
      )
      {
        SimpleParticle moved = movable;
        const Vec3 centreDiff = movable.getPos() - stationary.getPos();
        const double centreDist = centreDiff.norm();
        if (centreDist > 0.0) {
          const Vec3 newCentrePos = 
            stationary.getPos()
            +
            (centreDiff/(centreDist))*(stationary.getRad() + movable.getRad());
          moved.moveTo(newCentrePos);
        }
        return moved;
      }
      
      virtual SimpleParticle getFitParticle(
        const SimpleParticle &particle,
        const ParticleVector &neighbours,
        const Plane3D &plane
      )
      {
        incrGetFit();
        SimpleParticle newParticle = particle;
        if (neighbours.size() > 0) {
          if (
            (particle.getPos() - neighbours[0]->getPos()).norm()
            <
            (particle.getRad() + neighbours[0]->getRad())
          ) {
            newParticle = moveToSurface(*(neighbours[0]), particle);
          }
        }
        if (newParticle.isValid() && !particleFits(newParticle)) {
          newParticle = INVALID;
          incrFailedFit();
        } else if (newParticle.isValid()) {
          incrSuccessfulFit();
        }
        return newParticle;
      }
    };

    class ThreeDParticleFitter : public ParticleFitter
    {
    public:
      ThreeDParticleFitter(RandomBlockGenerator &blockGenerator)
        : ParticleFitter(blockGenerator)
      {
      }
      
      virtual std::string getName() const
      {
        return "    3D Fitter";
      }
      
      SimpleParticle findAFit(
        const SimpleParticle& Po,
        const ParticleVector &particleVector
      )
      {
        SimpleParticle fittedParticle = SimpleParticle::INVALID;
        Vec3 M;
        double r;
        int id=Po.getID();
        
        if (particleVector.size() > 3) {
          const bool foundFit = 
            Sphere3D::FillIn(
              particleVector[0]->getPos(),
              particleVector[1]->getPos(),
              particleVector[2]->getPos(),
              particleVector[3]->getPos(),
              particleVector[0]->getRad(),
              particleVector[1]->getRad(),
              particleVector[2]->getRad(),
              particleVector[3]->getRad(),
              M,
              r
            );
          if (foundFit) {
            fittedParticle = SimpleParticle(M, r, id);
          }  
        }  else {
          throw std::runtime_error("findAFit: particleVector argument has fewer than 4 elements.");
        }
        return fittedParticle;
      }      
      
      virtual SimpleParticle getFitParticle(
        const SimpleParticle &particle,
        const ParticleVector &neighbours,
        const Plane3D &plane
      )
      {
        incrGetFit();
        SimpleParticle newParticle = INVALID;
        if(neighbours.size() > 3){ // at least 4 neighbors 
          const SimpleParticle &Pi = *(neighbours[0]);
          const double ndist=(particle.getPos()-Pi.getPos()).norm();
          if (ndist > 0.0) {
            newParticle = particle;
            if (ndist < Pi.getRad()){ // if Po inside Pi -> move Po to the surface of Pi
              newParticle.moveTo(
                Pi.getPos()+((particle.getPos()-Pi.getPos())*(Pi.getRad()/ndist))
              );
            }
            const double dist_p = plane.sep(newParticle.getPos());
            const double dist_3 = (newParticle.getPos()-neighbours[3]->getPos()).norm()-neighbours[3]->getRad();
            if (dist_p > dist_3) {  // 4th particle closer than plane -> fit 4 particles
              newParticle = findAFit(newParticle, neighbours);
            }
            else {
              newParticle = INVALID;
            }
          }
        }
        if (newParticle.isValid() && !particleFits(newParticle)) {
          newParticle = INVALID;
          incrFailedFit();
        } else if (newParticle.isValid()) {
          incrSuccessfulFit();
        }
        return newParticle;
      }      
    };

    class TwoDParticleFitter : public ParticleFitter
    {
    public:
      TwoDParticleFitter(RandomBlockGenerator &blockGenerator)
        : ParticleFitter(blockGenerator)
      {
      }
      
      virtual std::string getName() const
      {
        return "    2D Fitter";
      }
      
      SimpleParticle findAFit(
        const SimpleParticle& Po,
        const ParticleVector &particleVector
      )
      {
        SimpleParticle fittedParticle = SimpleParticle::INVALID;
        Vec3 M;
        double r;
        int id=Po.getID();
        
        if (particleVector.size() > 2) {
          const bool foundFit = 
            Sphere2D::FillIn(
              particleVector[0]->getPos(),
              particleVector[1]->getPos(),
              particleVector[2]->getPos(),
              particleVector[0]->getRad(),
              particleVector[1]->getRad(),
              particleVector[2]->getRad(),
              M,
              r
            );
          if (foundFit) {
            fittedParticle = SimpleParticle(M, r, id);
          }  
        }  else {
          throw std::runtime_error("findAFit: particleVector argument has fewer than 3 elements.");
        }
        return fittedParticle;
      }

      virtual SimpleParticle getFitParticle(
        const SimpleParticle &particle,
        const ParticleVector &neighbours,
        const Plane3D &plane
      )
      {
        incrGetFit();
        SimpleParticle newParticle = INVALID;
        if(neighbours.size() > 2){ // at least 3 neighbors
          const SimpleParticle &Pi = *(neighbours[0]);
          const double ndist=(particle.getPos()-Pi.getPos()).norm();
          if (ndist > 0.0) {
            newParticle = particle;
            if (ndist < Pi.getRad()){ // if Po inside Pi -> move Po to the surface of Pi
              newParticle.moveTo(
                Pi.getPos()+((particle.getPos()-Pi.getPos())*(Pi.getRad()/ndist))
              );
            }
            const double dist_p = plane.sep(newParticle.getPos());
            const double dist_2 = (newParticle.getPos()-neighbours[2]->getPos()).norm()-neighbours[2]->getRad();
            if (dist_p > dist_2) {  // 4th particle closer than plane -> fit 4 particles
              newParticle = findAFit(newParticle, neighbours);
            }
            else {
              newParticle = INVALID;
            }
          }
        }
        if (newParticle.isValid() && !particleFits(newParticle)) {
          newParticle = INVALID;
          incrFailedFit();
        } else if (newParticle.isValid()) {
          incrSuccessfulFit();
        }
        return newParticle;
      }      
    };

    class TwoDPlaneParticleFitter : public ParticleFitter
    {
    public:
      TwoDPlaneParticleFitter(RandomBlockGenerator &blockGenerator)
        : ParticleFitter(blockGenerator)
      {
      }

      virtual std::string getName() const
      {
        return "     2D Plane";
      }
      
      SimpleParticle findAFit(
        const SimpleParticle& particle,
        const ParticleVector &particleVector,
        const Plane3D& plane
      )
      {
        SimpleParticle fittedParticle = SimpleParticle::INVALID;
        Vec3 M;
        double r;
        const int id = particle.getID();
  
        if (particleVector.size() > 1) {
          const bool foundFit =
            Sphere2D::FillInWP(
              particleVector[0]->getPos(),
              particleVector[1]->getPos(),
              plane.GetO(),
              plane.GetW(),
              particleVector[0]->getRad(),
              particleVector[1]->getRad(),
              M,
              r
          );
          if (foundFit) {
            fittedParticle = SimpleParticle(M, r, id);
          }
        } else {
          throw std::runtime_error("findAFit: particleVector vector contains less than 2 particles.");
        }
  
        return fittedParticle;   
      }

      virtual SimpleParticle getFitParticle(
        const SimpleParticle &particle,
        const ParticleVector &neighbours,
        const Plane3D &plane
      )
      {
        incrGetFit();
        SimpleParticle newParticle = INVALID;
        if (neighbours.size() >= 2) { // 2 neighbors  -> try  2 particles + plane
          const SimpleParticle &Pi = *(neighbours[0]);
          const double ndist=(particle.getPos()-Pi.getPos()).norm();
          if (
            (ndist > 0.0)
            && 
            (
              (neighbours.size() == 2)
              ||
              (
                plane.sep(particle.getPos())
                <
                (particle.getPos()-neighbours[2]->getPos()).norm()-neighbours[2]->getRad()
              )
            )
          ) {
            newParticle = particle;
            if (ndist < Pi.getRad()) { // if Po inside Pi -> move Po to the surface of Pi
              newParticle.moveTo(
                Pi.getPos() + ((particle.getPos() - Pi.getPos()) * (Pi.getRad()/ndist))
              );
            }
            newParticle = findAFit(newParticle, neighbours, plane);
          }
        }
        if (newParticle.isValid() && !particleFits(newParticle)) {
          newParticle = INVALID;
          incrFailedFit();
        } else if (newParticle.isValid()) {
          incrSuccessfulFit();
        }
        return newParticle;
      }
    };

    class ThreeDPlaneParticleFitter : public ParticleFitter
    {
    public:
      ThreeDPlaneParticleFitter(RandomBlockGenerator &blockGenerator)
        : ParticleFitter(blockGenerator)
      {
      }

      virtual std::string getName() const
      {
        return "     3D Plane";
      }
      
      SimpleParticle findAFit(
        const SimpleParticle& particle,
        const ParticleVector &particleVector,
        const Plane3D& plane
      )
      {
        SimpleParticle fittedParticle = SimpleParticle::INVALID;
        Vec3 M;
        double r;
        const int id = particle.getID();
  
        if (particleVector.size() > 2) {
          const bool foundFit =
            Sphere3D::FillInWP(
              particleVector[0]->getPos(),
              particleVector[1]->getPos(),
              particleVector[2]->getPos(),
              plane.GetO(),
              plane.GetW(),
              particleVector[0]->getRad(),
              particleVector[1]->getRad(),
              particleVector[2]->getRad(),
              M,
              r
          );
          if (foundFit) {          
            fittedParticle = SimpleParticle(M, r, id);
          }
        } else {
          throw std::runtime_error("findAFit: particleVector vector contains less than 3 particles.");
        }
  
        return fittedParticle;   
      }

      virtual SimpleParticle getFitParticle(
        const SimpleParticle &particle,
        const ParticleVector &neighbours,
        const Plane3D &plane
      )
      {
        incrGetFit();
        SimpleParticle newParticle = INVALID;
        if (neighbours.size() >= 3) { // 3 neighbors  -> try  3 particles + plane
          const SimpleParticle &Pi = *(neighbours[0]);
          const double ndist=(particle.getPos()-Pi.getPos()).norm();
          if (
            (ndist > 0.0)
            && 
            (
              (neighbours.size() == 3)
              || 
              (
                plane.sep(particle.getPos())
                <
                (particle.getPos()-neighbours[3]->getPos()).norm()-neighbours[3]->getRad()
              )
            )
          ) {
            newParticle = particle;
            if (ndist < Pi.getRad()) { // if Po inside Pi -> move Po to the surface of Pi
              newParticle.moveTo(
                Pi.getPos() + ((particle.getPos() - Pi.getPos()) * (Pi.getRad()/ndist))
              );
            }
            newParticle = findAFit(newParticle, neighbours, plane);
          }
        }
        if (newParticle.isValid() && !particleFits(newParticle)) {
          newParticle = INVALID;
          incrFailedFit();
        } else if (newParticle.isValid()) {
          incrSuccessfulFit();
        }
        return newParticle;
      }
    };
  }
}

#endif
