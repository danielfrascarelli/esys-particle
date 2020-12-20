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


#ifndef ESYS_LSMSPHEREFITTER_H
#define ESYS_LSMSPHEREFITTER_H

#include "Geometry/Sphere3d.h"
#include "Geometry/Sphere2d.h"
#include "Foundation/BoundingSphere.h"

#include <stdexcept>
#include <sstream>

namespace esys
{
  namespace lsm
  {
    template <typename TmplFitTraits>
    class SphereFitter
    {
    public:
      typedef TmplFitTraits                      FitTraits;
      typedef typename FitTraits::Validator      Validator;
      typedef typename FitTraits::Particle       Particle;
      typedef typename Validator::ParticleVector ParticleVector;
      typedef typename FitTraits::Plane3D          Plane3D;

      SphereFitter(const std::string &name, Validator &validator)
        : m_pValidator(&validator),
          m_successfulFitCount(0),
          m_getFitCount(0),
          m_failedFitCount(0),
          m_name(name)
      {
      }

      virtual ~SphereFitter()
      {
      }

      virtual Particle getFitParticle(
        const Particle &particle,
        const ParticleVector &neighbours,
        const Plane3D &plane
      ) = 0;

      static Particle getInvalidParticle()
      {
        return Particle::INVALID;
      }

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

      const std::string &getName() const
      {
        return m_name;
      }
      
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

      bool particleIsValid(const Particle &particle) const
      {
        return getValidator().particleIsValid(particle);
      }

    protected:
      Validator &getValidator()
      {
        return *m_pValidator;
      }

      const Validator &getValidator() const
      {
        return *m_pValidator;
      }
    private:
      Validator   *m_pValidator;
      int         m_successfulFitCount;
      int         m_getFitCount;
      int         m_failedFitCount;
      std::string m_name;
    };

    template <typename TmplFitTraits>    
    inline std::ostream &operator<<(
      std::ostream &oStream,
      const SphereFitter<TmplFitTraits> &fitter
    )
    {
      fitter.write(oStream);
      return oStream;
    }

    template <typename TmplFitTraits>
    class MoveToSurfaceFitter : public SphereFitter<TmplFitTraits>
    {
    public:
      typedef SphereFitter<TmplFitTraits>        Inherited;
      typedef typename Inherited::Validator      Validator;
      typedef typename Inherited::Particle       Particle;
      typedef typename Inherited::ParticleVector ParticleVector;
      typedef typename Inherited::Plane3D          Plane3D;
      MoveToSurfaceFitter(Validator &validator)
        : Inherited("Mv to Surface", validator)
      {
      }
      
      Particle moveToSurface(
        const Particle &stationary,
        const Particle &movable
      )
      {
        Particle moved = movable;
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
      
      virtual Particle getFitParticle(
        const Particle &particle,
        const ParticleVector &neighbours,
        const Plane3D &plane
      )
      {
        this->incrGetFit();
        Particle newParticle = particle;
        if (neighbours.size() > 0) {
          if (
            (particle.getPos() - neighbours[0]->getPos()).norm()
            <
            (particle.getRad() + neighbours[0]->getRad())
          ) {
            newParticle = moveToSurface(*(neighbours[0]), particle);
          }
        }
        if (newParticle.isValid() && !this->particleIsValid(newParticle)) {
          newParticle = this->getInvalidParticle();
          this->incrFailedFit();
        } else if (newParticle.isValid()) {
          this->incrSuccessfulFit();
        }
        return newParticle;
      }
    };

    template <typename TmplFitTraits>
    class ThreeDSphereFitter : public SphereFitter<TmplFitTraits>
    {
    public:
      typedef SphereFitter<TmplFitTraits>        Inherited;
      typedef typename Inherited::Validator      Validator;
      typedef typename Inherited::Particle       Particle;
      typedef typename Inherited::ParticleVector ParticleVector;
      typedef typename Inherited::Plane3D          Plane3D;

      ThreeDSphereFitter(Validator &validator)
        : Inherited("    3D Fitter", validator)
      {
      }

      Particle findAFit(
        const Particle& Po,
        const ParticleVector &particleVector
      )
      {
        Particle fittedParticle = this->getInvalidParticle();
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
            fittedParticle = Particle(M, r, id);
          }  
        }  else {
          throw std::runtime_error("findAFit: particleVector argument has fewer than 4 elements.");
        }
        return fittedParticle;
      }      
      
      virtual Particle getFitParticle(
        const Particle &particle,
        const ParticleVector &neighbours,
        const Plane3D &plane
      )
      {
        this->incrGetFit();
        Particle newParticle = this->getInvalidParticle();
        if(neighbours.size() > 3){ // at least 4 neighbors 
          const Particle &Pi = *(neighbours[0]);
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
              newParticle = this->getInvalidParticle();
            }
          }
        }
        if (newParticle.isValid() && !this->particleIsValid(newParticle)) {
          newParticle = this->getInvalidParticle();
          this->incrFailedFit();
        } else if (newParticle.isValid()) {
          this->incrSuccessfulFit();
        }
        return newParticle;
      }      
    };

    template <typename TmplFitTraits>
    class TwoDSphereFitter : public SphereFitter<TmplFitTraits>
    {
    public:
      typedef SphereFitter<TmplFitTraits>        Inherited;
      typedef typename Inherited::Validator      Validator;
      typedef typename Inherited::Particle       Particle;
      typedef typename Inherited::ParticleVector ParticleVector;
      typedef typename Inherited::Plane3D          Plane3D;

      TwoDSphereFitter(Validator &validator)
        : Inherited("    2D Fitter", validator)
      {
      }
      
      Particle findAFit(
        const Particle& Po,
        const ParticleVector &particleVector
      )
      {
        Particle fittedParticle = this->getInvalidParticle();
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
            fittedParticle = Particle(M, r, id);
          }  
        }  else {
          throw std::runtime_error("findAFit: particleVector argument has fewer than 3 elements.");
        }
        return fittedParticle;
      }

      virtual Particle getFitParticle(
        const Particle &particle,
        const ParticleVector &neighbours,
        const Plane3D &plane
      )
      {
        this->incrGetFit();
        Particle newParticle = this->getInvalidParticle();
        if(neighbours.size() > 2){ // at least 3 neighbors
          const Particle &Pi = *(neighbours[0]);
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
            if (dist_p > dist_2) {  // 3rd particle closer than plane -> fit 3 particles
              newParticle = findAFit(newParticle, neighbours);
            }
            else {
              newParticle = this->getInvalidParticle();
            }
          }
        }
        if (newParticle.isValid() && !this->particleIsValid(newParticle)) {
          newParticle = this->getInvalidParticle();
          this->incrFailedFit();
        } else if (newParticle.isValid()) {
          this->incrSuccessfulFit();
        }
        return newParticle;
      }      
    };

    template <typename TmplFitTraits>
    class TwoDSphereSphereFitter : public SphereFitter<TmplFitTraits>
    {
    public:
      typedef SphereFitter<TmplFitTraits>        Inherited;
      typedef typename Inherited::Validator      Validator;
      typedef typename Inherited::Particle       Particle;
      typedef typename Inherited::ParticleVector ParticleVector;
      typedef typename Inherited::Plane3D          Plane3D;

      TwoDSphereSphereFitter(
        Validator &validator,
        const BoundingSphere &sphere
      )
        : Inherited(" 2D Sph Fitter", validator),
          m_sphere(sphere)
      {
      }
      
      Particle findAFit(
        const Particle& Po,
        const ParticleVector &particleVector
      )
      {
        Particle fittedParticle = this->getInvalidParticle();
        Vec3 M;
        double r;
        int id=Po.getID();
        
        if (particleVector.size() >= 2) {
          const bool foundFit = 
            Sphere2D::FillIn(
              m_sphere.getCentre(),
              particleVector[0]->getPos(),
              particleVector[1]->getPos(),
              -m_sphere.getRadius(),
              particleVector[0]->getRad(),
              particleVector[1]->getRad(),
              M,
              r
            );
          if (foundFit) {
            fittedParticle = Particle(M, r, id);
          }
        }  else {
          throw std::runtime_error(
            "TwoDSphereSphereFitter::findAFit: particleVector "
            "argument has fewer than 2 elements."
          );
        }
        return fittedParticle;
      }

      virtual Particle getFitParticle(
        const Particle &particle,
        const ParticleVector &neighbours,
        const Plane3D &plane
      )
      {
        this->incrGetFit();
        Particle newParticle = this->getInvalidParticle();
        if (neighbours.size() >= 2) { // 2 neighbours: try  2 particles + bounding-sphere
          const Particle &Pi = *(neighbours[0]);
          const double ndist=(particle.getPos()-Pi.getPos()).norm();
          if (
            (ndist > 0.0)
            && 
            (
              (neighbours.size() == 2)
              ||
              (
                fabs((m_sphere.getCentre()-particle.getPos()).norm()-m_sphere.getRadius())
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
            newParticle = findAFit(newParticle, neighbours);
          }
        }
        if (newParticle.isValid() && !this->particleIsValid(newParticle)) {
          newParticle = this->getInvalidParticle();
          this->incrFailedFit();
        } else if (newParticle.isValid()) {
          this->incrSuccessfulFit();
        }
        return newParticle;
      }

    private:
      BoundingSphere m_sphere;
    };

    template <typename TmplFitTraits>
    class ThreeDSphereSphereFitter : public SphereFitter<TmplFitTraits>
    {
    public:
      typedef SphereFitter<TmplFitTraits>        Inherited;
      typedef typename Inherited::Validator      Validator;
      typedef typename Inherited::Particle       Particle;
      typedef typename Inherited::ParticleVector ParticleVector;
      typedef typename Inherited::Plane3D          Plane3D;

      ThreeDSphereSphereFitter(
        Validator &validator,
        const BoundingSphere &sphere
      )
        : Inherited(" 3D Sph Fitter", validator),
          m_sphere(sphere)
      {
      }
      
      Particle findAFit(
        const Particle& Po,
        const ParticleVector &particleVector
      )
      {
        Particle fittedParticle = this->getInvalidParticle();
        Vec3 M;
        double r;
        int id=Po.getID();
        
        if (particleVector.size() > 2) {
          const bool foundFit = 
            Sphere3D::FillIn(
              m_sphere.getCentre(),
              particleVector[0]->getPos(),
              particleVector[1]->getPos(),
              particleVector[2]->getPos(),
              -m_sphere.getRadius(),
              particleVector[0]->getRad(),
              particleVector[1]->getRad(),
              particleVector[2]->getRad(),
              M,
              r
            );
          if (foundFit) {
            fittedParticle = Particle(M, r, id);
          }
        } else {
          throw
            std::runtime_error(
              "ThreeDSphereSphereFitter::findAFit: particleVector argument"
              " has fewer than 3 elements."
            );
        }
        return fittedParticle;
      }

      virtual Particle getFitParticle(
        const Particle &particle,
        const ParticleVector &neighbours,
        const Plane3D &plane
      )
      {
        this->incrGetFit();
        Particle newParticle = this->getInvalidParticle();
        if (neighbours.size() >= 3) { // 3 neighbours: try  3 particles + bounding-sphere
          const Particle &Pi = *(neighbours[0]);
          const double ndist=(particle.getPos()-Pi.getPos()).norm();
          if (
            (ndist > 0.0)
            && 
            (
              (neighbours.size() == 3)
              ||
              (
                fabs((m_sphere.getCentre()-particle.getPos()).norm()-m_sphere.getRadius())
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
            newParticle = findAFit(newParticle, neighbours);
          }
        }
        if (newParticle.isValid() && !this->particleIsValid(newParticle)) {
          newParticle = this->getInvalidParticle();
          this->incrFailedFit();
        } else if (newParticle.isValid()) {
          this->incrSuccessfulFit();
        }
        return newParticle;
      }

    private:
      BoundingSphere m_sphere;
    };

    template <typename TmplFitTraits>
    class TwoDPlaneSphereFitter : public SphereFitter<TmplFitTraits>
    {
    public:
      typedef SphereFitter<TmplFitTraits>        Inherited;
      typedef typename Inherited::Validator      Validator;
      typedef typename Inherited::Particle       Particle;
      typedef typename Inherited::ParticleVector ParticleVector;
      typedef typename Inherited::Plane3D          Plane3D;

      TwoDPlaneSphereFitter(Validator &validator)
        : Inherited("     2D Plane", validator)
      {
      }

      Particle findAFit(
        const Particle& particle,
        const ParticleVector &particleVector,
        const Plane3D& plane
      )
      {
        Particle fittedParticle = this->getInvalidParticle();
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
            fittedParticle = Particle(M, r, id);
          }
        } else {
          throw std::runtime_error("findAFit: particleVector vector contains less than 2 particles.");
        }
  
        return fittedParticle;   
      }

      virtual Particle getFitParticle(
        const Particle &particle,
        const ParticleVector &neighbours,
        const Plane3D &plane
      )
      {
        this->incrGetFit();
        Particle newParticle = this->getInvalidParticle();
        if (neighbours.size() >= 2) { // 2 neighbors  -> try  2 particles + plane
          const Particle &Pi = *(neighbours[0]);
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
        if (newParticle.isValid() && !this->particleIsValid(newParticle)) {
          newParticle = this->getInvalidParticle();
          this->incrFailedFit();
        } else if (newParticle.isValid()) {
          this->incrSuccessfulFit();
        }
        return newParticle;
      }
    };

    template <typename TmplFitTraits>
    class ThreeDPlaneSphereFitter : public SphereFitter<TmplFitTraits>
    {
    public:
      typedef SphereFitter<TmplFitTraits>        Inherited;
      typedef typename Inherited::Validator      Validator;
      typedef typename Inherited::Particle       Particle;
      typedef typename Inherited::ParticleVector ParticleVector;
      typedef typename Inherited::Plane3D          Plane3D;

      ThreeDPlaneSphereFitter(Validator &validator)
        : Inherited("     3D Plane", validator)
      {
      }

      Particle findAFit(
        const Particle& particle,
        const ParticleVector &particleVector,
        const Plane3D& plane
      )
      {
        Particle fittedParticle = this->getInvalidParticle();
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
            fittedParticle = Particle(M, r, id);
          }
        } else {
          throw std::runtime_error("findAFit: particleVector vector contains less than 3 particles.");
        }
  
        return fittedParticle;   
      }

      virtual Particle getFitParticle(
        const Particle &particle,
        const ParticleVector &neighbours,
        const Plane3D &plane
      )
      {
        this->incrGetFit();
        Particle newParticle = this->getInvalidParticle();
        if (neighbours.size() >= 3) { // 3 neighbors  -> try  3 particles + plane
          const Particle &Pi = *(neighbours[0]);
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
        if (newParticle.isValid() && !this->particleIsValid(newParticle)) {
          newParticle = this->getInvalidParticle();
          this->incrFailedFit();
        } else if (newParticle.isValid()) {
          this->incrSuccessfulFit();
        }
        return newParticle;
      }
    };
  }
}

#include "Geometry/SphereFitter.hpp"

#endif
