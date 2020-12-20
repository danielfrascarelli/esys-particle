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


#ifndef ESYS_LSMBODYFORCEGROUP_H
#define ESYS_LSMBODYFORCEGROUP_H

#include "Model/InteractionGroup.h"
#include "Foundation/vec3.h"

template <class TmplParticle> class ParallelParticleArray;

namespace esys
{
  namespace lsm
  {
    class BodyForceIGP : public AIGParam
    {
    public:

      BodyForceIGP();

      BodyForceIGP(const std::string &name, const Vec3 &acceleration);

      virtual ~BodyForceIGP();

      const Vec3 &getAcceleration() const;

      const std::string &getName() const;

      virtual void packInto(CVarMPIBuffer *pBuffer) const;

      static BodyForceIGP extract(CVarMPIBuffer *pBuffer);
      
      virtual std::string getTypeString() const {return "BodyForce";}

    protected:
      Vec3 m_acceleration;
    };

    class GravityIGP : public BodyForceIGP
    {
    public:

      GravityIGP() : BodyForceIGP()
      {
      }

      GravityIGP(const std::string &name, const Vec3 &acceleration) : BodyForceIGP(name, acceleration)
      {
      }

      virtual std::string getTypeString() const {return "Gravity";}

    private:
    };

    class BuoyancyIGP : public AIGParam
    {
    public:

      BuoyancyIGP()
      {
      }

      BuoyancyIGP(const std::string &name, const Vec3 &acceleration, const double &fluidDensity, const double &fluidHeight);

      virtual std::string getTypeString() const {return "Buoyancy";}

      virtual void packInto(CVarMPIBuffer *pBuffer) const;

      static BuoyancyIGP extract(CVarMPIBuffer *pBuffer);

      const double &getFluidDensity () const;

      const double &getFluidHeight () const;

      const Vec3 &getAcceleration() const;

    private:
      Vec3 m_acceleration;
      double m_fluidDensity, m_fluidHeight;
    };

    /**
     * Objects of this class apply a gravitational body-acceleration
     * to individual particles.
     */
    template <class TmplParticle>
    class BodyForceGroup : public AInteractionGroup<TmplParticle>
    {
    public:
      typedef ParallelParticleArray<TmplParticle> ParticleArray;
      typedef typename ParticleArray::ParticleListIterator ParticleIterator;

      BodyForceGroup(const BodyForceIGP &prms, ParticleArray &particleArray);

      ~BodyForceGroup();

      /**
       * Returns the force which would be applied to a
       * particle of the specified mass.
       *
       * @param mass A mass ("units" assumed to be same as
       *             the acceleration units).
       */
      Vec3 getForce(double mass) const;

      /**
       * Applies body force to the specified particle.
       *
       * @param particle Force applied to this particle using
       *                 a call to particle.applyForce(...).
       */
      void applyForce(TmplParticle &particle) const;

      virtual void Update(ParallelParticleArray<TmplParticle> *particleArray);

      /**
       * Null op, time step size not required.
       */
      virtual void setTimeStepSize(double dt)
      {
      }
      
      virtual void calcForces();

    private:
      Vec3          m_acceleration;
      ParticleArray *m_pParticleArray;
    };

    template <class TmplParticle>
    class BuoyancyForceGroup : public AInteractionGroup<TmplParticle>
    {
    public:
      typedef ParallelParticleArray<TmplParticle> ParticleArray;
      typedef typename ParticleArray::ParticleListIterator ParticleIterator;

      BuoyancyForceGroup(const BuoyancyIGP &prms, ParticleArray &particleArray);

      ~BuoyancyForceGroup();

      /**
       * Returns the force which would be applied to a
       * particle of the specified mass.
       *
       * @param mass A mass ("units" assumed to be same as
       *             the acceleration units).
       */
      Vec3 getForce(double volume) const;

      /**
       * Applies body force to the specified particle.
       *
       * @param particle Force applied to this particle using
       *                 a call to particle.applyForce(...).
       */
      void applyForce(TmplParticle &particle) const;

      virtual void Update(ParallelParticleArray<TmplParticle> *particleArray);

      /**
       * Null op, time step size not required.
       */
      virtual void setTimeStepSize(double dt)
      {
      }

      virtual void calcForces();

    private:
      Vec3          m_acceleration;
      ParticleArray *m_pParticleArray;
      double m_fluidDensity, m_fluidHeight;
    };
  }
}

#include "Model/BodyForceGroup.hpp"

#endif
