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


#include "ppa/src/pp_array.h"
#include <math.h>

namespace esys
{
  namespace lsm
  {
    template <class TmplParticle>
    BodyForceGroup<TmplParticle>::BodyForceGroup(
      const BodyForceIGP &prms,
      ParticleArray &particleArray
    )
      : m_acceleration(prms.getAcceleration()),
        m_pParticleArray(&particleArray)
    {
    }

    template <class TmplParticle>
    BodyForceGroup<TmplParticle>::~BodyForceGroup()
    {
    }

    template <class TmplParticle>
    void BodyForceGroup<TmplParticle>::Update(ParallelParticleArray<TmplParticle> *particleArray)
    {
    }

    template <class TmplParticle>
    Vec3 BodyForceGroup<TmplParticle>::getForce(double mass) const
    {
      return m_acceleration*mass;
    }

    template <class TmplParticle>
    void BodyForceGroup<TmplParticle>::applyForce(TmplParticle &particle) const
    {
      particle.applyForce(getForce(particle.getMass()), particle.getPos());
    }

    template <class TmplParticle>
    void BodyForceGroup<TmplParticle>::calcForces()
    {
      typename ParticleArray::ParticleListHandle plh = m_pParticleArray->getAllParticles();
      for (
        ParticleIterator it = plh->begin();
        it != plh->end();
        it++
      )
      {
        applyForce(*(*it));
      }
    }

    template <class TmplParticle>
    BuoyancyForceGroup<TmplParticle>::BuoyancyForceGroup(
      const BuoyancyIGP &prms,
      ParticleArray &particleArray
    )
      : m_acceleration(prms.getAcceleration()),
        m_pParticleArray(&particleArray) , 
        m_fluidDensity(prms.getFluidDensity()),
        m_fluidHeight(prms.getFluidHeight())
    {
    }

    template <class TmplParticle>
    BuoyancyForceGroup<TmplParticle>::~BuoyancyForceGroup()
    {
    }

    template <class TmplParticle>
    Vec3 BuoyancyForceGroup<TmplParticle>::getForce(double volume) const
    {
      Vec3 force;
      
      force = -1.0*m_acceleration*m_fluidDensity*volume;

      return force;
    }

    template <class TmplParticle>
    void BuoyancyForceGroup<TmplParticle>::applyForce(TmplParticle &particle) const
    {
      double volume, radius, height;
      double theta, d;

      radius = particle.getRad();
      height = -1.0*particle.getPos() * m_acceleration.unit();
      if(particle.getDo2dCalculations()){
         if (m_fluidHeight >= height + radius) {
	    // completely submerged
            volume = M_PI*radius*radius;
         }
         else if (m_fluidHeight >= height - radius) {
	    // partially submerged
            if (m_fluidHeight < height) {
               d = height - m_fluidHeight;
               theta = 2.0 * acos (d / radius);
               volume = 0.5*radius*radius*(theta - sin (theta));
            }
            else {
               d = m_fluidHeight - height;
               theta = 2.0 * acos (d / radius);
               volume = M_PI*radius*radius*(M_PI - 0.5*(theta - sin (theta)));
            }
         }
         else {
            // not submerged
            volume = 0.0;
         }
      }
      else {
         if (m_fluidHeight >= height + radius) {
            // completely submerged
            volume = 4.0*M_PI*radius*radius*radius/3.0;
         }
         else if (m_fluidHeight >= height - radius) {
            // partially submerged
            if (m_fluidHeight < height) {
               d = radius - (height - m_fluidHeight);
               volume = M_PI*d*d*(3.0*radius - d)/3.0;
            }
            else {
               d = radius - (m_fluidHeight - height);
               volume = M_PI*(4.0*radius*radius*radius - d*d*(3.0*radius - d))/3.0;
            }
         }
         else {
            // not submerged
            volume = 0.0;
         }
      }
      particle.applyForce(getForce(volume), particle.getPos());
    }

    template <class TmplParticle>
    void BuoyancyForceGroup<TmplParticle>::calcForces()
    {
      typename ParticleArray::ParticleListHandle plh = m_pParticleArray->getAllParticles();
      for (
        ParticleIterator it = plh->begin();
        it != plh->end();
        it++
      )
      {
        applyForce(*(*it));
      }
    }

    template <class TmplParticle>
    void BuoyancyForceGroup<TmplParticle>::Update(ParallelParticleArray<TmplParticle> *particleArray)
    {
    }
  }
}
