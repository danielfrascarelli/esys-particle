/////////////////////////////////////////////////////////////
//                                                         //
// Copyright (c) 2003-2014 by The University of Queensland //
// Centre for Geoscience Computing                         //
// http://earth.uq.edu.au/centre-geoscience-computing      //
//                                                         //
// Primary Business: Brisbane, Queensland, Australia       //
// Licensed under the Open Software License version 3.0    //
// http://www.apache.org/licenses/LICENSE-2.0          //
//                                                         //
/////////////////////////////////////////////////////////////


#ifndef __FLUID_INTERACTION_H
#define __FLUID_INTERACTION_H


#include "Model/FluidCell.h"
#include "Model/Particle.h"
#include "Foundation/vec3.h"
#include "Foundation/triplet.h"

/*!
  \class CFluidInteraction
  \brief base class for fluid interactions

  \author Qi Shao
  $Revision$
  $Date$
*/

class CFluidInteraction {
 protected:
  CFluidCell *m_cell;
  CParticle *m_p;
  Vec3 m_drag;
  Vec3 m_buoyancy;
  Vec3 m_force;

 public:
  CFluidInteraction (CParticle*,CFluidCell*);
  virtual ~CFluidInteraction(){};

  void calcForces();

  inline Vec3 getParticlePos() const {return m_p->getPos();};
  inline Vec3 getDrag() const {return m_drag;};
  inline Vec3 getBuoyancy() const {return m_buoyancy;};
  inline Vec3 getForce() const {return m_force;};
  inline double getVbsDrag() const {return m_drag.norm();};
  inline double getVbsBuoyancy() const {return m_buoyancy.norm();};
  inline double getVbsForce() const {return m_force.norm();};
  inline double getParticleID() const {return m_p->getID();};
  inline esys::lsm::triplet<int,Vec3,double> getParticleData() const
  {
    return
      esys::lsm::triplet<int,Vec3,double>(
        m_p->getID(),
        m_p->getPos(),
        m_p->getRad()
      );
  };

  typedef double (CFluidInteraction::* ScalarFieldFunction)() const;
  typedef Vec3 (CFluidInteraction::* VectorFieldFunction)() const;

  static ScalarFieldFunction getScalarFieldFunction(const string&);
  static VectorFieldFunction getVectorFieldFunction(const string&);

};

#endif //__FLUID_INTERACTION_H
