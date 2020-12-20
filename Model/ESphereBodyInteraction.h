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

#ifndef __ESPHEREBODYINTERACTION_H
#define __ESPHEREBODYINTERACTION_H

#include "Model/SphereBodyInteraction.h"
#include "Model/SphereBody.h"
#include "Model/Particle.h"
#include "Model/RotParticle.h"

/*!
  \class CElasticSphereBodyInteraction
  \brief unbonded elastic interaction between a particle and a wall

  \author Steffen Abe
  $Revision$
  $Date$  
*/
template <class T>
class CElasticSphereBodyInteraction : public ASphereBodyInteraction<T>
{
protected:
  double m_k;//!< spring constant
public:
  CElasticSphereBodyInteraction();
  CElasticSphereBodyInteraction(T*,CSphereBody*,double,bool);
  virtual ~CElasticSphereBodyInteraction(){};

  virtual void calcForces();
  virtual Vec3 getForce();
  virtual void setPP(const vector<T*>){};
  virtual double getStiffness();
};

#include "ESphereBodyInteraction.hpp"

#endif //__ESPHEREBODYINTERACTION_H
