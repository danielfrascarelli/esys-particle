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

#ifndef __EWALLINTERACTION_H
#define __EWALLINTERACTION_H

#include "Model/WallInteraction.h"
#include "Model/Wall.h"
#include "Model/Particle.h"
#include "Model/RotParticle.h"

/*!
  \class CElasticWallInteraction
  \brief unbonded elastic interaction between a particle and a wall

  \author Steffen Abe
  $Revision$
  $Date$  
*/
template <class T>
class CElasticWallInteraction : public AWallInteraction<T>
{
protected:
  double m_k;//!< spring constant
public:
  CElasticWallInteraction();
  CElasticWallInteraction(T*,CWall*,double,bool);
  virtual ~CElasticWallInteraction(){};

  virtual void calcForces();
  virtual Vec3 getForce();
  virtual void setPP(const vector<T*>){};
  virtual double getStiffness();
};

#include "EWallInteraction.hpp"

#endif //__EWALLINTERACTION_H
