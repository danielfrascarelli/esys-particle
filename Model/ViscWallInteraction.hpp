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


#ifndef MODEL_VISCWALLINTERACTION_HPP
#define MODEL_VISCWALLINTERACTION_HPP

#include "Model/ViscWallInteraction.h"

/*!
  constructor with parameters

  \param p pointer to the particle
  \param w pointer to the wall
  \param k spring constant for the elastic interaction
  \param nu viscosity
  \param iflag inner flag
*/
template <class T>
CViscWallInteraction<T>::CViscWallInteraction(T* p,CWall* w,double nu,bool iflag):
  AWallInteraction<T>(p,w,iflag)
{
  m_nu=nu;
}

/*!
  calculate and apply viscous forces
  F=visc*vol*dv 

  \warning Hack - currently using mass instead of volume
*/
template <class T>
void CViscWallInteraction<T>::calcForces()
{
  Vec3 vp=this->m_p->getVel();
  Vec3 vw=this->m_wall->getVel();

  Vec3 force=m_nu*this->m_p->getMass()*(vp-vw);
  Vec3 pos=this->m_p->getPos();

  this->m_p->applyForce(force,pos);
  if(this->m_inner_flag) this->m_wall->addForce(force);
}

/*!
  calculate and return, but don't apply the viscous force
*/
template <class T>
Vec3 CViscWallInteraction<T>::getForce()
{
  Vec3 vp=this->m_p->getVel();
  Vec3 vw=this->m_wall->getVel();

  Vec3 force=m_nu*this->m_p->getMass()*(vw-vp);

  return force;
}

#endif
