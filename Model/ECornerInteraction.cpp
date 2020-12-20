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

#include "ECornerInteraction.h"

/*!
  default constructor
*/
ECornerInteraction::ECornerInteraction()
{
  m_p=NULL;
  m_corner=NULL;
  m_k=0.0;
  m_inner_flag=false;
}

/*!
  constructor with parameters

  \param p a pointer to the particle
  \param c a pointer to the corner
  \param param the interaction parameters
  \param iflag
*/
ECornerInteraction::ECornerInteraction(CParticle* p,Corner* c,ETriMeshIP param,bool iflag)
{
  m_p=p;
  m_corner=c;
  // scale elastic param
  double f=1.0;
  if(!CParticle::getDo2dCalculations()){
    f*=3.141592654*this->m_p->getRad();
  }
  m_k=f*param.k;
  m_inner_flag=iflag;
}

/*!
  destructor
*/
ECornerInteraction::~ECornerInteraction()
{}

/*!
  calculate & apply forces
*/
void ECornerInteraction::calcForces()
{
  Vec3 ppos=m_p->getPos();
  if(m_corner->isValidContact(ppos)){ // if no contact to adjacent edges or triangles
    double sep=m_corner->sep(ppos);
    if(sep<m_p->getRad()){
      Vec3 force=m_k*(m_p->getRad()-sep)*m_corner->getDirectionFromPoint(ppos);
      m_p->applyForce(force,ppos);
      if(m_inner_flag) m_corner->applyForce(-1.0*force);
    }
  }
}
