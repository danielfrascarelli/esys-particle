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
#include "ECorner2DInteraction.h"

/*!
  default constructor
*/
ECorner2DInteraction::ECorner2DInteraction()
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
ECorner2DInteraction::ECorner2DInteraction(CParticle* p,Corner2D* c,ETriMeshIP param,bool iflag)
{
  m_p=p;
  m_corner=c;
  m_k=param.k;
  m_inner_flag=iflag;
}

/*!
  destructor
*/
ECorner2DInteraction::~ECorner2DInteraction()
{}

/*!
  calculate & apply forces
*/
void ECorner2DInteraction::calcForces()
{
  Vec3 ppos=m_p->getPos();
  if(m_corner->isValidContact(ppos)){ // if no contact to adjacent edges or triangles
    double sep=m_corner->sep(ppos);
    if(sep<m_p->getRad()){
      Vec3 force=m_k*(m_p->getRad()-sep)*m_corner->getDirectionFromPoint(ppos);
      m_p->applyForce(force,ppos);
    }
  }
}
