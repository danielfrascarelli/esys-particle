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

#include "ETriangleInteraction.h"

// -- STL includes --
#include <utility>
using std::pair;

/*!
  default constructor
*/
ETriangleInteraction::ETriangleInteraction()
{
  m_p=NULL;
  m_t=NULL;
  m_k=0.0;
  m_inner_flag=false;
}

/*!
  constructor with parameters

  \param p
  \param t
  \param param
  \param iflag
*/
ETriangleInteraction::ETriangleInteraction(CParticle* p,Triangle* t,ETriMeshIP param,bool iflag)
{
  m_p=p;
  m_t=t;
  // scale elastic param
  double f=1.0;
  if(!CParticle::getDo2dCalculations()){
    f*=3.141592654*this->m_p->getRad();
  }
  m_k=f*param.k;
  m_inner_flag=iflag;

  m_normal = m_t->getNormal();

  //check if triangle normal is in direction of particle; swap if req'd:
  Vec3 p0 = (m_t->getP0()).second;
  Vec3 dp = (m_p->getPos() - p0);
  if (dot(dp,m_normal) < 0.0) {
    m_normal *= -1.0;
  }
}

/*!
  destructor
*/
ETriangleInteraction::~ETriangleInteraction()
{}

/*!
  calculate & apply forces
*/
void ETriangleInteraction::calcForces()
{
  pair<bool,double> dist=m_t->dist(m_p->getPos());
  double pdist = dist.second*(m_t->getNormal()*m_normal);
  //if(dist.first && (dist.second<m_p->getRad())){
  if(dist.first && (pdist<m_p->getRad())){
    /*
    Vec3 force=m_k*(m_p->getRad()-dist.second)*m_t->getNormal();
    Vec3 pos=m_p->getPos()-dist.second*m_t->getNormal();
    */
    Vec3 force=m_k*(m_p->getRad()-pdist)*m_normal;
    Vec3 pos=m_p->getPos()-pdist*m_normal;
    m_p->applyForce(force,pos);
    if(m_inner_flag) m_t->applyForce(-1.0*force);
  }
}
