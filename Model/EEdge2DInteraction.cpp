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
#include "EEdge2DInteraction.h"
#include "Foundation/vec3.h"

// -- STL includes --
#include <utility>
using std::pair;

/*!
  default constructor
*/
EEdge2DInteraction::EEdge2DInteraction()
{
  m_p=NULL;
  m_edge=NULL;
  m_k=0.0;
  m_inner_flag=false;
}

/*!
  constructor with parameters

  \param p
  \param e
  \param param
  \param iflag
*/
EEdge2DInteraction::EEdge2DInteraction(CParticle* p,Edge2D* e,ETriMeshIP param ,bool iflag)
{
  m_p=p;
  m_edge=e;
  m_k=param.k;
  m_inner_flag=iflag;
}

/*!
  destructor
*/
EEdge2DInteraction::~EEdge2DInteraction()
{}


/*!
  calculate & apply forces
*/
void EEdge2DInteraction::calcForces()
{
  Vec3 ppos=m_p->getPos();
  pair<bool,double> dist=m_edge->dist(ppos);
  if(dist.first && (dist.second<m_p->getRad())){
    Vec3 force=m_k*(m_p->getRad()-dist.second)*m_edge->getNormal();
    Vec3 pos=m_p->getPos()-dist.second*m_edge->getNormal();
    m_p->applyForce(force,pos);
  }
}
